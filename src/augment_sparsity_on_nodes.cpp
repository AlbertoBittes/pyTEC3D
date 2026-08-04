/**
 * @file augment_sparsity_on_nodes.cpp
 * @brief AugmentSparsityOnNodes: CP-node pair detection and ghosting-functor implementation.
 */
#include "augment_sparsity_on_nodes.h"

#include <algorithm>
#include <cstring>
#include <unordered_set>
#include <utility>

#include <libmesh/elem.h>

#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace
{
    bool add_pair(std::vector<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>> &pairs,
                  libMesh::dof_id_type first,
                  libMesh::dof_id_type second)
    {
        const auto candidate = std::make_pair(first, second);
        if (std::find(pairs.begin(), pairs.end(), candidate) != pairs.end())
            return false;
        pairs.push_back(candidate);
        return true;
    }

} // namespace

AugmentSparsityOnNodes::AugmentSparsityOnNodes(
    libMesh::MeshBase &mesh,
    vtkUnstructuredGrid *vtk_grid)
    : _mesh(mesh), _vtk_grid(vtk_grid)
{
    rebuild_elem_neighbors();
}

void AugmentSparsityOnNodes::rebuild_node_pairs()
{
    cp_pairs_by_cp.clear();
    node_pairs.clear();

    if (!_vtk_grid)
        return;

    vtkPointData *point_data = _vtk_grid->GetPointData();
    if (!point_data)
        return;

    std::vector<vtkDataArray *> cp_arrays;
    cp_arrays.reserve(static_cast<std::size_t>(point_data->GetNumberOfArrays()));

    for (int i = 0; i < point_data->GetNumberOfArrays(); ++i)
    {
        vtkDataArray *array = vtkDataArray::SafeDownCast(point_data->GetArray(i));
        if (!array)
            continue;

        const char *name = array->GetName();
        if (!name || std::strncmp(name, "CP_", 3) != 0)
            continue;

        if (array->GetNumberOfComponents() != 1)
            continue;

        cp_arrays.push_back(array);
    }

    if (cp_arrays.empty())
        return;

    std::unordered_map<int, std::unordered_set<libMesh::dof_id_type>> cp_nodes;
    const vtkIdType n_points = _vtk_grid->GetNumberOfPoints();

    // Match legacy behavior: only build CP pairs from nodes known on this mesh
    // partition (local + currently ghosted), not from all global VTK points.
    for (const auto *node_ptr : _mesh.node_ptr_range())
    {
        if (!node_ptr)
            continue;

        const libMesh::dof_id_type node_id = static_cast<libMesh::dof_id_type>(node_ptr->id());
        const vtkIdType vtk_node_id = static_cast<vtkIdType>(node_id);

        for (vtkDataArray *cp_array : cp_arrays)
        {
            const int cp_value = static_cast<int>(cp_array->GetComponent(vtk_node_id, 0));
            if (cp_value <= 0)
                continue;

            cp_nodes[cp_value].insert(node_id);
        }
    }

    for (const auto &entry : cp_nodes)
    {
        const int cp_id = entry.first;
        const auto &nodes = entry.second;
        if (nodes.empty())
            continue;

        const libMesh::dof_id_type master_node = *std::min_element(nodes.begin(), nodes.end());
        auto &cp_pairs = cp_pairs_by_cp[cp_id];
        for (libMesh::dof_id_type node_id : nodes)
        {
            if (node_id != master_node)
            {
                cp_pairs.emplace_back(master_node, node_id);
                node_pairs.emplace_back(master_node, node_id);
            }
        }

        std::sort(cp_pairs.begin(), cp_pairs.end());
        cp_pairs.erase(std::unique(cp_pairs.begin(), cp_pairs.end()), cp_pairs.end());
    }

    std::sort(node_pairs.begin(), node_pairs.end());
    node_pairs.erase(std::unique(node_pairs.begin(), node_pairs.end()), node_pairs.end());
}

void AugmentSparsityOnNodes::rebuild_elem_neighbors()
{
    rebuild_node_pairs();

    elem_pairs.clear();
    slave_to_master_elems.clear();

    std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::dof_id_type>> node_to_elems;
    for (const auto *elem : _mesh.element_ptr_range())
    {
        if (!elem)
            continue;
        for (unsigned int i = 0; i < elem->n_nodes(); ++i)
            node_to_elems[static_cast<libMesh::dof_id_type>(elem->node_id(i))].push_back(elem->id());
    }

    for (auto &entry : node_to_elems)
    {
        auto &elem_ids = entry.second;
        std::sort(elem_ids.begin(), elem_ids.end());
        elem_ids.erase(std::unique(elem_ids.begin(), elem_ids.end()), elem_ids.end());
    }

    for (const auto &pair : node_pairs)
    {
        const auto master_elem_it = node_to_elems.find(pair.first);
        const auto slave_elem_it = node_to_elems.find(pair.second);
        if (master_elem_it == node_to_elems.end() || slave_elem_it == node_to_elems.end())
            continue;

        const auto &master_elem_ids = master_elem_it->second;
        const auto &slave_elem_ids = slave_elem_it->second;
        if (master_elem_ids.empty() || slave_elem_ids.empty())
            continue;

        for (const libMesh::dof_id_type master_elem_id : master_elem_ids)
        {
            for (const libMesh::dof_id_type slave_elem_id : slave_elem_ids)
            {
                add_pair(elem_pairs, master_elem_id, slave_elem_id);

                auto &master_ids = slave_to_master_elems[slave_elem_id];
                if (std::find(master_ids.begin(), master_ids.end(), master_elem_id) == master_ids.end())
                    master_ids.push_back(master_elem_id);
            }
        }
    }

    for (auto &entry : slave_to_master_elems)
    {
        auto &master_ids = entry.second;
        std::sort(master_ids.begin(), master_ids.end());
        master_ids.erase(std::unique(master_ids.begin(), master_ids.end()), master_ids.end());
    }
}

void AugmentSparsityOnNodes::operator()(const libMesh::MeshBase::const_element_iterator &range_begin,
                                        const libMesh::MeshBase::const_element_iterator &range_end,
                                        libMesh::processor_id_type p,
                                        libMesh::GhostingFunctor::map_type &coupled_elements)
{
    const libMesh::CouplingMatrix *const null_mat = nullptr;

    for (const auto &elem : libMesh::as_range(range_begin, range_end))
    {
        if (!elem)
            continue;

        if (elem->processor_id() != p)
            continue;

        const auto master_it = slave_to_master_elems.find(elem->id());
        if (master_it == slave_to_master_elems.end())
            continue;

        for (const libMesh::dof_id_type master_elem_id : master_it->second)
        {
            const libMesh::Elem *master_elem = _mesh.elem_ptr(master_elem_id);
            if (master_elem && master_elem->processor_id() != p)
                coupled_elements.emplace(master_elem, null_mat);
        }
    }
}

void AugmentSparsityOnNodes::mesh_reinit()
{
    rebuild_elem_neighbors();
}

void AugmentSparsityOnNodes::redistribute()
{
    this->mesh_reinit();
}

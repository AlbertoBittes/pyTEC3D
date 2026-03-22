#include "augment_sparsity_on_nodes.h"

#include <algorithm>
#include <cstring>
#include <unordered_set>
#include <utility>

#include <libmesh/elem.h>

#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

using namespace libMesh;

namespace
{
    bool add_pair(std::vector<std::pair<dof_id_type, dof_id_type>> &pairs,
                  dof_id_type first,
                  dof_id_type second)
    {
        const auto candidate = std::make_pair(first, second);
        if (std::find(pairs.begin(), pairs.end(), candidate) != pairs.end())
            return false;
        pairs.push_back(candidate);
        return true;
    }

} // namespace

AugmentSparsityOnNodes::AugmentSparsityOnNodes(
    MeshBase &mesh,
    vtkUnstructuredGrid *vtk_grid)
    : _mesh(mesh), _vtk_grid(vtk_grid)
{
    rebuild_elem_neighbors();
}

void AugmentSparsityOnNodes::rebuild_node_pairs()
{
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

    vtkDataArray *node_id_array =
        vtkDataArray::SafeDownCast(point_data->GetAbstractArray("libmesh_node_id"));
    const vtkIdType n_points = _vtk_grid->GetNumberOfPoints();

    std::unordered_map<int, std::unordered_set<dof_id_type>> cp_nodes;

    for (vtkDataArray *cp_array : cp_arrays)
        for (vtkIdType i = 0; i < n_points; ++i)
        {
            const int cp_value = static_cast<int>(cp_array->GetComponent(i, 0));
            if (cp_value <= 0)
                continue;

            const dof_id_type node_id =
                node_id_array
                    ? static_cast<dof_id_type>(node_id_array->GetComponent(i, 0))
                    : static_cast<dof_id_type>(i);

            cp_nodes[cp_value].insert(node_id);
        }

    for (const auto &entry : cp_nodes)
    {
        const auto &nodes = entry.second;
        if (nodes.empty())
            continue;

        const dof_id_type master_node = *std::min_element(nodes.begin(), nodes.end());
        for (dof_id_type node_id : nodes)
            if (node_id != master_node)
                node_pairs.emplace_back(master_node, node_id);
    }

    std::sort(node_pairs.begin(), node_pairs.end());
    node_pairs.erase(std::unique(node_pairs.begin(), node_pairs.end()), node_pairs.end());
}

void AugmentSparsityOnNodes::rebuild_elem_neighbors()
{
    rebuild_node_pairs();

    elem_pairs.clear();
    slave_to_master_elem.clear();

    std::unordered_map<dof_id_type, dof_id_type> node_to_elem;
    for (const auto *elem : _mesh.element_ptr_range())
    {
        if (!elem)
            continue;
        for (unsigned int i = 0; i < elem->n_nodes(); ++i)
            node_to_elem.emplace(static_cast<dof_id_type>(elem->node_id(i)), elem->id());
    }

    for (const auto &pair : node_pairs)
    {
        const auto master_elem_it = node_to_elem.find(pair.first);
        const auto slave_elem_it = node_to_elem.find(pair.second);
        if (master_elem_it == node_to_elem.end() || slave_elem_it == node_to_elem.end())
            continue;

        if (!add_pair(elem_pairs, master_elem_it->second, slave_elem_it->second))
            continue;

        slave_to_master_elem.emplace(slave_elem_it->second, master_elem_it->second);
    }
}

void AugmentSparsityOnNodes::operator()(const MeshBase::const_element_iterator &range_begin,
                                        const MeshBase::const_element_iterator &range_end,
                                        processor_id_type p,
                                        map_type &coupled_elements)
{
    const CouplingMatrix *const null_mat = nullptr;

    for (const auto &elem : as_range(range_begin, range_end))
    {
        if (!elem)
            continue;

        if (elem->processor_id() != p)
            continue;

        const auto master_it = slave_to_master_elem.find(elem->id());
        if (master_it == slave_to_master_elem.end())
            continue;

        const Elem *master_elem = _mesh.elem_ptr(master_it->second);
        if (master_elem && master_elem->processor_id() != p)
            coupled_elements.emplace(master_elem, null_mat);
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

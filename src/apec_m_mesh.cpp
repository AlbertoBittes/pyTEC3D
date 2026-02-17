#include "apec_m_mesh.h"

#include "vtkio_updated.h"
#include "element_types_factory.h"

#include <libmesh/elem.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/node.h>
#include <libmesh/point.h>
#include <libmesh/vector_value.h>

#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include "fem_data_repository.h"

std::shared_ptr<ApecMMesh> ApecMMesh::read_mesh(libMesh::Parallel::Communicator &comm,
                                                FEMDataRepository &mesh_repository,
                                                const std::string &mesh_filepath)
{
    libMesh::out << "Reading mesh and importing field data..." << std::endl;

    if (mesh_filepath.empty())
        libmesh_error_msg("Mesh filepath not provided for ApecMMesh::read_mesh.");

    auto mesh = mesh_repository.get_object<ApecMMesh>();
    if (!mesh)
    {
        mesh = std::make_shared<ApecMMesh>(comm);
        mesh_repository.set_object<ApecMMesh>(mesh);
    }
    VTKIO_Updated vtkio(*mesh, mesh_filepath);
    MeshCommunication().broadcast(*mesh);
    MeshTools::Modification::orient_elements(*mesh);

    return mesh;
}

std::vector<std::string> ApecMMesh::collect_node_integer_arrays(const std::string &prefix) const
{
    std::vector<std::string> arrays;
    for (const auto &name : VTKIO_Updated::VTK_DATA_ARRAY_NAMES)
    {
        if (name.rfind(prefix, 0) != 0)
            continue;

        if (!name.empty() && name.back() == '_')
        {
            for (unsigned int idx = 0;; ++idx)
            {
                std::string generated = name + std::to_string(idx);
                if (!this->has_node_integer(generated))
                    break;
                arrays.push_back(generated);
            }
        }
        else if (this->has_node_integer(name))
        {
            arrays.push_back(name);
        }
    }
    return arrays;
}

bool ApecMMesh::generate_nodal_couplings()
{
    std::vector<unsigned int> cp_integer_indices;
    std::vector<std::string> cp_arrays = collect_node_integer_arrays("CP_");
    for (const auto &name : cp_arrays)
        cp_integer_indices.push_back(this->get_node_integer_index(name));

    if (cp_integer_indices.empty())
    {
        libMesh::out << "No CP_* sets found. Skipping nodal coupling generation." << std::endl;
        return false;
    }

    std::unordered_map<int, std::unordered_set<unsigned int>> cp_nodes;

    cp_pairs.clear();

    for (const auto &node_ptr : this->node_ptr_range())
    {
        if (!node_ptr)
            continue;

        for (unsigned int data_index : cp_integer_indices)
        {
            const int cp_value = static_cast<int>(node_ptr->get_extra_integer(data_index));
            if (cp_value <= 0)
                continue;

            cp_nodes[cp_value].insert(static_cast<unsigned int>(node_ptr->id()));
        }
    }

    for (const auto &entry : cp_nodes)
    {
        const int cp_id = entry.first;
        const auto &nodes = entry.second;
        if (nodes.empty())
            continue;

        unsigned int master_node = *std::min_element(nodes.begin(), nodes.end());
        auto &pairs = cp_pairs[cp_id];
        pairs.clear();

        for (unsigned int node_id : nodes)
        {
            if (node_id == master_node)
                continue;
            pairs.emplace_back(master_node, node_id);
        }
    }

    return true;
}

const std::unordered_map<int, std::vector<std::pair<unsigned int, unsigned int>>> &ApecMMesh::get_cp_pairs() const
{
    return cp_pairs;
}

const std::map<dof_id_type, Real> &ApecMMesh::edge_midside_orientation_map()
{
    if (!edge_midside_orientation_cache_built)
        build_edge_midside_orientation_map();

    return edge_midside_orientation_cache;
}

void ApecMMesh::build_edge_midside_orientation_map()
{
    edge_midside_orientation_cache.clear();

    for (const auto &elem_ptr : this->active_local_element_ptr_range())
    {
        if (!elem_ptr)
            continue;

        const auto n_edges = elem_ptr->n_edges();
        for (unsigned int edge = 0; edge < n_edges; ++edge)
        {
            const Real orientation = elem_ptr->positive_edge_orientation(edge) ? 1.0 : -1.0;
            const auto n_nodes = elem_ptr->n_nodes();
            for (unsigned int node_i = 0; node_i < n_nodes; ++node_i)
            {
                if (elem_ptr->is_vertex(node_i))
                    continue;
                if (!elem_ptr->is_node_on_edge(node_i, edge))
                    continue;

                const Node *node_ptr = elem_ptr->node_ptr(node_i);
                if (!node_ptr)
                    continue;

                const dof_id_type node_id = node_ptr->id();
                edge_midside_orientation_cache.emplace(node_id, orientation);
            }
        }
    }

    edge_midside_orientation_cache_built = true;
}

ApecMMesh::AngleDataIndices ApecMMesh::detect_angle_indices() const
{
    AngleDataIndices info;
    bool all_found = true;
    for (unsigned int i = 0; i < 3; ++i)
    {
        const std::string name = "angles-" + std::to_string(i);
        if (this->has_node_integer(name))
            info.indices(i) = this->get_node_integer_index(name);
        else
        {
            all_found = false;
            break;
        }
    }
    if (all_found)
        info.available = true;
    return info;
}

CartesianCoordinateSystem ApecMMesh::node_coordinate_system(dof_id_type node_id) const
{
    const Node *node_ptr = this->query_node_ptr(node_id);
    const AngleDataIndices angle_data = detect_angle_indices();
    if (node_ptr && angle_data.available)
    {
        VectorValue<Real> angles_deg(0., 0., 0.);
        for (unsigned int i = 0; i < 3; ++i)
            angles_deg(i) = node_ptr->get_extra_datum<Real>(angle_data.indices(i));
        return CartesianCoordinateSystem(Point(0.0, 0.0, 0.0),
                                         angles_deg(0),
                                         angles_deg(1),
                                         angles_deg(2));
    }

    return CartesianCoordinateSystem();
}

CartesianCoordinateSystem ApecMMesh::node_coordinate_system(const Node &node) const
{
    const AngleDataIndices angle_data = detect_angle_indices();
    if (!angle_data.available)
        return CartesianCoordinateSystem();

    VectorValue<Real> angles_deg(0., 0., 0.);
    for (unsigned int i = 0; i < 3; ++i)
        angles_deg(i) = node.get_extra_datum<Real>(angle_data.indices(i));

    return CartesianCoordinateSystem(Point(0.0, 0.0, 0.0),
                                     angles_deg(0),
                                     angles_deg(1),
                                     angles_deg(2));
}

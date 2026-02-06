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

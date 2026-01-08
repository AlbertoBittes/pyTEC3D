#include "apec_m_mesh.h"

#include "inf_elem_builder_updated.h"
#include "vtkio_updated.h"
#include "element_types_factory.h"

#include <libmesh/elem.h>
#include <libmesh/boundary_info.h>
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

void ApecMMesh::apply_boundary_info()
{
    libMesh::out << "Assigning boundary information..." << std::endl;

    BoundaryInfo &boundary_info = this->get_boundary_info();
    boundary_info.clear();

    std::vector<std::string> boundary_arrays = collect_node_integer_arrays("boundary_id_");

    for (const auto &array : boundary_arrays)
    {
        unsigned int node_array_index = this->get_node_integer_index(array);

        for (const Node *node_ptr : this->node_ptr_range())
        {
            unsigned int boundary_id_mesh_int = node_ptr->get_extra_integer(node_array_index);
            boundary_id_type boundary_id_mesh = static_cast<boundary_id_type>(boundary_id_mesh_int);
            if (boundary_id_mesh == 0)
                continue;
            boundary_info.add_node(node_ptr, boundary_id_mesh);
        }

        boundary_info.build_side_list_from_node_list();
    }

    boundary_info.print_info();
    boundary_info.print_summary();
}

void ApecMMesh::generate_infinite_elements(const std::shared_ptr<ElementTypesFactory> &element_types_factory)
{
    libMesh::out << "Generating infinite elements..." << std::endl;

    if (!element_types_factory)
    {
        libMesh::out << "No ElementTypesFactory attached; skipping infinite element generation." << std::endl;
        return;
    }

    this->apply_boundary_info();

    BoundaryInfo &boundary_info = this->get_boundary_info();
    // boundary_info.clear();

    const std::set<libMesh::boundary_id_type> &side_boundary_ids = boundary_info.get_side_boundary_ids();
    if (side_boundary_ids.empty())
        return;

    ElementTypesMap *types_map = element_types_factory->get_element_types_map();
    for (const auto &id : side_boundary_ids)
    {
        auto it = types_map->find(id);
        if (it == types_map->end())
            continue;

        ElementTypeEnum elem_enum = static_cast<ElementTypeEnum>(it->second->get_keyopts().front());
        if (elem_enum != ElementTypeEnum::INFIN110 &&
            elem_enum != ElementTypeEnum::INFIN111 &&
            elem_enum != ElementTypeEnum::INFIN112)
            continue;

        InfElemBuilderUpdated builder(*this);
        builder.set_infinite_elem_side_boundary_id(id);
        builder.build_inf_elem(Point(0.0, 0.0, 0.0));

        for (auto &elem : this->element_ptr_range())
            if (elem->infinite())
                elem->subdomain_id() = id;
    }
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

    std::unordered_set<unsigned int> nodeelem_nodes;
    std::size_t created = 0;

    auto ensure_nodeelem = [&](unsigned int node_id)
    {
        if (!nodeelem_nodes.insert(node_id).second)
            return;

        Node *node = this->query_node_ptr(node_id);
        if (!node)
        {
            libMesh::out << "Warning: node " << node_id
                         << " not found while creating nodal couplings." << std::endl;
            nodeelem_nodes.erase(node_id);
            return;
        }

        auto node_elem = this->add_elem(Elem::build(NODEELEM));
        node_elem->set_node(0, node);
        ++created;
    };

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
            ensure_nodeelem(node_id);
        }

        ensure_nodeelem(master_node);
    }

    libMesh::out << "Generated " << created << " nodeelems." << std::endl;
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

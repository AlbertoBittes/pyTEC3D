#include "apec_m_mesh.h"

#include "globals.h"
#include "inf_elem_builder_updated.h"
#include "vtkio_updated.h"
#include "element_types_factory.h"

#include <libmesh/boundary_info.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/node.h>
#include <libmesh/point.h>

#include "fem_data_repository.h"

std::shared_ptr<ApecMMesh> ApecMMesh::read(libMesh::Parallel::Communicator &comm,
                                           FEMDataRepository &model_data,
                                           const std::string &mesh_filepath)
{
    libMesh::out << "Reading mesh and importing field data..." << std::endl;

    if (mesh_filepath.empty())
        libmesh_error_msg("Mesh filepath not provided for ApecMMesh::read.");

    auto mesh = model_data.get_object<ApecMMesh>();
    if (!mesh)
    {
        mesh = std::make_shared<ApecMMesh>(comm);
        model_data.set_object<ApecMMesh>(mesh);
    }

    VTKIO_Updated vtkio(*mesh, mesh_filepath);
    MeshCommunication().broadcast(*mesh);
    MeshTools::Modification::orient_elements(*mesh);

    if (auto source_repo = vtkio.get_field_data_repo())
    {
        model_data.move_all_from(source_repo);
    }

    return mesh;
}

void ApecMMesh::apply_boundary_info()
{
    libMesh::out << "Assigning boundary information..." << std::endl;

    BoundaryInfo &boundary_info = this->get_boundary_info();
    boundary_info.clear();

    std::vector<std::string> boundary_arrays;
    for (const auto &name : Globals::VTK_DATA_ARRAY_NAMES)
    {
        if (name.find("boundary_id_") == 0)
        {
            for (unsigned int idx = 0;; ++idx)
            {
                std::string generated = name + std::to_string(idx);
                if (!this->has_node_integer(generated))
                    break;
                boundary_arrays.push_back(generated);
            }
        }
    }

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
    boundary_info.clear();

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

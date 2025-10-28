#include "apec_m_mesh.h"

#include "globals.h"
#include "inf_elem_builder_updated.h"
#include "vtkio_updated.h"

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

void ApecMMesh::apply_boundary_info(const std::set<std::string> &bc_arrays)
{
    libMesh::out << "Assigning boundary information..." << std::endl;

    BoundaryInfo &boundary_info = this->get_boundary_info();
    boundary_info.clear();

    for (const auto &array : bc_arrays)
    {
        unsigned int node_array_index = this->get_node_integer_index(array);
        Globals::BoundaryConditionType bc_enum = Globals::get_key_by_predicate_one(
            Globals::BOUNDARY_CONDITION_INFO_MAP,
            [array](const Globals::BoundaryConditionInfo &info)
            {
                return info.label == array;
            });

        for (const Node *node_ptr : this->node_ptr_range())
        {
            unsigned int boundary_id_mesh_int = node_ptr->get_extra_integer(node_array_index);
            boundary_id_type boundary_id_mesh = static_cast<boundary_id_type>(boundary_id_mesh_int);
            if (boundary_id_mesh == 0)
                continue;

            switch (bc_enum)
            {
            case Globals::BoundaryConditionType::INFINITE_BOUNDARY:
            case Globals::BoundaryConditionType::PERIODIC_BOUNDARY:
                boundary_info.add_node(node_ptr, boundary_id_mesh);
                break;
            case Globals::BoundaryConditionType::FLUX_NORMAL:
                boundary_info.add_node(node_ptr, Globals::FLUX_NORMAL_ID_BC);
                break;
            case Globals::BoundaryConditionType::FLUX_PARALLEL:
                boundary_info.add_node(node_ptr, Globals::FLUX_PARALLEL_ID_BC);
                break;
            default:
                break;
            }
        }

        boundary_info.build_side_list_from_node_list();
    }

    boundary_info.print_info();
    boundary_info.print_summary();
}

void ApecMMesh::generate_infinite_elements()
{
    libMesh::out << "Generating infinite elements..." << std::endl;

    this->apply_boundary_info({"infinite_boundary"});

    if (!this->has_node_integer("infinite_boundary"))
        return;

    BoundaryInfo &boundary_info = this->get_boundary_info();
    const std::set<libMesh::boundary_id_type> &side_boundary_ids = boundary_info.get_side_boundary_ids();

    unsigned int etype_num_index = this->get_elem_integer_index("ansys_etype");
    unsigned int material_num_index = this->get_elem_integer_index("ansys_material_type");
    unsigned int real_constants_num_index = this->get_elem_integer_index("ansys_real_constant");

    for (const auto &id : side_boundary_ids)
    {
        InfElemBuilderUpdated builder(*this);
        builder.set_infinite_elem_side_boundary_id(id);
        const libMesh::Point origin(0.0, 0.0, 0.0);
        builder.build_inf_elem(origin);

        for (auto &elem : this->element_ptr_range())
            if (elem->infinite())
            {
                elem->subdomain_id() = id;
                elem->set_extra_integer(etype_num_index, id);
                elem->set_extra_integer(material_num_index, 0);
                elem->set_extra_integer(real_constants_num_index, 0);
            }
    }
}

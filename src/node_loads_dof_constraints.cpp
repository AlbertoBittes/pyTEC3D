#include "node_loads_dof_constraints.h"
#include "apec_m_system_objects.h"
#include <libmesh/numeric_vector.h>
#include <libmesh/node.h>
#include "globals.h"

using namespace libMesh;

struct SystemDataRepo
{
    unsigned int sys_number;
    NumericVector<Number> *vec;
    NumericVector<Number> *old_vec;
    Globals::SystemInfo *system_info; // pointer to system_info attribute
};

void NodeLoadsDOFConstraints::read_data(std::ifstream &file)
{
    libMesh::out << "Reading node loads data from an open file stream..." << std::endl;

    // Build system contexts for all relevant systems
    std::vector<SystemDataRepo> systems;

    // DOF system
    {
        ApecMLinearSystem *system = &data_map->get_system<ApecMLinearSystem>("ApecMLinearSystem");
        systems.push_back({system->number(),
                           system->request_vector("node_dof_constraints"),
                           system->request_vector("old_node_dof_constraints"),
                           &system->system_info});
    }

    // Nodal loads system
    {
        TransientSystem<ExplicitSystem> *nodal_loads_system = &data_map->get_system<TransientSystem<ExplicitSystem>>("NodalLoadsSystem");
        // Cast to ApecMSystem to access system_info
        auto *apecm_loads = dynamic_cast<ApecMLoadsSystem *>(nodal_loads_system);
        systems.push_back({nodal_loads_system->number(),
                           nodal_loads_system->current_local_solution.get(),
                           nodal_loads_system->old_local_solution,
                           &apecm_loads->system_info});
    }

    // Initialize all vectors: copy current to old, set current to NaN
    for (auto &ctx : systems)
    {
        *ctx.old_vec = *ctx.vec;
        *ctx.vec = std::numeric_limits<Real>::quiet_NaN();
    }

    dof_id_type num_cell_data = 0;
    file.read(reinterpret_cast<char *>(&num_cell_data), sizeof(dof_id_type));
    if (!file)
    {
        libmesh_error_msg("Error: Failed to read the number of loads from the binary stream.");
    }

    const MeshBase &mesh = data_map->get_mesh();

    dof_id_type prev_node_load_label_number = static_cast<dof_id_type>(-1);
    dof_id_type prev_component_number = static_cast<dof_id_type>(-1);
    Globals::VariableCategory cached_category;
    std::string cached_variable_label;
    unsigned int cached_var_num = 0;
    unsigned int cached_sys_number = 0;

    for (dof_id_type i = 0; i < num_cell_data; ++i)
    {
        dof_id_type node_number = 0;
        dof_id_type node_load_label_number = 0;
        dof_id_type component_number = 0;
        Real load_value = 0.0;

        file.read(reinterpret_cast<char *>(&node_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&node_load_label_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&component_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&load_value), sizeof(Real));

        const Node *node_ptr = mesh.query_node_ptr(static_cast<unsigned int>(node_number));
        if (!node_ptr)
            continue;

        if (node_load_label_number != prev_node_load_label_number || component_number != prev_component_number)
        {
            // Map the load label number to the corresponding API enum
            Globals::APIVariableName api_var_enum = static_cast<Globals::APIVariableName>(static_cast<unsigned int>(node_load_label_number));

            Globals::GlobalVariableName load_label_global = Globals::get_key_by_predicate_one(
                Globals::GLOBAL_VARIABLE_INFO_MAP,
                [api_var_enum, component_number](const Globals::VariableInfo &info)
                {
                    return info.api_name == api_var_enum && info.component_index == static_cast<unsigned int>(component_number);
                });

            // Retrieve the variable category from the VariableInfo
            cached_category = Globals::get_info_by_key(
                                  Globals::GLOBAL_VARIABLE_INFO_MAP, load_label_global)
                                  .category;

            // Retrieve the variable name in string
            cached_variable_label = Globals::get_info_by_key(
                                        Globals::GLOBAL_VARIABLE_INFO_MAP, load_label_global)
                                        .global_label;

            // Cache var_num and sys_number for the current system
            for (auto &ctx : systems)
            {
                if (cached_category == ctx.system_info->category)
                {
                    cached_var_num = data_map->get_system(ctx.sys_number).variable_number(cached_variable_label);
                    cached_sys_number = ctx.sys_number;
                    break;
                }
            }

            prev_node_load_label_number = node_load_label_number;
            prev_component_number = component_number;
        }

        for (auto &ctx : systems)
        {
            if (cached_category == ctx.system_info->category)
            {
                // Use cached_var_num and cached_sys_number
                dof_id_type dof_id = node_ptr->dof_number(cached_sys_number, cached_var_num, 0);
                ctx.vec->set(dof_id, load_value);
                break;
            }
        }
    }
}

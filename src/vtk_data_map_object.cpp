#include "globals.h"
#include "variables_map.h"
#include <libmesh/numeric_vector.h>
#include <libmesh/dof_map.h>

#include "vtk_data_map_object.h"

void VTKDataMapFEDataObject::read_data(std::ifstream &file)
{
    libMesh::out << "Reading VTK data map from binary file..." << std::endl;

    // Initialize all vectors: copy current to old, set current to NaN
    for (auto &ctx : systems)
    {
        if (ctx.old_vec)
            ctx.old_vec->close();
        if (ctx.vec)
            ctx.vec->close();
        if (ctx.old_vec && ctx.vec)
            *ctx.old_vec = *ctx.vec;
        if (ctx.vec)
            *ctx.vec = ctx.initial_val; // Set current vector to NaN or initial value
    }

    dof_id_type num_data = 0;
    file.read(reinterpret_cast<char *>(&num_data), sizeof(dof_id_type));
    if (!file)
    {
        libmesh_error_msg("Error: Failed to read the number of loads from the binary stream.");
    }

    const MeshBase &mesh = data_map->get_mesh();

    dof_id_type prev_label_number = static_cast<dof_id_type>(-1);
    dof_id_type prev_component_number = static_cast<dof_id_type>(-1);
    Globals::SystemType cached_sys_type;
    std::string cached_variable_label;
    unsigned int cached_var_num = 0;
    unsigned int cached_sys_number = 0;
    std::vector<dof_id_type> dof_indices;

    for (dof_id_type i = 0; i < num_data; ++i)
    {
        dof_id_type entity_number = 0;
        dof_id_type label_number = 0;
        dof_id_type component_number = 0;
        Real value = 0.0;

        file.read(reinterpret_cast<char *>(&entity_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&label_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&component_number), sizeof(dof_id_type));
        file.read(reinterpret_cast<char *>(&value), sizeof(Real));

        const Node *node_ptr = nullptr;
        const Elem *elem_ptr = nullptr;

        if (is_element_data)
            elem_ptr = mesh.query_elem_ptr(static_cast<unsigned int>(entity_number));
        else
            node_ptr = mesh.query_node_ptr(static_cast<unsigned int>(entity_number));

        // if the node/element is not found, skip (for now)
        if ((is_element_data && !elem_ptr) || (!is_element_data && !node_ptr))
            continue;

        if (label_number != prev_label_number || component_number != prev_component_number)
        {
            APIVariableName api_var_enum = static_cast<APIVariableName>(static_cast<unsigned int>(label_number));

            auto global_label = VariablesMap::get_key_by_predicate_one(
                [api_var_enum, component_number](const VariableInfo &info)
                {
                    return info.api_name == api_var_enum &&
                           info.component_index == static_cast<unsigned int>(component_number);
                });

            auto var_info = VariablesMap::get_info_by_key(global_label);

            cached_sys_type = var_info.system_type;

            cached_variable_label = var_info.global_label;

            for (auto &ctx : systems)
            {
                if (cached_sys_type == *ctx.system_type)
                {
                    cached_var_num = data_map->get_system(ctx.sys_number).variable_number(cached_variable_label);
                    cached_sys_number = ctx.sys_number;
                    break;
                }
            }

            prev_label_number = label_number;
            prev_component_number = component_number;
        }

        for (auto &ctx : systems)
        {
            if (cached_sys_type == *ctx.system_type)
            {
                dof_id_type dof_id = is_element_data ? elem_ptr->dof_number(cached_sys_number, cached_var_num, 0) : node_ptr->dof_number(cached_sys_number, cached_var_num, 0);
                ctx.vec->set(dof_id, value);
                break;
            }
        }
    }
}

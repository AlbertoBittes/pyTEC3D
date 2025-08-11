#include "fe_parent_object.h"
#include "global_constraints_object.h"
#include "node_dof_constraints.h"
#include "globals.h"
#include "fem_data_repository.h"
#include "fe_data_object.h"
#include <libmesh/dof_map.h>
#include <libmesh/node.h>

void GlobalConstraintsObject::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.
    // Retrieve the system by name
    MeshBase &mesh = this->parent_object->get_mesh();
    DofMap &dof_map = this->parent_object->get_dof_map();
    if (!this->ctx_view)
    {
        return;
    }
    auto node_bc_object = this->context<NodeDOFConstraints>();
    if (!node_bc_object)
    {
        return;
    }
    auto data_map = node_bc_object->get_data_map();
    std::string string_variable_label;
    for (const auto &bc_pair : *(data_map))
    {
        // This will not work for edge elements
        // Access the node pointer using the node_id, you will have to replace this with something taht can apply bcs for edge elements
        const Node *node_ptr = mesh.query_node_ptr(bc_pair.first);

        if (node_ptr)
        {
            const auto &bc_values = bc_pair.second;

            for (const auto &bc_value : bc_values)
            {
                string_variable_label = Globals::enum_to_string<Globals::GlobalVariableName>(bc_value.first, Globals::ENUM_TO_STRING);
                // if (string_variable_label == "VOLT" && !this->context<TimeStep>()->get_transient_effect())
                // {
                //     continue; // Skip this variable
                // }
                if (parent_object->has_variable(string_variable_label))
                {
                    // if (mesh.comm().rank() != 0)
                    // {
                    unsigned int variable_number = parent_object->variable_number(string_variable_label);
                    // std::cout << "Applying boundary condition for variable: " << string_variable_label
                    //           << " on node: " << node_ptr->id() << std::endl;
                    // std::cout << "Value: " << bc_value.second << std::endl;
                    // // Apply the constraint for the variable
                    // std::cout << "Variable number: " << variable_number << std::endl;

                    dof_id_type dof_id = node_ptr->dof_number(parent_object->number(), variable_number, 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(dof_id, constraint_row, bc_value.second, true);
                    // }
                }
                else
                    libmesh_error_msg("Variable " + string_variable_label + " not found");
            }
        }
        else
        {
            continue;
        }
    }

    // If we have VOLT but transient effects are not active, bound V to 0 to avoid crash of the solver
    // if (!this->context<TimeStep>()->get_transient_effect() && parent_object->has_variable("VOLT"))
    // {
    //     // Apply the constraint for the time derivative of the magnetic vector potential
    //     std::vector<dof_id_type> dof_indices;
    //     dof_map.local_variable_indices(dof_indices, mesh, parent_object->variable_number("VOLT"));
    //     for (const auto &dof : dof_indices)
    //     {
    //         // Apply the constraint
    //         DofConstraintRow constraint_row;
    //         dof_map.add_constraint_row(dof, constraint_row, true);
    //     }
    // }
}

void GlobalConstraintsObject::check_update_boundary_conditions()
{
    this->update_boundary_conditions = false;
    auto prev_node_bc_object = this->old_context<NodeDOFConstraints>();
    auto curr_node_bc_object = this->context<NodeDOFConstraints>();

    if (!prev_node_bc_object || !curr_node_bc_object)
    {
        this->update_boundary_conditions = true;
        return;
    }

    if (*(prev_node_bc_object->get_data_map()) == *(curr_node_bc_object->get_data_map()))
        this->update_boundary_conditions = false;
    else
        this->update_boundary_conditions = true;
}

void GlobalConstraintsObject::set_boundary_conditions()
{
    this->check_update_boundary_conditions();
    if (this->update_boundary_conditions)
    {
        parent_object->reinit_constraints();
    }
}
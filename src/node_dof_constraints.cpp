#include "libmesh/equation_systems.h"
#include "libmesh/node.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h" // Include the header file for LinearImplicitSystem
#include "../include/node_dof_constraints.h"

using namespace libMesh;
// BCs are probably broken now

void NodalBCConstraint::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.
    // Retrieve the system by name
    MeshBase &mesh = system.get_mesh();
    DofMap &dof_map = system.get_dof_map();

    dof_id_type dof_id;
    const Node *node_ptr = nullptr;
    for (const auto &bc_pair : *(active_bc_map))
    {
        const auto &bc_values = bc_pair.second;

        // This will not work for edge elements
        // Access the node pointer using the node_id, you will have to replace this with something taht can apply bcs for edge elements
        node_ptr = mesh.node_ptr(bc_pair.first);
        if (node_ptr)
        {
            for (const auto &bc_value : bc_values)
            {
                if (system.has_variable(bc_value.first))
                {
                    dof_id_type dof_id = node_ptr->dof_number(system.number(), system.variable_number(bc_value.first), 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(dof_id, constraint_row, bc_value.second, true);
                }
                else
                    libmesh_error_msg("Variable " + bc_value.first + " not found in system " + system.name());
            }
        }
    }

    // If we have VOLT but transient effects are not active, bound V to 0 to avoid crash of the solver
    if (!this->time_step->get_transient_effect() && system.has_variable("V"))
    {
        // Apply the constraint for the time derivative of the magnetic vector potential
        std::vector<dof_id_type> dof_indices;
        dof_map.local_variable_indices(dof_indices, system.get_mesh(), system.variable_number("V"));
        for (const auto &dof : dof_indices)
        {
            // Apply the constraint
            DofConstraintRow constraint_row;
            dof_map.add_constraint_row(dof, constraint_row, true);
        }
    }
}

void NodalBCConstraint::check_boundary_conditions()
{
    // Example implementation
    // This function should implement the logic to check the boundary conditions applied to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be checked.

    // Retrieve the system by name
    bool update_bc = false;

    if (old_bc_map == nullptr)
        update_bc = true;
    else
    {
        if (*old_bc_map == *active_bc_map)
            update_bc = false;
        else
            update_bc = true;
    }
    system.comm().max(update_bc);
    update_boundary_conditions = update_bc;
}
#include "libmesh/equation_systems.h"
#include "libmesh/node.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h" // Include the header file for LinearImplicitSystem
#include "../include/node_dof_constraints.h"
#include <vtkPointData.h>
#include <vtkDataArray.h>

using namespace libMesh;

void NodalBCConstraint::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.

    // Retrieve the system by name
    const MeshBase &mesh = _sys.get_mesh();
    std::string array_name = "A_" + std::to_string(time_stepper.get_current_time_step() + 1);
    vtkDataArray *A_Array = _pointData->GetArray(array_name.c_str());
    DofMap &dof_map = _sys.get_dof_map();
    if (!A_Array)
    {
        // Handle the case where the "Az" array does not exist
        libmesh_error_msg("Error: " + array_name + " Boundary Conditions array not found in vtkPointData.");
    }
    vtkDataArray *V_Array = nullptr;
    if (_sys.has_variable("V") && time_stepper.get_time_steps()[time_stepper.get_current_time_step()].get_transient_effect())
    {
        array_name = "V_" + std::to_string(time_stepper.get_current_time_step() + 1);
        vtkDataArray *V_Array = _pointData->GetArray(array_name.c_str());
        if (!V_Array)
        {
            // Handle the case where the "Az" array does not exist
            libmesh_error_msg("Error: V Boundary Conditions array not found in vtkPointData.");
        }
    }
    for (const auto &node : mesh.local_node_ptr_range())
    {
        if (_dim == 2)
        {
            // Get the scalar value corresponding to node->id()
            Real az_bc = A_Array->GetComponent(node->id(), 0);
            // Now you can use the variable val as needed
            // Check if val is NaN
            if (!std::isnan(az_bc))
            {
                // Apply the row constraint
                dof_id_type dof_id = node->dof_number(_sys.number(), _sys.variable_number("AZ"), 0);
                DofConstraintRow constraint_row;
                dof_map.add_constraint_row(dof_id, constraint_row, az_bc, false);
            }
        }
        else if (_dim == 3)
        {
            // Get the 3-component values corresponding to node->id()
            Real ax_bc = A_Array->GetComponent(node->id(), 0);
            Real ay_bc = A_Array->GetComponent(node->id(), 1);
            Real az_bc = A_Array->GetComponent(node->id(), 2);
            // Now you can use the variables val0, val1, and val2 as needed
            // Apply the row constraint for each component

            if (!std::isnan(ax_bc))
            {
                DofConstraintRow constraint_row;
                dof_id_type dof_id_ax = node->dof_number(_sys.number(), _sys.variable_number("AX"), 0);
                dof_map.add_constraint_row(dof_id_ax, constraint_row, ax_bc, false);
            }
            if (!std::isnan(ay_bc))
            {
                DofConstraintRow constraint_row;
                dof_id_type dof_id_ay = node->dof_number(_sys.number(), _sys.variable_number("AY"), 0);
                dof_map.add_constraint_row(dof_id_ay, constraint_row, ay_bc, false);
            }
            if (!std::isnan(az_bc))
            {
                DofConstraintRow constraint_row;
                dof_id_type dof_id_az = node->dof_number(_sys.number(), _sys.variable_number("AZ"), 0);
                dof_map.add_constraint_row(dof_id_az, constraint_row, az_bc, false);
            }
            if (_sys.has_variable("V") && time_stepper.get_time_steps()[time_stepper.get_current_time_step()].get_transient_effect())
            {
                Real v_bc = V_Array->GetComponent(node->id(), 0);

                if (!std::isnan(v_bc))
                {
                    // Apply the row constraint
                    dof_id_type dof_id = node->dof_number(_sys.number(), _sys.variable_number("V"), 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(dof_id, constraint_row, v_bc, false);
                }
            }
        }
    }
    // If we have VOLT but transient effects are not active, bound V to 0 to avoid crash of the solver
    if (!time_stepper.get_time_steps()[time_stepper.get_current_time_step()].get_transient_effect() && _sys.has_variable("V"))
    {
        // Apply the constraint for the time derivative of the magnetic vector potential
        std::vector<dof_id_type> dof_indices;
        dof_map.local_variable_indices(dof_indices, _sys.get_mesh(), _sys.variable_number("V"));
        for (const auto &dof : dof_indices)
        {
            // Apply the constraint
            DofConstraintRow constraint_row;
            dof_map.add_constraint_row(dof, constraint_row, false);
        }
    }
}

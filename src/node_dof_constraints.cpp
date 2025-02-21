#include "libmesh/equation_systems.h"
#include "libmesh/node.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h" // Include the header file for LinearImplicitSystem
#include "../include/node_dof_constraints.h"
#include <vtkPointData.h>
#include <vtkDataArray.h>

using namespace libMesh;
// BCs are probably broken now

void NodalBCConstraint::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.

    // Retrieve the system by name
    const MeshBase &mesh = _sys.get_mesh();
    std::string array_name = "A_" + std::to_string(timeStepper->get_current_time_step() + 1);
    vtkDataArray *A_Array = _pointData->GetArray(array_name.c_str());
    vtkDataArray *orig_id = _pointData->GetArray("VTKorigID");
    DofMap &dof_map = _sys.get_dof_map();
    if (!A_Array)
    {
        // Handle the case where the "Az" array does not exist
        libmesh_error_msg("Error: " + array_name + " Boundary Conditions array not found in vtkPointData.");
    }
    vtkDataArray *V_Array = nullptr;
    if (_sys.has_variable("V") && timeStepper->get_time_steps()[timeStepper->get_current_time_step()].get_transient_effect())
    {
        array_name = "V_" + std::to_string(timeStepper->get_current_time_step() + 1);
        V_Array = _pointData->GetArray(array_name.c_str());
        if (!V_Array)
        {
            // Handle the case where the "Az" array does not exist
            libmesh_error_msg("Error: V Boundary Conditions array not found in vtkPointData.");
        }
    }
    vtkIdType size_vtk = A_Array->GetNumberOfTuples();
    dof_id_type orig_id_num;
    const Node *node_ptr = nullptr;
    for (vtkIdType j = 0; j < size_vtk; j++)
    {
        orig_id_num = static_cast<dof_id_type>(orig_id->GetTuple1(j));
        if (_dim == 2)
        {
            // Get the scalar value corresponding to node->id()
            Real az_bc = A_Array->GetComponent(j, 0);
            // Now you can use the variable val as needed
            // Check if val is NaN
            if (!std::isnan(az_bc))
            {
                // Apply the row constraint
                node_ptr = mesh.node_ptr(orig_id_num);
                dof_id_type dof_id = node_ptr->dof_number(_sys.number(), _sys.variable_number("A"), 0);
                DofConstraintRow constraint_row;
                dof_map.add_constraint_row(dof_id, constraint_row, az_bc, true);
            }
        }
        else if (_dim == 3)
        {
            // Get the 3-component values corresponding to node->id()
            Real ax_bc = A_Array->GetComponent(j, 0);
            Real ay_bc = A_Array->GetComponent(j, 1);
            Real az_bc = A_Array->GetComponent(j, 2);
            // Now you can use the variables val0, val1, and val2 as needed
            // Apply the row constraint for each component

            if (!std::isnan(ax_bc))
            {
                node_ptr = mesh.node_ptr(orig_id_num);

                DofConstraintRow constraint_row;
                dof_id_type dof_id_ax = node_ptr->dof_number(_sys.number(), _sys.variable_number("A"), 0);
                dof_map.add_constraint_row(dof_id_ax, constraint_row, ax_bc, true);
            }
            if (!std::isnan(ay_bc))
            {
                node_ptr = mesh.node_ptr(orig_id_num);
                DofConstraintRow constraint_row;
                dof_id_type dof_id_ay = node_ptr->dof_number(_sys.number(), _sys.variable_number("A"), 1);
                dof_map.add_constraint_row(dof_id_ay, constraint_row, ay_bc, true);
            }
            if (!std::isnan(az_bc))
            {
                node_ptr = mesh.node_ptr(orig_id_num);

                DofConstraintRow constraint_row;
                dof_id_type dof_id_az = node_ptr->dof_number(_sys.number(), _sys.variable_number("A"), 2);
                dof_map.add_constraint_row(dof_id_az, constraint_row, az_bc, true);
            }
            if (_sys.has_variable("V") && timeStepper->get_time_steps()[timeStepper->get_current_time_step()].get_transient_effect())
            {

                Real v_bc = V_Array->GetTuple1(j);

                if (!std::isnan(v_bc))
                {
                    // Apply the row constraint
                    node_ptr = mesh.node_ptr(orig_id_num);
                    dof_id_type dof_id = node_ptr->dof_number(_sys.number(), _sys.variable_number("V"), 0);
                    DofConstraintRow constraint_row;
                    dof_map.add_constraint_row(dof_id, constraint_row, v_bc, true);
                }
            }
        }
    }
    // If we have VOLT but transient effects are not active, bound V to 0 to avoid crash of the solver
    if (!timeStepper->get_time_steps()[timeStepper->get_current_time_step()].get_transient_effect() && _sys.has_variable("V"))
    {
        // Apply the constraint for the time derivative of the magnetic vector potential
        std::vector<dof_id_type> dof_indices;
        dof_map.local_variable_indices(dof_indices, _sys.get_mesh(), _sys.variable_number("V"));
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
    std::string array_name = "A_" + std::to_string(timeStepper->get_current_time_step() + 1);
    vtkDataArray *A_Array = _pointData->GetArray(array_name.c_str());
    array_name = "V_" + std::to_string(timeStepper->get_current_time_step() + 1);
    vtkDataArray *V_Array = _pointData->GetArray(array_name.c_str());
    array_name = "A_" + std::to_string(timeStepper->get_current_time_step());
    vtkDataArray *prev_A_array = _pointData->GetArray(array_name.c_str());
    array_name = "V_" + std::to_string(timeStepper->get_current_time_step());
    vtkDataArray *prev_V_array = _pointData->GetArray(array_name.c_str());
    Real A_bc;
    Real prev_A_bc;
    Real V_bc;
    Real prev_V_bc;
    vtkIdType size_vtk = A_Array->GetNumberOfTuples();
    char n_comps;

    if (_dim == 2)
    {
        n_comps = 1;
    }
    else
    {
        n_comps = 3;
    }

    if (!A_Array)
    {
        // Handle the case where the "Az" array does not exist
        libmesh_error_msg("Error: " + array_name + " Boundary Conditions array not found in vtkPointData.");
    }
    if (!V_Array)
    {
        // Handle the case where the "Az" array does not exist
        libmesh_error_msg("Error: " + array_name + " Boundary Conditions array not found in vtkPointData.");
    }
    if (!prev_A_array)
    {
        // Handle the case where the "Az" array does not exist
        libmesh_error_msg("Error: " + array_name + " Boundary Conditions array not found in vtkPointData.");
    }
    if (!prev_V_array)
    {
        // Handle the case where the "Az" array does not exist
        libmesh_error_msg("Error: " + array_name + " Boundary Conditions array not found in vtkPointData.");
    }

    for (vtkIdType j = 0; j < size_vtk; j++)
    {
        for (char i = 0; i < n_comps; i++)
        {
            prev_A_bc = prev_A_array->GetComponent(j, i);
            A_bc = A_Array->GetComponent(j, i);
            if ((A_bc != prev_A_bc && !(std::isnan(A_bc) && std::isnan(prev_A_bc))) ||
                (std::isnan(A_bc) && !std::isnan(prev_A_bc)) ||
                (!std::isnan(A_bc) && std::isnan(prev_A_bc)))
            {
                update_bc = true;
            }
        }
        prev_V_bc = prev_V_array->GetTuple1(j);
        V_bc = V_Array->GetTuple1(j);
        if ((V_bc != prev_V_bc && !(std::isnan(V_bc) && std::isnan(prev_V_bc))) ||
            (std::isnan(V_bc) && !std::isnan(prev_V_bc)) ||
            (!std::isnan(V_bc) && std::isnan(prev_V_bc)))
        {
            update_bc = true;
        }
        if (update_bc)
        {
            break;
        }
    }

    _sys.get_mesh().comm().max(update_bc);
    update_boundary_conditions = update_bc;
}
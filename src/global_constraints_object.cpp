#include "fe_parent_object.h"
#include "global_constraints_object.h"
#include "node_loads_dof_constraints.h"
#include "apec_m_system_objects.h"
#include "globals.h"
#include "fem_data_repository.h"
#include "fe_data_object.h"
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>

void GlobalConstraintsObject::constrain()
{
    // Example implementation
    // This function should implement the logic to apply constraints to the nodes of the system
    // identified by system_name within the provided equation_systems object.
    // The actual implementation would depend on the specific constraints to be applied.
    // Retrieve the system by name
    DofMap &dof_map = this->parent_object->get_dof_map();
    if (!this->ctx_view)
    {
        return;
    }
    auto node_bc_object = this->context<NodeLoadsDOFConstraints>();
    if (!node_bc_object)
    {
        return;
    }
    ApecMEquationSystems *data_map = node_bc_object->get_data_map();
    ApecMLinearSystem &system = data_map->get_system<ApecMLinearSystem>(this->parent_object->number());
    NumericVector<Number> &node_dof_constraints = *system.request_vector("node_dof_constraints");
    // std::string string_variable_label;
    for (numeric_index_type k = node_dof_constraints.first_local_index(); k <= node_dof_constraints.last_local_index(); ++k)
    {
        Real value = node_dof_constraints(k);
        if (value != std::numeric_limits<Real>::quiet_NaN())
        {
            DofConstraintRow constraint_row;
            dof_map.add_constraint_row(k, constraint_row, value, true);
        }
    }
}

void GlobalConstraintsObject::check_update_boundary_conditions()
{
    this->update_boundary_conditions = false;
    auto prev_node_bc_object = this->old_context<NodeLoadsDOFConstraints>();
    auto curr_node_bc_object = this->context<NodeLoadsDOFConstraints>();

    if (!prev_node_bc_object || !curr_node_bc_object)
    {
        this->update_boundary_conditions = true;
        return;
    }

    ApecMEquationSystems *data_map = curr_node_bc_object->get_data_map();
    ApecMLinearSystem &system = data_map->get_system<ApecMLinearSystem>(this->parent_object->number());
    NumericVector<Number> &node_dof_constraints = *system.request_vector("node_dof_constraints");
    NumericVector<Number> &old_node_dof_constraints = *system.request_vector("old_node_dof_constraints");
    this->update_boundary_conditions = false;
    for (numeric_index_type k = node_dof_constraints.first_local_index(); k <= node_dof_constraints.last_local_index(); ++k)
    {
        if (node_dof_constraints(k) != old_node_dof_constraints(k))
        {
            this->update_boundary_conditions = true;
            break;
        }
    }
    return;
}

void GlobalConstraintsObject::set_boundary_conditions()
{
    this->check_update_boundary_conditions();
    if (this->update_boundary_conditions)
    {
        parent_object->reinit_constraints();
    }
}
#include "fe_parent_object.h"
#include "global_constraints_object.h"
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
    // at system init
    DofMap &dof_map = this->parent_object->get_dof_map();
    System &system = *this->parent_object;
    // MeshBase &mesh = this->parent_object->get_mesh();
    EquationSystems &equation_systems = system.get_equation_systems();
    if (!this->ctx_view)
    {
        return;
    }
    libMesh::out << "Applying constraints to the system." << std::endl;

    const MeshBase &mesh = equation_systems.get_mesh();
    libMesh::out << "Setting transient effects constraints..." << std::endl;
    auto &time_step = *this->context<TimeStep>();

    // If we have VOLT but transient effects are not active, bound V to 0 to avoid crash of the solver
    if (!time_step.get_transient_effect() && system.has_variable("VOLT"))
    {
        // Apply the constraint for the time derivative of the magnetic vector potential
        std::vector<dof_id_type> dof_indices;
        dof_map.local_variable_indices(dof_indices, mesh, system.variable_number("VOLT"));
        for (const auto &dof : dof_indices)
        {
            // Apply the constraint
            DofConstraintRow constraint_row;
            if (!dof_map.is_constrained_dof(dof))
                dof_map.add_constraint_row(dof, constraint_row, true);
        }
    }

    NumericVector<Number> &node_dof_constraints = *this->get_parent_object()->request_vector("ramped_node_dof_constraints");

    for (numeric_index_type k = node_dof_constraints.first_local_index(); k < node_dof_constraints.last_local_index(); ++k)
    {
        Real value = node_dof_constraints(k);
        if (!std::isnan(value))
        {
            DofConstraintRow constraint_row;
            if (!dof_map.is_constrained_dof(k))
                dof_map.add_constraint_row(k, constraint_row, value, true);
        }
    }
}

void GlobalConstraintsObject::check_update_boundary_conditions()
{
    this->update_boundary_conditions = false;
    DofMap &dof_map = this->parent_object->get_dof_map();
    const auto &comm = dof_map.comm();

    // Swap to bring the old constraints to the front
    dof_map.swap_dof_constraints();
    const auto constraints_map_old = dof_map.get_dof_constraints();

    // Swap back to bring the new constraints to the front
    dof_map.swap_dof_constraints();
    const auto &constraints_map_new = dof_map.get_dof_constraints();

    const auto keys_differ = [constraints_map_new, constraints_map_old]() -> bool
    {
        const bool has_old = !constraints_map_old.empty();
        const bool has_new = !constraints_map_new.empty();

        if (!has_old && !has_new)
            return false;

        if (has_old != has_new)
            return true;

        auto size_1 = constraints_map_old.size();
        auto size_2 = constraints_map_new.size();
        if (constraints_map_old.size() != constraints_map_new.size())
            return true;

        for (const auto &entry : constraints_map_new)
        {
            if (!constraints_map_old.count(entry.first))
                return true;
        }
        return false;
    };

    if (keys_differ())
    {
        this->update_boundary_conditions = true;
        libMesh::out << "Boundary conditions update required: different constrained DOFs" << std::endl;
    }

    comm.max(this->update_boundary_conditions);

    return;
}

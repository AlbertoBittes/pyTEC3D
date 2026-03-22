#include "global_constraints_object.h"
#include "fe_parent_object.h"
#include "apec_m_system.h"

#include "variables_map.h"
#include <libmesh/dof_map.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

void GlobalConstraintsObject::constrain()
{
    DofMap &dof_map = this->parent_object->get_dof_map();
    if (!this->ctx_view)
    {
        return;
    }
    libMesh::out << "Applying constraints to the system." << std::endl;

    libMesh::out << "Setting transient effects constraints..." << std::endl;

    std::map<dof_id_type, ConstraintEntry> constraints;
    auto time_integration_constraints = build_time_integration_constraints();
    constraints.insert(time_integration_constraints.begin(), time_integration_constraints.end());

    libMesh::out << "Setting point constraints..." << std::endl;
    auto point_constraints = build_point_constraints();
    for (const auto &entry : point_constraints)
        constraints.emplace(entry.first, entry.second);

    libMesh::out << "Setting coupling constraints..." << std::endl;
    auto coupling_constraints = build_coupling_constraints();
    for (const auto &entry : coupling_constraints)
        constraints.emplace(entry.first, entry.second);

    for (const auto &entry : constraints)
    {
        if (dof_map.is_constrained_dof(entry.first))
            continue;
        dof_map.add_constraint_row(entry.first, entry.second.row, entry.second.value, true);
    }
}

bool GlobalConstraintsObject::check_update_boundary_conditions()
{
    libMesh::out << "Checking constrained DOFs..." << std::endl;

    bool update_boundary_conditions = false;
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
        update_boundary_conditions = true;
        libMesh::out << "Boundary conditions update required: different constrained DOFs" << std::endl;
    }

    comm.max(update_boundary_conditions);

    if (update_boundary_conditions)
    {
        ApecMLinearSystem &system = static_cast<ApecMLinearSystem &>(*this->get_parent_object());
        system.delete_rows_cols_constrained_dofs();
    }

    return update_boundary_conditions;
}

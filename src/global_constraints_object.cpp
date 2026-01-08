#include "global_constraints_object.h"
#include "fe_parent_object.h"
#include "apec_m.h"
#include "apec_m_mesh.h"
#include "apec_m_system.h"
#include "apec_m_tensor_value.h"
#include "dof_couplings_factory.h"
#include "coordinate_system.h"

#include "fem_data_repository.h"
#include "fe_data_object.h"
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

std::map<dof_id_type, GlobalConstraintsObject::ConstraintEntry>
GlobalConstraintsObject::build_time_integration_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    System *system = this->get_parent_object();
    DofMap &dof_map = system->get_dof_map();
    auto time_step = this->context<TimeStep>();

    EquationSystems &equation_systems = system->get_equation_systems();
    const MeshBase &mesh = equation_systems.get_mesh();
    if (time_step->get_transient_effect())
        return constraints;

    const std::array<std::string, 2> variable_labels = {"VOLT", "EMF"};
    std::vector<dof_id_type> dof_indices;
    for (const auto &label : variable_labels)
    {
        if (!system->has_variable(label))
            continue;
        dof_indices.clear();
        dof_map.local_variable_indices(dof_indices, mesh, system->variable_number(label));
        for (const auto &dof : dof_indices)
            constraints.emplace(dof, ConstraintEntry{});
    }
    return constraints;
}

std::map<dof_id_type, GlobalConstraintsObject::ConstraintEntry>
GlobalConstraintsObject::build_point_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;
    if (!this->parent_object)
        return constraints;

    System *system = this->get_parent_object();
    ApecMEquationSystems &equation_systems = dynamic_cast<ApecMEquationSystems &>(system->get_equation_systems());
    ApecMMesh *mesh = dynamic_cast<ApecMMesh *>(&equation_systems.get_mesh());

    auto node_dof_constraints = system->request_vector("ramped_node_dof_constraints");
    if (!node_dof_constraints)
        return constraints;

    NumericVector<Number> &node_dof_constraints_ref = *node_dof_constraints;
    const auto first_local = node_dof_constraints_ref.first_local_index();
    const auto last_local = node_dof_constraints_ref.last_local_index();

    std::map<std::string, std::vector<VariablesMap::VectorComponent>> vector_variables =
        VariablesMap::get_vector_variables_for_system(*system);
    std::set<dof_id_type> constrained_dofs;

    for (const auto &node_ptr : mesh->local_node_ptr_range())
    {

        bool rotation_ready = false;
        ApecMTensorValue rotation;

        for (const auto &var_entry : vector_variables)
        {
            const auto &components = var_entry.second;
            const unsigned int max_components = std::min<unsigned int>(3, components.size());

            VectorValue<unsigned int> dofs(libMesh::invalid_uint,
                                           libMesh::invalid_uint,
                                           libMesh::invalid_uint);
            VectorValue<Real> values(std::numeric_limits<Real>::quiet_NaN(),
                                     std::numeric_limits<Real>::quiet_NaN(),
                                     std::numeric_limits<Real>::quiet_NaN());

            for (unsigned int i = 0; i < max_components; ++i)
            {
                const auto &comp = components[i];
                dof_id_type dof = libMesh::invalid_uint;
                try
                {
                    dof = node_ptr->dof_number(system->number(), comp.var_number, 0);
                }
                catch (...)
                {
                    continue;
                }
                if (dof < first_local || dof >= last_local)
                    continue;

                dofs(comp.component_index) = dof;
                values(comp.component_index) = node_dof_constraints_ref(dof);
            }

            bool any_constraint = false;
            for (unsigned int i = 0; i < max_components; ++i)
            {
                if (!std::isnan(values(i)))
                {
                    any_constraint = true;
                    break;
                }
            }
            if (!any_constraint)
                continue;

            bool full_vector_constraint = (max_components == 3);
            for (unsigned int i = 0; i < max_components && full_vector_constraint; ++i)
            {
                if (dofs(i) == libMesh::invalid_uint || std::isnan(values(i)))
                    full_vector_constraint = false;
            }

            if (!rotation_ready)
            {
                const auto coord_system = mesh->node_coordinate_system(*node_ptr);
                rotation = coord_system.global_to_local_rotation(*node_ptr);
            }

            if (full_vector_constraint)
            {
                const VectorValue<Real> global_values = rotation * values;

                for (unsigned int i = 0; i < max_components; ++i)
                {
                    ConstraintEntry entry;
                    entry.value = global_values(i);
                    if (constraints.emplace(dofs(i), entry).second)
                        constrained_dofs.insert(dofs(i));
                }
                continue;
            }

            std::vector<unsigned int> constrained_components;
            constrained_components.reserve(max_components);
            for (unsigned int i = 0; i < max_components; ++i)
            {
                if (!std::isnan(values(i)))
                    constrained_components.push_back(i);
            }

            if (constrained_components.size() == 2)
            {
                const unsigned int eq0 = constrained_components[0];
                const unsigned int eq1 = constrained_components[1];
                const ApecMVector2Value rhs(values(eq0), values(eq1));

                std::vector<unsigned int> candidate_dofs;
                candidate_dofs.reserve(max_components);
                for (unsigned int j = 0; j < max_components; ++j)
                {
                    if (dofs(j) == libMesh::invalid_uint)
                        continue;
                    if (constrained_dofs.count(dofs(j)) != 0)
                        continue;
                    candidate_dofs.push_back(j);
                }

                const Real tol = 1e-12;
                bool handled = false;
                if (candidate_dofs.size() >= 2)
                {
                    Real best_det = 0.;
                    unsigned int s0 = 0;
                    unsigned int s1 = 0;
                    for (size_t a = 0; a < candidate_dofs.size(); ++a)
                    {
                        for (size_t b = a + 1; b < candidate_dofs.size(); ++b)
                        {
                            const unsigned int j0 = candidate_dofs[a];
                            const unsigned int j1 = candidate_dofs[b];
                            const ApecMTensor2Value sub(rotation, j0, j1, eq0, eq1);
                            const Real det = sub.det();
                            const Real abs_det = std::abs(det);
                            if (abs_det > best_det)
                            {
                                best_det = abs_det;
                                s0 = j0;
                                s1 = j1;
                            }
                        }
                    }

                    if (best_det >= tol)
                    {
                        const ApecMTensor2Value sub(rotation, s0, s1, eq0, eq1);
                        ApecMTensor2Value inv;
                        if (!sub.inverse(inv, tol))
                            continue;
                        const Real inv00 = inv(0, 0);
                        const Real inv01 = inv(0, 1);
                        const Real inv10 = inv(1, 0);
                        const Real inv11 = inv(1, 1);

                        const unsigned int slaves[2] = {s0, s1};
                        for (unsigned int k = 0; k < 2; ++k)
                        {
                            const unsigned int slave = slaves[k];
                            const Real inv0 = (k == 0) ? inv00 : inv10;
                            const Real inv1 = (k == 0) ? inv01 : inv11;

                            ConstraintEntry entry;
                            for (unsigned int j = 0; j < max_components; ++j)
                            {
                                if (j == s0 || j == s1)
                                    continue;
                                if (dofs(j) == libMesh::invalid_uint)
                                    continue;
                                const Real coeff = -(inv0 * rotation(j, eq0) + inv1 * rotation(j, eq1));
                                entry.row[dofs(j)] = coeff;
                            }

                            entry.value = inv0 * rhs(0) + inv1 * rhs(1);
                            const dof_id_type slave_dof = dofs(slave);
                            if (constraints.emplace(slave_dof, entry).second)
                                constrained_dofs.insert(slave_dof);
                        }

                        handled = true;
                    }
                }

                if (handled)
                    continue;
            }

            for (unsigned int i = 0; i < max_components; ++i)
            {
                const VectorValue<Real> column = rotation.column(i);

                const Real constraint_value = values(i);
                if (std::isnan(constraint_value))
                    continue;

                unsigned int pivot_index = i;
                Real pivot_coeff = column(pivot_index);

                if (dofs(pivot_index) == libMesh::invalid_uint || std::abs(pivot_coeff) < 1e-12)
                {
                    Real best = 0.;
                    for (unsigned int j = 0; j < max_components; ++j)
                    {
                        if (dofs(j) == libMesh::invalid_uint)
                            continue;
                        const Real abs_coeff = std::abs(column(j));
                        if (abs_coeff > best)
                        {
                            best = abs_coeff;
                            pivot_coeff = column(j);
                            pivot_index = j;
                        }
                    }
                }

                if (dofs(pivot_index) == libMesh::invalid_uint || std::abs(pivot_coeff) < 1e-12)
                    continue;

                ConstraintEntry entry;
                for (unsigned int j = 0; j < max_components; ++j)
                {
                    if (j == pivot_index)
                        continue;
                    if (dofs(j) == libMesh::invalid_uint)
                        continue;
                    entry.row[dofs(j)] = -column(j) / pivot_coeff;
                }

                entry.value = constraint_value / pivot_coeff;
                dof_id_type master_dof = dofs(pivot_index);
                if (constraints.emplace(master_dof, entry).second)
                    constrained_dofs.insert(master_dof);
            }
        }
    }

    for (numeric_index_type k = first_local; k < last_local; ++k)
    {
        Real value = node_dof_constraints_ref(k);
        if (!std::isnan(value) && constrained_dofs.count(k) == 0)
        {
            ConstraintEntry entry;
            entry.value = value;
            constraints.emplace(k, entry);
        }
    }

    return constraints;
}

std::map<dof_id_type, GlobalConstraintsObject::ConstraintEntry>
GlobalConstraintsObject::build_coupling_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    System *system = this->get_parent_object();
    ApecMEquationSystems *equation_systems = dynamic_cast<ApecMEquationSystems *>(&system->get_equation_systems());
    ApecMMesh *mesh = dynamic_cast<ApecMMesh *>(&equation_systems->get_mesh());

    auto *program = equation_systems->get_parent_object();

    auto repo = program->model_data();

    auto couplings_object = repo->get_object<DofCouplingsObject>();

    const auto &entries = couplings_object->get_couplings();
    if (entries.empty())
        return constraints;

    const auto &pairs_by_cp = mesh->get_cp_pairs();
    if (pairs_by_cp.empty())
        return constraints;

    std::map<std::string, std::vector<VariablesMap::VectorComponent>> vector_variables =
        VariablesMap::get_vector_variables_for_system(*system);
    // TODO: the coupling part should be done after serializing the mesh
    for (const auto &entry : entries)
    {
        if (entry.variable.empty())
        {
            libMesh::out << "Warning: DOF coupling entry for CP set " << entry.cp_id
                         << " has no variable defined." << std::endl;
            continue;
        }

        auto pairs_it = pairs_by_cp.find(entry.cp_id);
        if (pairs_it == pairs_by_cp.end() || pairs_it->second.empty())
            continue;

        std::set<std::string> coupled_labels;
        for (const auto &label : VariablesMap::get_global_labels_for_label(entry.variable))
            coupled_labels.insert(label);
        if (coupled_labels.empty())
            continue;

        std::map<std::string, std::vector<VariablesMap::VectorComponent>> coupled_vectors;
        for (const auto &vec_entry : vector_variables)
        {
            const auto &components = vec_entry.second;
            if (components.size() < 3)
                continue;

            bool all_present = true;
            for (const auto &component : components)
            {
                if (!coupled_labels.count(component.global_label))
                {
                    all_present = false;
                    break;
                }
            }

            if (!all_present)
                continue;

            coupled_vectors.emplace(vec_entry.first, components);
            for (const auto &component : components)
                coupled_labels.erase(component.global_label);
        }

        for (const auto &pair : pairs_it->second)
        {
            const Node *node_first = mesh->node_ptr(pair.first);
            const Node *node_second = mesh->node_ptr(pair.second);

            if (!node_first || !node_second)
                continue;

            const auto coord_first = mesh->node_coordinate_system(*node_first);
            const auto coord_second = mesh->node_coordinate_system(*node_second);

            for (const auto &vec_entry : coupled_vectors)
            {
                const auto &components = vec_entry.second;
                VectorValue<unsigned int> dofs_first(libMesh::invalid_uint,
                                                     libMesh::invalid_uint,
                                                     libMesh::invalid_uint);
                VectorValue<unsigned int> dofs_second(libMesh::invalid_uint,
                                                      libMesh::invalid_uint,
                                                      libMesh::invalid_uint);

                bool complete = true;
                for (const auto &component : components)
                {
                    if (component.component_index >= 3)
                        continue;

                    try
                    {
                        dofs_first(component.component_index) =
                            node_first->dof_number(system->number(), component.var_number, 0);
                        dofs_second(component.component_index) =
                            node_second->dof_number(system->number(), component.var_number, 0);
                    }
                    catch (...)
                    {
                        complete = false;
                    }
                }

                for (unsigned int i = 0; i < 3; ++i)
                {
                    if (dofs_first(i) == libMesh::invalid_uint || dofs_second(i) == libMesh::invalid_uint)
                    {
                        complete = false;
                        break;
                    }
                }

                if (!complete)
                    continue;

                const Point &point_first = *node_first;
                const Point &point_second = *node_second;
                TensorValue<Real> transform =
                    coord_second.global_to_local_rotation(point_second) *
                    coord_first.local_to_global_rotation(point_first);

                for (unsigned int i = 0; i < 3; ++i)
                {
                    const dof_id_type dof_slave = dofs_second(i);

                    ConstraintEntry constraint_entry;
                    for (unsigned int j = 0; j < 3; ++j)
                        constraint_entry.row[dofs_first(j)] = transform(i, j);
                    constraints.emplace(dof_slave, constraint_entry);
                }
            }

            for (const auto &label : coupled_labels)
            {
                const unsigned int var_num = system->variable_number(label);
                dof_id_type dof_first = libMesh::invalid_uint;
                dof_id_type dof_second = libMesh::invalid_uint;
                try
                {
                    dof_first = node_first->dof_number(system->number(), var_num, 0);
                    dof_second = node_second->dof_number(system->number(), var_num, 0);
                }
                catch (...)
                {
                    continue;
                }

                if (dof_second == libMesh::invalid_uint)
                    continue;
                if (dof_first == libMesh::invalid_uint)
                    continue;

                ConstraintEntry constraint_entry;
                constraint_entry.row[dof_first] = 1.0;
                constraints.emplace(dof_second, constraint_entry);
            }
        }
    }

    return constraints;
}

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

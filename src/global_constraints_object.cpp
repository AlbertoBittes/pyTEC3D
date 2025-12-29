#include "global_constraints_object.h"
#include "fe_parent_object.h"
#include "apec_m.h"
#include "apec_m_system_objects.h"
#include "dof_couplings_factory.h"
#include "globals.h"
#include "fem_data_repository.h"
#include "fe_data_object.h"
#include "variables_map.h"
#include <libmesh/dof_map.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/tensor_value.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace
{
    struct ConstraintEntry
    {
        DofConstraintRow row;
        Real value = 0.0;
    };

    struct VectorComponent
    {
        unsigned int var_number = libMesh::invalid_uint;
        unsigned int component_index = libMesh::invalid_uint;
        std::string api_label;
        std::string global_label;
    };

    struct AngleDataIndices
    {
        std::array<unsigned int, 3> indices{{libMesh::invalid_uint,
                                             libMesh::invalid_uint,
                                             libMesh::invalid_uint}};
        bool available = false;
    };

    AngleDataIndices detect_angle_indices(const MeshBase &mesh)
    {
        AngleDataIndices info;
        bool all_found = true;
        for (unsigned int i = 0; i < 3; ++i)
        {
            const std::string name = "angles-" + std::to_string(i);
            if (mesh.has_node_integer(name))
                info.indices[i] = mesh.get_node_integer_index(name);
            else
            {
                all_found = false;
                break;
            }
        }
        if (all_found)
            info.available = true;
        return info;
    }

    Real snap_trig_value(Real value)
    {
        const Real eps = 1e-12;
        if (std::abs(value) < eps)
            return 0.;
        if (std::abs(value - 1.) < eps)
            return 1.;
        if (std::abs(value + 1.) < eps)
            return -1.;
        return value;
    }

    libMesh::TensorValue<Real> build_rotation_matrix(const std::array<Real, 3> &angles_deg)
    {
        const Real to_rad = libMesh::pi / 180.0;
        const Real ax = angles_deg[0] * to_rad;
        const Real ay = angles_deg[1] * to_rad;
        const Real az = angles_deg[2] * to_rad;

        const Real cx = snap_trig_value(std::cos(ax));
        const Real sx = snap_trig_value(std::sin(ax));
        const Real cy = snap_trig_value(std::cos(ay));
        const Real sy = snap_trig_value(std::sin(ay));
        const Real cz = snap_trig_value(std::cos(az));
        const Real sz = snap_trig_value(std::sin(az));

        // Rotation order Z (yaw), Y (pitch), X (roll)
        libMesh::TensorValue<Real> rot;
        rot(0, 0) = cz * cy;
        rot(0, 1) = cz * sy * sx - sz * cx;
        rot(0, 2) = cz * sy * cx + sz * sx;
        rot(1, 0) = sz * cy;
        rot(1, 1) = sz * sy * sx + cz * cx;
        rot(1, 2) = sz * sy * cx - cz * sx;
        rot(2, 0) = -sy;
        rot(2, 1) = cy * sx;
        rot(2, 2) = cy * cx;
        return rot;
    }

    GlobalVariableName resolve_global_variable(const std::string &var_name)
    {
        auto keys = VariablesMap::get_keys_by_predicate(
            VariablesMap::GLOBAL_VARIABLE_INFO_MAP,
            [&var_name](const VariableInfo &info)
            {
                return info.global_label == var_name;
            });
        if (keys.empty())
            return GlobalVariableName::UNKNOWN;
        if (keys.size() > 1)
            libmesh_error_msg("Variable name maps to multiple global variables: " + var_name);
        return keys.front();
    }

    std::map<dof_id_type, ConstraintEntry> build_volt_constraints(const MeshBase &mesh,
                                                                  System &system,
                                                                  DofMap &dof_map,
                                                                  const TimeStep &time_step)
    {
        std::map<dof_id_type, ConstraintEntry> constraints;
        if (time_step.get_transient_effect() || !system.has_variable("VOLT"))
            return constraints;

        std::vector<dof_id_type> dof_indices;
        dof_map.local_variable_indices(dof_indices, mesh, system.variable_number("VOLT"));
        for (const auto &dof : dof_indices)
        {
            constraints.emplace(dof, ConstraintEntry{});
        }
        return constraints;
    }

    std::map<dof_id_type, ConstraintEntry> build_point_constraints(const MeshBase &mesh,
                                                                   System &system,
                                                                   NumericVector<Number> &node_dof_constraints)
    {
        std::map<dof_id_type, ConstraintEntry> constraints;
        const auto first_local = node_dof_constraints.first_local_index();
        const auto last_local = node_dof_constraints.last_local_index();

        std::map<std::string, std::vector<VectorComponent>> vector_variables;
        for (unsigned int var_num = 0; var_num < system.n_vars(); ++var_num)
        {
            const std::string &var_name = system.variable_name(var_num);
            GlobalVariableName global_var = resolve_global_variable(var_name);
            if (global_var == GlobalVariableName::UNKNOWN)
                continue;

            const auto &var_info = VariablesMap::get_info_by_key(global_var);
            auto components_info = VariablesMap::get_variable_components_for_system(system, global_var);
            if (components_info.is_scalar)
                continue;

            vector_variables[var_info.api_label].push_back(
                {var_num, var_info.component_index, var_info.api_label, var_info.global_label});
        }

        for (auto &entry : vector_variables)
        {
            auto &components = entry.second;
            std::sort(components.begin(), components.end(),
                      [](const VectorComponent &a, const VectorComponent &b)
                      { return a.component_index < b.component_index; });
        }

        AngleDataIndices angle_data = detect_angle_indices(mesh);
        auto read_angle = [&](const Node &node, unsigned int idx) -> Real
        {
            if (!angle_data.available || idx == libMesh::invalid_uint)
                return 0.;
            return node.get_extra_datum<Real>(idx);
        };

        std::set<dof_id_type> constrained_dofs;

        for (const auto &node_ptr : mesh.local_node_ptr_range())
        {
            const Node &node = *node_ptr;

            bool rotation_ready = false;
            std::array<Real, 3> angles_deg{{0., 0., 0.}};
            libMesh::TensorValue<Real> rotation;

            for (const auto &var_entry : vector_variables)
            {
                const auto &components = var_entry.second;
                const unsigned int max_components = std::min<unsigned int>(3, components.size());

                std::array<dof_id_type, 3> dofs{libMesh::invalid_uint,
                                                libMesh::invalid_uint,
                                                libMesh::invalid_uint};
                std::array<Real, 3> values{std::numeric_limits<Real>::quiet_NaN(),
                                           std::numeric_limits<Real>::quiet_NaN(),
                                           std::numeric_limits<Real>::quiet_NaN()};

                for (unsigned int i = 0; i < max_components; ++i)
                {
                    const auto &comp = components[i];
                    dof_id_type dof = libMesh::invalid_uint;
                    try
                    {
                        dof = node.dof_number(system.number(), comp.var_number, 0);
                    }
                    catch (...)
                    {
                        continue;
                    }
                    if (dof < first_local || dof >= last_local)
                        continue;

                    dofs[comp.component_index] = dof;
                    values[comp.component_index] = node_dof_constraints(dof);
                }

                bool any_constraint = false;
                for (unsigned int i = 0; i < max_components; ++i)
                {
                    if (!std::isnan(values[i]))
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
                    if (dofs[i] == libMesh::invalid_uint || std::isnan(values[i]))
                        full_vector_constraint = false;
                }

                if (!rotation_ready)
                {
                    if (angle_data.available)
                    {
                        for (unsigned int i = 0; i < 3; ++i)
                            angles_deg[i] = read_angle(node, angle_data.indices[i]);
                    }
                    rotation = build_rotation_matrix(angles_deg).transpose();
                    rotation_ready = true;
                }

                if (full_vector_constraint)
                {
                    std::array<Real, 3> global_values{
                        rotation(0, 0) * values[0] + rotation(0, 1) * values[1] + rotation(0, 2) * values[2],
                        rotation(1, 0) * values[0] + rotation(1, 1) * values[1] + rotation(1, 2) * values[2],
                        rotation(2, 0) * values[0] + rotation(2, 1) * values[1] + rotation(2, 2) * values[2],
                    };

                    for (unsigned int i = 0; i < max_components; ++i)
                    {
                        ConstraintEntry entry;
                        entry.value = global_values[i];
                        if (constraints.emplace(dofs[i], entry).second)
                            constrained_dofs.insert(dofs[i]);
                    }
                    continue;
                }

                for (unsigned int i = 0; i < max_components; ++i)
                {
                    std::array<Real, 3> column{{rotation(0, i),
                                                rotation(1, i),
                                                rotation(2, i)}};

                    const Real constraint_value = values[i];
                    if (std::isnan(constraint_value))
                        continue;

                    unsigned int pivot_index = i;
                    Real pivot_coeff = column[pivot_index];

                    if (dofs[pivot_index] == libMesh::invalid_uint || std::abs(pivot_coeff) < 1e-12)
                    {
                        Real best = 0.;
                        for (unsigned int j = 0; j < max_components; ++j)
                        {
                            if (dofs[j] == libMesh::invalid_uint)
                                continue;
                            const Real abs_coeff = std::abs(column[j]);
                            if (abs_coeff > best)
                            {
                                best = abs_coeff;
                                pivot_coeff = column[j];
                                pivot_index = j;
                            }
                        }
                    }

                    if (dofs[pivot_index] == libMesh::invalid_uint || std::abs(pivot_coeff) < 1e-12)
                        continue;

                    ConstraintEntry entry;
                    for (unsigned int j = 0; j < max_components; ++j)
                    {
                        if (j == pivot_index)
                            continue;
                        if (dofs[j] == libMesh::invalid_uint)
                            continue;
                        entry.row[dofs[j]] = -column[j] / pivot_coeff;
                    }

                    entry.value = constraint_value / pivot_coeff;
                    dof_id_type master_dof = dofs[pivot_index];
                    if (constraints.emplace(master_dof, entry).second)
                        constrained_dofs.insert(master_dof);
                }
            }
        }

        for (numeric_index_type k = first_local; k < last_local; ++k)
        {
            Real value = node_dof_constraints(k);
            if (!std::isnan(value) && constrained_dofs.count(k) == 0)
            {
                ConstraintEntry entry;
                entry.value = value;
                constraints.emplace(k, entry);
            }
        }

        return constraints;
    }

    std::shared_ptr<DofCouplingsObject> resolve_dof_couplings(const EquationSystems &equation_systems)
    {
        auto *apec_systems = dynamic_cast<const ApecMEquationSystems *>(&equation_systems);
        if (!apec_systems)
            return nullptr;

        auto *program = apec_systems->get_parent_object();
        if (!program)
            return nullptr;

        auto repo = program->model_data();
        if (!repo)
            return nullptr;

        return repo->get_object<DofCouplingsObject>();
    }

    std::map<dof_id_type, ConstraintEntry> build_coupling_constraints(const MeshBase &mesh,
                                                                      System &system,
                                                                      const DofCouplingsObject &couplings,
                                                                      numeric_index_type first_local,
                                                                      numeric_index_type last_local)
    {
        std::map<dof_id_type, ConstraintEntry> constraints;
        const auto &entries = couplings.get_couplings();
        if (entries.empty())
            return constraints;

        std::map<std::string, std::vector<VectorComponent>> vector_variables;
        for (unsigned int var_num = 0; var_num < system.n_vars(); ++var_num)
        {
            const std::string &var_name = system.variable_name(var_num);
            GlobalVariableName global_var = resolve_global_variable(var_name);
            if (global_var == GlobalVariableName::UNKNOWN)
                continue;

            const auto &var_info = VariablesMap::get_info_by_key(global_var);
            auto components_info = VariablesMap::get_variable_components_for_system(system, global_var);
            if (components_info.is_scalar)
                continue;

            vector_variables[var_info.api_label].push_back(
                {var_num, var_info.component_index, var_info.api_label, var_info.global_label});
        }

        for (auto &entry : vector_variables)
        {
            auto &components = entry.second;
            std::sort(components.begin(), components.end(),
                      [](const VectorComponent &a, const VectorComponent &b)
                      { return a.component_index < b.component_index; });
        }

        AngleDataIndices angle_data = detect_angle_indices(mesh);
        auto read_angle = [&](const Node &node, unsigned int idx) -> Real
        {
            if (!angle_data.available || idx == libMesh::invalid_uint)
                return 0.;
            return static_cast<Real>(node.get_extra_datum<dof_id_type>(idx));
        };

        auto node_rotation = [&](const Node &node) -> libMesh::TensorValue<Real>
        {
            std::array<Real, 3> angles_deg{{0., 0., 0.}};
            if (angle_data.available)
            {
                for (unsigned int i = 0; i < 3; ++i)
                    angles_deg[i] = read_angle(node, angle_data.indices[i]);
            }
            return build_rotation_matrix(angles_deg).transpose();
        };

        for (const auto &entry : entries)
        {
            std::set<std::string> coupled_labels(entry.coupled_variables.begin(),
                                                 entry.coupled_variables.end());
            if (coupled_labels.empty())
                continue;

            std::map<std::string, std::vector<VectorComponent>> coupled_vectors;
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

            for (const auto &pair : entry.pairs)
            {
                const Node *node_first = mesh.node_ptr(pair.first);
                const Node *node_second = mesh.node_ptr(pair.second);
                if (!node_first || !node_second)
                    continue;

                const auto rotation_first = node_rotation(*node_first);
                const auto rotation_second = node_rotation(*node_second);

                for (const auto &vec_entry : coupled_vectors)
                {
                    const auto &components = vec_entry.second;
                    std::array<dof_id_type, 3> dofs_first{libMesh::invalid_uint,
                                                          libMesh::invalid_uint,
                                                          libMesh::invalid_uint};
                    std::array<dof_id_type, 3> dofs_second{libMesh::invalid_uint,
                                                           libMesh::invalid_uint,
                                                           libMesh::invalid_uint};

                    bool complete = true;
                    for (const auto &component : components)
                    {
                        if (component.component_index >= 3)
                            continue;

                        try
                        {
                            dofs_first[component.component_index] =
                                node_first->dof_number(system.number(), component.var_number, 0);
                            dofs_second[component.component_index] =
                                node_second->dof_number(system.number(), component.var_number, 0);
                        }
                        catch (...)
                        {
                            complete = false;
                        }
                    }

                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        if (dofs_first[i] == libMesh::invalid_uint ||
                            dofs_second[i] == libMesh::invalid_uint)
                        {
                            complete = false;
                            break;
                        }
                    }

                    if (!complete)
                        continue;

                    Real transform[3][3];
                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        for (unsigned int j = 0; j < 3; ++j)
                        {
                            transform[i][j] =
                                rotation_second(0, i) * rotation_first(0, j) +
                                rotation_second(1, i) * rotation_first(1, j) +
                                rotation_second(2, i) * rotation_first(2, j);
                        }
                    }

                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        const dof_id_type dof_slave = dofs_second[i];
                        if (dof_slave < first_local || dof_slave >= last_local)
                            continue;

                        ConstraintEntry constraint_entry;
                        for (unsigned int j = 0; j < 3; ++j)
                            constraint_entry.row[dofs_first[j]] = transform[i][j];
                        constraints.emplace(dof_slave, constraint_entry);
                    }
                }

                for (const auto &label : coupled_labels)
                {
                    if (!system.has_variable(label))
                        continue;

                    const unsigned int var_num = system.variable_number(label);
                    dof_id_type dof_first = libMesh::invalid_uint;
                    dof_id_type dof_second = libMesh::invalid_uint;
                    try
                    {
                        dof_first = node_first->dof_number(system.number(), var_num, 0);
                        dof_second = node_second->dof_number(system.number(), var_num, 0);
                    }
                    catch (...)
                    {
                        continue;
                    }

                    if (dof_second == libMesh::invalid_uint ||
                        dof_second < first_local || dof_second >= last_local)
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
} // namespace

void GlobalConstraintsObject::constrain()
{
    DofMap &dof_map = this->parent_object->get_dof_map();
    System &system = *this->parent_object;
    EquationSystems &equation_systems = system.get_equation_systems();
    if (!this->ctx_view)
    {
        return;
    }
    libMesh::out << "Applying constraints to the system." << std::endl;

    const MeshBase &mesh = equation_systems.get_mesh();
    libMesh::out << "Setting transient effects constraints..." << std::endl;
    auto &time_step = *this->context<TimeStep>();

    NumericVector<Number> &node_dof_constraints = *this->get_parent_object()->request_vector("ramped_node_dof_constraints");
    NumericVector<Number> *solution_vec = system.current_local_solution.get();
    if (!solution_vec)
        solution_vec = system.solution.get();

    std::map<dof_id_type, ConstraintEntry> constraints;
    auto volt_constraints = build_volt_constraints(mesh, system, dof_map, time_step);
    constraints.insert(volt_constraints.begin(), volt_constraints.end());

    auto point_constraints = build_point_constraints(mesh, system, node_dof_constraints);
    for (const auto &entry : point_constraints)
        constraints.emplace(entry.first, entry.second);

    auto couplings_object = resolve_dof_couplings(equation_systems);
    if (couplings_object && solution_vec)
    {
        auto coupling_constraints = build_coupling_constraints(
            mesh,
            system,
            *couplings_object,
            solution_vec->first_local_index(),
            solution_vec->last_local_index());
        for (const auto &entry : coupling_constraints)
            constraints.emplace(entry.first, entry.second);
    }

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

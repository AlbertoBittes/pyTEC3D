#include "global_constraints_object.h"
#include "fe_parent_object.h"
#include "apec_m_system.h"
#include "apec_m_equation_systems.h"
#include "apec_m_solver.h"
#include "augment_sparsity_on_nodes.h"

#include "variables_map.h"
#include <libmesh/dof_map.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace py = pybind11;

namespace
{
    using RotationCache = std::unordered_map<vtkIdType, TensorValue<Real>>;

    const py::object &coordinate_system_class()
    {
        static py::object cls = []() -> py::object
        {
            try
            {
                return py::module_::import("apec_m.coordinate_system").attr("CoordinateSystem");
            }
            catch (const py::error_already_set &)
            {
                // Backward-compat fallback for legacy path setups.
                return py::module_::import("coordinate_system").attr("CoordinateSystem");
            }
        }();
        return cls;
    }

    void build_rotation_cache_from_point_ids(
        vtkDataArray *angles_array,
        const std::vector<vtkIdType> &point_ids_for_rotation,
        RotationCache &local_to_global_cache)
    {
        if (!angles_array || point_ids_for_rotation.empty())
            return;

        py::gil_scoped_acquire gil;
        const py::object &coord_system = coordinate_system_class();
        const std::vector<py::ssize_t> angles_shape = {
            static_cast<py::ssize_t>(point_ids_for_rotation.size()),
            static_cast<py::ssize_t>(3)};
        py::array_t<Real> angles_np(angles_shape);
        auto angles_view = angles_np.mutable_unchecked<2>();
        for (py::ssize_t idx = 0; idx < static_cast<py::ssize_t>(point_ids_for_rotation.size()); ++idx)
        {
            const vtkIdType point_id = point_ids_for_rotation[static_cast<std::size_t>(idx)];
            angles_view(idx, 0) = static_cast<Real>(angles_array->GetComponent(point_id, 0));
            angles_view(idx, 1) = static_cast<Real>(angles_array->GetComponent(point_id, 1));
            angles_view(idx, 2) = static_cast<Real>(angles_array->GetComponent(point_id, 2));
        }

        const py::object coeffs_obj =
            coord_system.attr("build_local_to_global_coefficients")(angles_np);
        const py::array_t<Real, py::array::c_style | py::array::forcecast> coeffs_np =
            py::cast<py::array_t<Real, py::array::c_style | py::array::forcecast>>(coeffs_obj);
        const auto coeffs = coeffs_np.unchecked<2>();

        for (py::ssize_t idx = 0; idx < static_cast<py::ssize_t>(point_ids_for_rotation.size()); ++idx)
        {
            const vtkIdType point_id = point_ids_for_rotation[static_cast<std::size_t>(idx)];
            local_to_global_cache.emplace(
                point_id,
                TensorValue<Real>(
                    coeffs(idx, 0), coeffs(idx, 1), coeffs(idx, 2),
                    coeffs(idx, 3), coeffs(idx, 4), coeffs(idx, 5),
                    coeffs(idx, 6), coeffs(idx, 7), coeffs(idx, 8)));
        }
    }

} // namespace

std::map<dof_id_type, GlobalConstraintsObject::ConstraintEntry>
GlobalConstraintsObject::build_time_integration_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    System *system = this->get_parent_object();
    auto *linear_system = dynamic_cast<ApecMLinearSystem *>(system);

    if (py::cast<bool>(linear_system->time_step->time_step_dict["transient_effects"]))
        return constraints;

    auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&system->get_equation_systems());
    auto *solver = equation_systems ? equation_systems->get_parent_object() : nullptr;

    std::vector<unsigned int> option4_values;
    {
        py::gil_scoped_acquire gil;
        const py::object element_types = solver->get_apec_object().attr("element_types");
        option4_values = py::cast<std::vector<unsigned int>>(
            element_types.attr("__getitem__")("OPTION4").attr("to_numpy")());
    }

    bool has_timeintegrated_formulation = false;
    for (const unsigned int opt4 : option4_values)
    {
        if (opt4 == 1)
        {
            has_timeintegrated_formulation = true;
            break;
        }
    }
    if (!has_timeintegrated_formulation)
        return constraints;

    DofMap &dof_map = system->get_dof_map();
    const MeshBase &mesh = system->get_equation_systems().get_mesh();

    const std::array<std::string, 2> variable_labels = {"VOLT", "EMF"};
    std::vector<dof_id_type> dof_indices;
    for (const auto &label : variable_labels)
    {
        if (!system->has_variable(label))
            continue;

        dof_indices.clear();
        dof_map.local_variable_indices(dof_indices, mesh, system->variable_number(label));
        for (const dof_id_type dof : dof_indices)
            constraints.emplace(dof, ConstraintEntry{});
    }

    return constraints;
}

std::map<dof_id_type, GlobalConstraintsObject::ConstraintEntry>
GlobalConstraintsObject::build_point_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    System *system = this->get_parent_object();
    auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&system->get_equation_systems());
    auto *solver = equation_systems ? equation_systems->get_parent_object() : nullptr;
    vtkUnstructuredGrid *vtk_grid = solver ? solver->active_vtk_grid : nullptr;

    vtkPointData *point_data = vtk_grid->GetPointData();

    vtkDataArray *angles_array = vtkDataArray::SafeDownCast(point_data->GetArray("angles"));

    auto vector_variables = VariablesMap::get_vector_variables_for_system(*system);
    std::unordered_map<std::string, vtkDataArray *> vector_arrays;
    std::unordered_set<unsigned int> vector_var_numbers;
    for (const auto &entry : vector_variables)
    {
        const std::string array_name = "_" + entry.first;
        vector_arrays.emplace(entry.first, vtkDataArray::SafeDownCast(point_data->GetArray(array_name.c_str())));
        for (const auto &component : entry.second)
            vector_var_numbers.insert(component.var_number);
    }

    struct ScalarArrayEntry
    {
        unsigned int var_number = libMesh::invalid_uint;
        vtkDataArray *array = nullptr;
        unsigned int component_index = 0;
    };

    std::vector<ScalarArrayEntry> scalar_arrays;
    const unsigned int n_vars = system->n_vars();
    scalar_arrays.reserve(n_vars);
    for (unsigned int var_num = 0; var_num < n_vars; ++var_num)
    {
        if (vector_var_numbers.count(var_num))
            continue;

        const std::string &var_name = system->variable_name(var_num);
        const GlobalVariableName global_var = VariablesMap::get_key_by_predicate_one(
            [&var_name](const VariableInfo &info)
            {
                return info.global_label == var_name;
            });
        const auto &var_info = VariablesMap::get_info_by_key(global_var);

        vtkDataArray *array = nullptr;
        std::string array_name = "_" + var_info.api_label;
        array = vtkDataArray::SafeDownCast(point_data->GetArray(array_name.c_str()));

        scalar_arrays.push_back({var_num, array, var_info.component_index});
    }

    const Real tol = std::numeric_limits<Real>::epsilon();
    std::unordered_set<dof_id_type> constrained_dofs;
    RotationCache local_to_global_cache;

    std::vector<vtkIdType> point_ids_for_rotation;
    std::unordered_set<vtkIdType> seen_point_ids;
    const vtkIdType n_angle_tuples = angles_array->GetNumberOfTuples();

    for (const auto &node_ptr : equation_systems->get_mesh().local_node_ptr_range())
    {
        const vtkIdType point_id = static_cast<vtkIdType>(node_ptr->id());
        // if (point_id < 0 || point_id >= n_angle_tuples)
        //     continue;

        bool needs_rotation = false;
        for (const auto &var_entry : vector_variables)
        {
            auto array_it = vector_arrays.find(var_entry.first);
            if (array_it == vector_arrays.end() || array_it->second == nullptr)
                continue;

            vtkDataArray *constraint_array = array_it->second;
            const auto &components = var_entry.second;
            const unsigned int max_components = std::min<unsigned int>(3, components.size());
            for (unsigned int i = 0; i < max_components; ++i)
            {
                const unsigned int comp_idx = components[i].component_index;
                const Real value = static_cast<Real>(constraint_array->GetComponent(point_id, comp_idx));
                if (!std::isnan(value))
                {
                    needs_rotation = true;
                    break;
                }
            }

            if (needs_rotation)
                break;
        }

        if (needs_rotation)
        {
            if (point_id < 0 || point_id >= n_angle_tuples)
                continue;
            if (seen_point_ids.insert(point_id).second)
                point_ids_for_rotation.push_back(point_id);
        }
    }

    build_rotation_cache_from_point_ids(angles_array, point_ids_for_rotation, local_to_global_cache);

    for (const auto &node_ptr : equation_systems->get_mesh().local_node_ptr_range())
    {
        const vtkIdType point_id = static_cast<vtkIdType>(node_ptr->id());

        bool rotation_ready = false;
        TensorValue<Real> global_to_local(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0);
        TensorValue<Real> local_to_global(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0);

        // Constraints for vector variables
        for (const auto &var_entry : vector_variables)
        {
            auto array_it = vector_arrays.find(var_entry.first);
            if (array_it == vector_arrays.end() || array_it->second == nullptr)
                continue;
            vtkDataArray *constraint_array = array_it->second;

            // Get the number of components for this variable
            const auto &components = var_entry.second;
            const unsigned int max_components = std::min<unsigned int>(3, components.size());

            // Initialize dof numbers and the values of the candidate constraint dofs
            std::array<dof_id_type, 3> dofs = {
                libMesh::invalid_uint,
                libMesh::invalid_uint,
                libMesh::invalid_uint};
            std::array<Real, 3> values = {
                std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN(),
                std::numeric_limits<Real>::quiet_NaN()};

            // Retrieve the value of the constraints
            for (unsigned int i = 0; i < max_components; ++i)
            {
                const auto &component = components[i];

                const unsigned int sys_num = system->number();
                if (component.var_number == libMesh::invalid_uint ||
                    component.var_number >= node_ptr->n_vars(sys_num) ||
                    node_ptr->n_dofs(sys_num, component.var_number) == 0)
                    continue;

                dof_id_type dof = node_ptr->dof_number(sys_num, component.var_number, 0);

                dofs[component.component_index] = dof;
                values[component.component_index] =
                    static_cast<Real>(constraint_array->GetComponent(point_id, component.component_index));
            }

            // Check if there is at least one component constrained.
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

            // If all 3 components are constrained, we can directly apply the constraint without worrying about linear dependence.
            bool full_vector_constraint = (max_components == 3);
            for (unsigned int i = 0; i < max_components && full_vector_constraint; ++i)
            {
                if (dofs[i] == libMesh::invalid_uint || std::isnan(values[i]))
                    full_vector_constraint = false;
            }

            // Apply rotation to constraint values
            if (!rotation_ready)
            {
                auto cache_it = local_to_global_cache.find(point_id);
                if (cache_it == local_to_global_cache.end())
                    libmesh_error_msg("Point constraint rotation matrix cache miss for node id " << point_id << ".");
                local_to_global = cache_it->second;
                global_to_local = TensorValue<Real>(local_to_global.transpose());
                rotation_ready = true;
            }

            // If all components are constrained, we can directly apply the rotated constraints
            // resulting from the explicit rotation of the original constraints, without worrying about linear dependence.
            if (full_vector_constraint)
            {
                VectorValue<Real> local_values(values[0], values[1], values[2]);
                const VectorValue<Real> global_values(local_to_global * local_values);

                for (unsigned int i = 0; i < 3; ++i)
                {
                    ConstraintEntry entry;
                    entry.value = global_values(i);
                    if (constraints.emplace(dofs[i], entry).second)
                        constrained_dofs.insert(dofs[i]);
                }
                continue;
            }

            // If not all components of the vector variable are constrained
            // we need to check for linear dependence and possibly reduce the constraints to a subset of independent ones.
            std::array<unsigned int, 3> constrained_components = {0, 0, 0};
            unsigned int n_constrained_components = 0;
            for (unsigned int i = 0; i < max_components; ++i)
                if (!std::isnan(values[i]))
                    constrained_components[n_constrained_components++] = i;

            if (n_constrained_components == 2)
            {
                const unsigned int eq0 = constrained_components[0];
                const unsigned int eq1 = constrained_components[1];
                const Real rhs0 = values[eq0];
                const Real rhs1 = values[eq1];

                std::array<unsigned int, 3> candidate_dofs = {0, 0, 0};
                unsigned int n_candidate_dofs = 0;
                for (unsigned int j = 0; j < max_components; ++j)
                {
                    if (dofs[j] == libMesh::invalid_uint)
                        continue;
                    if (constrained_dofs.count(dofs[j]) != 0)
                        continue;
                    candidate_dofs[n_candidate_dofs++] = j;
                }

                bool handled = false;
                if (n_candidate_dofs >= 2)
                {
                    Real best_det = 0.0;
                    unsigned int s0 = 0;
                    unsigned int s1 = 0;
                    for (unsigned int a = 0; a < n_candidate_dofs; ++a)
                    {
                        for (unsigned int b = a + 1; b < n_candidate_dofs; ++b)
                        {
                            const unsigned int j0 = candidate_dofs[a];
                            const unsigned int j1 = candidate_dofs[b];
                            const Real a00 = global_to_local(eq0, j0);
                            const Real a01 = global_to_local(eq0, j1);
                            const Real a10 = global_to_local(eq1, j0);
                            const Real a11 = global_to_local(eq1, j1);
                            const Real det = a00 * a11 - a01 * a10;
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
                        const Real a00 = global_to_local(eq0, s0);
                        const Real a01 = global_to_local(eq0, s1);
                        const Real a10 = global_to_local(eq1, s0);
                        const Real a11 = global_to_local(eq1, s1);
                        const Real det = a00 * a11 - a01 * a10;
                        if (std::abs(det) < tol)
                            continue;

                        const Real inv00 = a11 / det;
                        const Real inv01 = -a01 / det;
                        const Real inv10 = -a10 / det;
                        const Real inv11 = a00 / det;

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
                                if (dofs[j] == libMesh::invalid_uint)
                                    continue;
                                const Real coeff =
                                    -(inv0 * global_to_local(eq0, j) + inv1 * global_to_local(eq1, j));
                                entry.row[dofs[j]] = coeff;
                            }

                            entry.value = inv0 * rhs0 + inv1 * rhs1;
                            const dof_id_type slave_dof = dofs[slave];
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
                const Real constraint_value = values[i];
                if (std::isnan(constraint_value))
                    continue;

                unsigned int pivot_index = i;
                Real pivot_coeff = global_to_local(i, pivot_index);

                if (dofs[pivot_index] == libMesh::invalid_uint || std::abs(pivot_coeff) < tol)
                {
                    Real best = 0.0;
                    for (unsigned int j = 0; j < max_components; ++j)
                    {
                        if (dofs[j] == libMesh::invalid_uint)
                            continue;
                        const Real abs_coeff = std::abs(global_to_local(i, j));
                        if (abs_coeff > best)
                        {
                            best = abs_coeff;
                            pivot_coeff = global_to_local(i, j);
                            pivot_index = j;
                        }
                    }
                }

                if (dofs[pivot_index] == libMesh::invalid_uint || std::abs(pivot_coeff) < tol)
                    continue;

                ConstraintEntry entry;
                for (unsigned int j = 0; j < max_components; ++j)
                {
                    if (j == pivot_index)
                        continue;
                    if (dofs[j] == libMesh::invalid_uint)
                        continue;
                    entry.row[dofs[j]] = -global_to_local(i, j) / pivot_coeff;
                }

                entry.value = constraint_value / pivot_coeff;
                const dof_id_type constrained_dof = dofs[pivot_index];
                if (constraints.emplace(constrained_dof, entry).second)
                    constrained_dofs.insert(constrained_dof);
            }
        }

        // Handle scalar variable constraints
        for (const auto &scalar_entry : scalar_arrays)
        {
            const unsigned int sys_num = system->number();
            if (scalar_entry.var_number == libMesh::invalid_uint ||
                scalar_entry.var_number >= node_ptr->n_vars(sys_num) ||
                node_ptr->n_dofs(sys_num, scalar_entry.var_number) == 0)
                continue;

            dof_id_type dof = node_ptr->dof_number(sys_num, scalar_entry.var_number, 0);

            if (constrained_dofs.count(dof) != 0)
                continue;

            const Real value = static_cast<Real>(
                scalar_entry.array->GetComponent(point_id, scalar_entry.component_index));
            if (std::isnan(value))
                continue;

            ConstraintEntry entry;
            entry.value = value;
            if (constraints.emplace(dof, entry).second)
                constrained_dofs.insert(dof);
        }
    }

    return constraints;
}

std::map<dof_id_type, GlobalConstraintsObject::ConstraintEntry>
GlobalConstraintsObject::build_coupling_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    System *system = this->get_parent_object();
    auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&system->get_equation_systems());
    auto *solver = equation_systems ? equation_systems->get_parent_object() : nullptr;
    vtkUnstructuredGrid *vtk_grid = solver ? solver->active_vtk_grid : nullptr;

    vtkPointData *point_data = vtk_grid->GetPointData();

    vtkDataArray *angles_array = vtkDataArray::SafeDownCast(point_data->GetArray("angles"));
    std::vector<int> cp_id_values;
    std::vector<std::string> variable_values;
    {
        py::gil_scoped_acquire gil;
        py::object dof_couplings_df = solver->get_apec_object().attr("dof_couplings");
        cp_id_values = py::cast<std::vector<int>>(
            dof_couplings_df.attr("__getitem__")("CP_ID").attr("to_numpy")());
        variable_values = py::cast<std::vector<std::string>>(
            dof_couplings_df.attr("__getitem__")("VARIABLE").attr("to_numpy")());
    }
    const std::size_t n_entries = cp_id_values.size();
    if (n_entries == 0)
        return constraints;

    const auto &pairs_by_cp = solver->coupling_functor->cp_pairs_by_cp;
    if (pairs_by_cp.empty())
        return constraints;

    auto vector_variables = VariablesMap::get_vector_variables_for_system(*system);
    const MeshBase &mesh = equation_systems->get_mesh();
    const unsigned int proc_id = mesh.comm().rank();
    RotationCache local_to_global_cache;

    const auto build_coupled_groups =
        [&vector_variables](
            const std::string &variable,
            std::map<std::string, std::vector<VariablesMap::VectorComponent>> &coupled_vectors,
            std::unordered_set<std::string> &coupled_labels)
    {
        coupled_labels.clear();
        for (const auto &label : VariablesMap::get_global_labels_for_label(variable))
            coupled_labels.insert(label);

        coupled_vectors.clear();
        if (coupled_labels.empty())
            return;

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
    };

    struct LocalPair
    {
        const Node *node_first = nullptr;
        const Node *node_second = nullptr;
        vtkIdType point_id_first = -1;
        vtkIdType point_id_second = -1;
    };

    struct CouplingEntryPlan
    {
        std::map<std::string, std::vector<VariablesMap::VectorComponent>> coupled_vectors;
        std::unordered_set<std::string> coupled_labels;
        std::vector<LocalPair> local_pairs;
    };

    std::vector<CouplingEntryPlan> entry_plans;
    entry_plans.reserve(n_entries);
    std::vector<vtkIdType> point_ids_for_rotation;
    std::unordered_set<vtkIdType> seen_point_ids;
    const vtkIdType n_angle_tuples = angles_array ? angles_array->GetNumberOfTuples() : 0;

    for (std::size_t entry_idx = 0; entry_idx < n_entries; ++entry_idx)
    {
        const int cp_id = cp_id_values[entry_idx];
        const std::string &variable = variable_values[entry_idx];
        if (variable.empty())
            continue;

        auto pairs_it = pairs_by_cp.find(cp_id);
        if (pairs_it == pairs_by_cp.end() || pairs_it->second.empty())
            continue;

        CouplingEntryPlan plan;
        build_coupled_groups(variable, plan.coupled_vectors, plan.coupled_labels);
        if (plan.coupled_vectors.empty() && plan.coupled_labels.empty())
            continue;

        plan.local_pairs.reserve(pairs_it->second.size());
        const bool needs_rotation = !plan.coupled_vectors.empty();

        for (const auto &pair : pairs_it->second)
        {
            const Node *node_first = mesh.query_node_ptr(pair.first);
            const Node *node_second = mesh.query_node_ptr(pair.second);
            if (!node_second)
                continue;
            if (node_second->processor_id() != proc_id)
                continue;
            if (!node_first)
                libmesh_error_msg("Coupling constraint: slave node " << pair.second
                                                                     << " found without master node " << pair.first
                                                                     << ".");

            const vtkIdType point_id_first = static_cast<vtkIdType>(node_first->id());
            const vtkIdType point_id_second = static_cast<vtkIdType>(node_second->id());
            plan.local_pairs.push_back({node_first, node_second, point_id_first, point_id_second});

            if (needs_rotation && angles_array)
            {
                if (point_id_first >= 0 && point_id_first < n_angle_tuples && seen_point_ids.insert(point_id_first).second)
                    point_ids_for_rotation.push_back(point_id_first);
                if (point_id_second >= 0 && point_id_second < n_angle_tuples && seen_point_ids.insert(point_id_second).second)
                    point_ids_for_rotation.push_back(point_id_second);
            }
        }

        if (!plan.local_pairs.empty())
            entry_plans.emplace_back(std::move(plan));
    }

    build_rotation_cache_from_point_ids(angles_array, point_ids_for_rotation, local_to_global_cache);

    for (const auto &plan : entry_plans)
    {
        for (const auto &pair : plan.local_pairs)
        {
            const Node *node_first = pair.node_first;
            const Node *node_second = pair.node_second;
            const vtkIdType point_id_first = pair.point_id_first;
            const vtkIdType point_id_second = pair.point_id_second;

            for (const auto &vec_entry : plan.coupled_vectors)
            {
                const auto &components = vec_entry.second;
                std::array<dof_id_type, 3> dofs_first = {
                    libMesh::invalid_uint,
                    libMesh::invalid_uint,
                    libMesh::invalid_uint};
                std::array<dof_id_type, 3> dofs_second = {
                    libMesh::invalid_uint,
                    libMesh::invalid_uint,
                    libMesh::invalid_uint};

                bool complete = true;
                const unsigned int sys_num = system->number();
                for (const auto &component : components)
                {
                    if (component.var_number == libMesh::invalid_uint ||
                        component.var_number >= node_first->n_vars(sys_num) ||
                        component.var_number >= node_second->n_vars(sys_num) ||
                        node_first->n_dofs(sys_num, component.var_number) == 0 ||
                        node_second->n_dofs(sys_num, component.var_number) == 0)
                    {
                        complete = false;
                        break;
                    }

                    dofs_first[component.component_index] =
                        node_first->dof_number(sys_num, component.var_number, 0);
                    dofs_second[component.component_index] =
                        node_second->dof_number(sys_num, component.var_number, 0);
                }

                for (unsigned int i = 0; i < 3; ++i)
                {
                    if (dofs_first[i] == libMesh::invalid_uint || dofs_second[i] == libMesh::invalid_uint)
                    {
                        complete = false;
                        break;
                    }
                }
                if (!complete)
                    continue;

                TensorValue<Real> first_local_to_global(
                    1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0);
                auto first_cache_it = local_to_global_cache.find(point_id_first);
                if (first_cache_it == local_to_global_cache.end())
                    libmesh_error_msg("Coupling constraint master rotation matrix cache miss for node id " << point_id_first << ".");
                first_local_to_global = first_cache_it->second;
                const TensorValue<Real> first_global_to_local =
                    TensorValue<Real>(first_local_to_global.transpose());

                TensorValue<Real> second_local_to_global(
                    1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0);
                auto second_cache_it = local_to_global_cache.find(point_id_second);
                if (second_cache_it == local_to_global_cache.end())
                    libmesh_error_msg("Coupling constraint slave rotation matrix cache miss for node id " << point_id_second << ".");
                second_local_to_global = second_cache_it->second;

                const TensorValue<Real> transform =
                    TensorValue<Real>(second_local_to_global * first_global_to_local);

                for (unsigned int i = 0; i < 3; ++i)
                {
                    ConstraintEntry constraint_entry;
                    for (unsigned int j = 0; j < 3; ++j)
                        constraint_entry.row[dofs_first[j]] =
                            transform(i, j);
                    constraints.emplace(dofs_second[i], constraint_entry);
                }
            }

            for (const auto &label : plan.coupled_labels)
            {
                const unsigned int sys_num = system->number();
                const unsigned int var_num = system->variable_number(label);
                if (var_num >= node_first->n_vars(sys_num) ||
                    var_num >= node_second->n_vars(sys_num) ||
                    node_first->n_dofs(sys_num, var_num) == 0 ||
                    node_second->n_dofs(sys_num, var_num) == 0)
                    continue;

                const dof_id_type dof_first = node_first->dof_number(sys_num, var_num, 0);
                const dof_id_type dof_second = node_second->dof_number(sys_num, var_num, 0);

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
    if (!this->get_parent_object()->is_initialized())
        return;
    DofMap &dof_map = this->parent_object->get_dof_map();
    libMesh::out << "Applying constraints to the system." << std::endl;

    libMesh::out << "Setting transient effects constraints..." << std::endl;

    std::unordered_map<dof_id_type, ConstraintEntry> constraints;
    auto time_integration_constraints = build_time_integration_constraints();
    for (const auto &entry : time_integration_constraints)
        constraints.emplace(entry.first, entry.second);

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

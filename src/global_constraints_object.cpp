/**
 * @file global_constraints_object.cpp
 * @brief GlobalConstraintsObject: time-integration pins, Dirichlet BCs, and DOF-coupling constraints.
 *
 * Merges three constraint sources each time libMesh calls constrain():
 *   1. Time-integration homogeneous pins (static formulation, OPTION4==1).
 *   2. Point/nodal Dirichlet BCs from VTK point-data arrays (with optional local-to-global rotation).
 *   3. DOF-coupling constraints from CP_* periodic-boundary node pairs.
 */
#include "global_constraints_object.h"
#include "fe_parent_object.h"
#include "apec_m_system.h"
#include "apec_m_equation_systems.h"
#include "apec_m_model.h"
#include "augment_sparsity_on_nodes.h"
#include "coordinate_system.h"

#include "variables_map.h"
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_base.h>
#include <libmesh/node.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include "apec_m_logging.h"

namespace py = pybind11;

using libMesh::dof_id_type;
using libMesh::DofMap;
using libMesh::Elem;
using libMesh::MeshBase;
using libMesh::Node;
using libMesh::Number;
using libMesh::NumericVector;
using libMesh::Point;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::RealVectorValue;
using libMesh::SparseMatrix;
using libMesh::System;
using libMesh::TensorValue;
using libMesh::VectorValue;

namespace
{
    std::unordered_map<dof_id_type, unsigned int> build_dof_to_variable_index_map(const System &system)
    {
        std::unordered_map<dof_id_type, unsigned int> dof_to_var_idx;
        const auto &mesh = system.get_equation_systems().get_mesh();
        const unsigned int sys_num = system.number();
        const unsigned int n_vars = system.n_vars();

        // Map each local DOF to its variable index for O(1) per-constraint variable lookups.
        for (unsigned int var_num = 0; var_num < n_vars; ++var_num)
        {
            std::vector<dof_id_type> dof_indices;
            system.get_dof_map().local_variable_indices(dof_indices, mesh, var_num);
            for (const dof_id_type dof : dof_indices)
                dof_to_var_idx.emplace(dof, var_num);
        }

        return dof_to_var_idx;
    }

    TensorValue<Real> build_local_to_global_from_angles(
        vtkDataArray *angles_array,
        const vtkIdType point_id,
        const std::string &context)
    {
        if (!angles_array)
        {
            std::ostringstream msg;
            msg << context << ": missing 'angles' point-data array.";
            apec::logging::error_and_abort(msg.str());
        }

        const vtkIdType n_tuples = angles_array->GetNumberOfTuples();
        if (point_id < 0 || point_id >= n_tuples)
        {
            std::ostringstream msg;
            msg << context << ": invalid node id " << point_id
                << " for 'angles' array with " << n_tuples << " tuples.";
            apec::logging::error_and_abort(msg.str());
        }

        double rotation[9];
        apec_m::CoordinateSystem::build_local_to_global_matrix(
            static_cast<double>(angles_array->GetComponent(point_id, 0)),
            static_cast<double>(angles_array->GetComponent(point_id, 1)),
            static_cast<double>(angles_array->GetComponent(point_id, 2)),
            rotation);

        return TensorValue<Real>(
            rotation[0], rotation[1], rotation[2],
            rotation[3], rotation[4], rotation[5],
            rotation[6], rotation[7], rotation[8]);
    }

} // namespace

template <typename ImplSystem>
GlobalConstraintsObject<ImplSystem>::GlobalConstraintsObject(ApecMDOFSystem<ImplSystem> *sys)
    : FEParentObjectBase<ApecMDOFSystem<ImplSystem>>(sys)
{
    this->initialize_coupling_entries_cache();
}

template <typename ImplSystem>
void GlobalConstraintsObject<ImplSystem>::initialize_coupling_entries_cache()
{
    _coupling_cp_id_values.clear();
    _coupling_variable_values.clear();
    _az_coupling_cp_ids.clear();

    ApecMDOFSystem<ImplSystem> *dof_system = this->get_parent_object();
    auto &equation_systems = dof_system->get_apec_equation_systems();
    auto *solver = equation_systems.get_parent_object();
    if (!solver)
        return;

    {
        py::gil_scoped_acquire gil;
        py::object dof_couplings_df = solver->get_apec_object().attr("dof_couplings");
        _coupling_cp_id_values = py::cast<std::vector<int>>(
            dof_couplings_df.attr("__getitem__")("CP_ID").attr("to_numpy")());
        _coupling_variable_values = py::cast<std::vector<std::string>>(
            dof_couplings_df.attr("__getitem__")("VARIABLE").attr("to_numpy")());
    }

    const std::size_t n_entries =
        std::min(_coupling_cp_id_values.size(), _coupling_variable_values.size());
    for (std::size_t i = 0; i < n_entries; ++i)
        if (_coupling_variable_values[i] == "Az")
            _az_coupling_cp_ids.insert(_coupling_cp_id_values[i]);
}

template <typename ImplSystem>
std::map<dof_id_type, typename GlobalConstraintsObject<ImplSystem>::ConstraintEntry>
GlobalConstraintsObject<ImplSystem>::build_time_integration_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    ApecMDOFSystem<ImplSystem> *dof_system = this->get_parent_object();
    System &system = *dof_system;

    const auto *time_step = dof_system->current_time_step();
    if (time_step && time_step->transient_effects)
        return constraints;

    DofMap &dof_map = system.get_dof_map();
    const MeshBase &mesh = system.get_equation_systems().get_mesh();

    const std::array<std::string, 2> variable_labels = {"V", "EMF"};
    std::vector<dof_id_type> dof_indices;
    for (const auto &label : variable_labels)
    {
        if (!system.has_variable(label))
            continue;

        dof_indices.clear();
        dof_map.local_variable_indices(dof_indices, mesh, system.variable_number(label));
        for (const dof_id_type dof : dof_indices)
            constraints.emplace(dof, ConstraintEntry{});
    }

    return constraints;
}

template <typename ImplSystem>
std::map<dof_id_type, typename GlobalConstraintsObject<ImplSystem>::ConstraintEntry>
GlobalConstraintsObject<ImplSystem>::build_point_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;
    this->_point_constraint_owner_cache.clear();

    ApecMDOFSystem<ImplSystem> *dof_system = this->get_parent_object();
    System &system = *dof_system;
    auto &equation_systems = dof_system->get_apec_equation_systems();
    auto *solver = equation_systems.get_parent_object();
    vtkUnstructuredGrid *vtk_grid = solver ? solver->get_active_vtk_grid() : nullptr;

    vtkPointData *point_data = vtk_grid->GetPointData();

    vtkDataArray *angles_array = vtkDataArray::SafeDownCast(point_data->GetArray("angles"));

    auto vector_variables = VariablesMap::get_vector_variables_for_system(system);
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
    const unsigned int n_vars = system.n_vars();
    scalar_arrays.reserve(n_vars);
    for (unsigned int var_num = 0; var_num < n_vars; ++var_num)
    {
        if (vector_var_numbers.count(var_num))
            continue;

        const std::string &var_name = system.variable_name(var_num);
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
    for (const auto &node_ptr : equation_systems.get_mesh().local_node_ptr_range())
    {
        const vtkIdType point_id = static_cast<vtkIdType>(node_ptr->id());
        const dof_id_type node_id = node_ptr->id();

        bool rotation_ready = false;
        TensorValue<Real> global_to_local(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0);
        TensorValue<Real> local_to_global(
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0);

        // Apply Dirichlet constraints for vector-valued variables.
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

                const unsigned int sys_num = system.number();
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
                local_to_global = build_local_to_global_from_angles(
                    angles_array,
                    point_id,
                    "Point constraint rotation");
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
                    {
                        this->_point_constraint_owner_cache[dofs[i]] = DofOwnerInfo{node_id, system.variable_name(components[i].var_number)};
                        constrained_dofs.insert(dofs[i]);
                    }
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
                            {
                                this->_point_constraint_owner_cache[slave_dof] = DofOwnerInfo{node_id, system.variable_name(components[slave].var_number)};
                                constrained_dofs.insert(slave_dof);
                            }
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
                {
                    this->_point_constraint_owner_cache[constrained_dof] = DofOwnerInfo{node_id, system.variable_name(components[pivot_index].var_number)};
                    constrained_dofs.insert(constrained_dof);
                }
            }
        }

        // Handle scalar variable constraints
        for (const auto &scalar_entry : scalar_arrays)
        {
            const unsigned int sys_num = system.number();
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
            {
                this->_point_constraint_owner_cache[dof] = DofOwnerInfo{node_id, system.variable_name(scalar_entry.var_number)};
                constrained_dofs.insert(dof);
            }
        }
    }

    return constraints;
}

template <typename ImplSystem>
std::map<dof_id_type, typename GlobalConstraintsObject<ImplSystem>::ConstraintEntry>
GlobalConstraintsObject<ImplSystem>::build_coupling_constraints()
{
    std::map<dof_id_type, ConstraintEntry> constraints;

    ApecMDOFSystem<ImplSystem> *dof_system = this->get_parent_object();
    System &system = *dof_system;
    auto &equation_systems = dof_system->get_apec_equation_systems();
    auto *solver = equation_systems.get_parent_object();
    vtkUnstructuredGrid *vtk_grid = solver ? solver->get_active_vtk_grid() : nullptr;

    vtkPointData *point_data = vtk_grid->GetPointData();

    vtkDataArray *angles_array = vtkDataArray::SafeDownCast(point_data->GetArray("angles"));
    const std::size_t n_entries =
        std::min(_coupling_cp_id_values.size(), _coupling_variable_values.size());
    if (n_entries == 0)
        return constraints;

    const auto *coupling_functor = solver->get_coupling_functor();
    const auto &pairs_by_cp = coupling_functor->cp_pairs_by_cp;
    if (pairs_by_cp.empty())
        return constraints;

    auto vector_variables = VariablesMap::get_vector_variables_for_system(system);
    const MeshBase &mesh = equation_systems.get_mesh();
    const unsigned int proc_id = mesh.comm().rank();
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

    for (std::size_t entry_idx = 0; entry_idx < n_entries; ++entry_idx)
    {
        const int cp_id = _coupling_cp_id_values[entry_idx];
        const std::string &variable = _coupling_variable_values[entry_idx];
        if (variable.empty())
            continue;

        auto pairs_it = pairs_by_cp.find(cp_id);
        if (pairs_it == pairs_by_cp.end() || pairs_it->second.empty())
            continue;

        std::map<std::string, std::vector<VariablesMap::VectorComponent>> coupled_vectors;
        std::unordered_set<std::string> coupled_labels;
        build_coupled_groups(variable, coupled_vectors, coupled_labels);
        if (coupled_vectors.empty() && coupled_labels.empty())
            continue;

        for (const auto &pair : pairs_it->second)
        {
            const Node *node_first = mesh.query_node_ptr(pair.first);
            const Node *node_second = mesh.query_node_ptr(pair.second);
            if (!node_second)
                continue;
            if (node_second->processor_id() != proc_id)
                continue;
            if (!node_first)
            {
                std::ostringstream msg;
                msg << "Coupling constraint: slave node " << pair.second
                    << " found without master node " << pair.first
                    << ".";
                apec::logging::error_and_abort(msg.str());
            }

            const vtkIdType point_id_first = static_cast<vtkIdType>(node_first->id());
            const vtkIdType point_id_second = static_cast<vtkIdType>(node_second->id());

            for (const auto &vec_entry : coupled_vectors)
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
                const unsigned int sys_num = system.number();
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

                const TensorValue<Real> first_local_to_global =
                    build_local_to_global_from_angles(
                        angles_array,
                        point_id_first,
                        "Coupling constraint master rotation");
                const TensorValue<Real> first_global_to_local =
                    TensorValue<Real>(first_local_to_global.transpose());

                const TensorValue<Real> second_local_to_global =
                    build_local_to_global_from_angles(
                        angles_array,
                        point_id_second,
                        "Coupling constraint slave rotation");

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

            for (const auto &label : coupled_labels)
            {
                const unsigned int sys_num = system.number();
                const unsigned int var_num = system.variable_number(label);
                if (var_num >= node_first->n_vars(sys_num) ||
                    var_num >= node_second->n_vars(sys_num) ||
                    node_first->n_dofs(sys_num, var_num) == 0 ||
                    node_second->n_dofs(sys_num, var_num) == 0)
                    continue;

                const dof_id_type dof_first = node_first->dof_number(sys_num, var_num, 0);
                const dof_id_type dof_second = node_second->dof_number(sys_num, var_num, 0);

                ConstraintEntry constraint_entry;
                bool same_edge_orientation = true;
                if (label == "Az")
                {
                    if (_edge_orientations_map)
                    {
                        const auto orientation_it = _edge_orientations_map->find(
                            std::make_pair(static_cast<dof_id_type>(node_first->id()),
                                           static_cast<dof_id_type>(node_second->id())));
                        if (orientation_it != _edge_orientations_map->end())
                            same_edge_orientation = orientation_it->second;
                    }
                    constraint_entry.row[dof_first] = same_edge_orientation ? 1.0 : -1.0;
                }
                else
                    constraint_entry.row[dof_first] = 1.0;
                constraints.emplace(dof_second, constraint_entry);
            }
        }
    }

    return constraints;
}

template <typename ImplSystem>
void GlobalConstraintsObject<ImplSystem>::build_edge_orientations_map()
{
    ApecMDOFSystem<ImplSystem> *dof_system = this->get_parent_object();
    System &system = *dof_system;
    if (!system.has_variable("Az") || _az_coupling_cp_ids.empty())
    {
        _edge_orientations_map.reset();
        return;
    }

    if (!_edge_orientations_map)
        _edge_orientations_map =
            std::make_unique<std::map<std::pair<dof_id_type, dof_id_type>, bool>>();
    else
        _edge_orientations_map->clear();

    auto &equation_systems = dof_system->get_apec_equation_systems();
    auto *solver = equation_systems.get_parent_object();

    const auto *coupling_functor = solver->get_coupling_functor();
    const auto &pairs_by_cp = coupling_functor->cp_pairs_by_cp;
    if (pairs_by_cp.empty())
        return;

    const MeshBase &mesh = equation_systems.get_mesh();
    const unsigned int proc_id = mesh.comm().rank();

    struct EdgeInfo
    {
        dof_id_type corner0 = libMesh::invalid_uint;
        dof_id_type corner1 = libMesh::invalid_uint;
    };
    std::unordered_map<dof_id_type, EdgeInfo> edge_info_by_midside_node;

    // Build edge-corner lookup on the distributed mesh view (owned + ghosted elements).
    for (const auto *elem : mesh.element_ptr_range())
    {
        for (const auto edge_idx : elem->edge_index_range())
        {
            const unsigned int corner0_local = elem->local_edge_node(edge_idx, 0);
            const unsigned int corner1_local = elem->local_edge_node(edge_idx, 1);
            const dof_id_type corner0 = static_cast<dof_id_type>(elem->node_id(corner0_local));
            const dof_id_type corner1 = static_cast<dof_id_type>(elem->node_id(corner1_local));

            const auto local_nodes_on_edge = elem->nodes_on_edge(edge_idx);
            for (const unsigned int local_node_idx : local_nodes_on_edge)
            {
                if (elem->is_vertex(local_node_idx))
                    continue;

                const auto node_id =
                    static_cast<dof_id_type>(elem->node_id(local_node_idx));
                edge_info_by_midside_node.emplace(node_id, EdgeInfo{corner0, corner1});
            }
        }
    }

    // Build a per-node CP label map from the already assembled CP pairs.
    // Missing labels are treated as 0 during orientation inference.
    std::unordered_map<dof_id_type, int> node_cp_label;
    for (const auto &cp_entry : pairs_by_cp)
    {
        const int cp_label = cp_entry.first;
        for (const auto &cp_pair : cp_entry.second)
        {
            node_cp_label.emplace(cp_pair.first, cp_label);
            node_cp_label.emplace(cp_pair.second, cp_label);
        }
    }

    const auto cp_of = [&node_cp_label](const dof_id_type node_id) -> int
    {
        const auto it = node_cp_label.find(node_id);
        return it == node_cp_label.end() ? 0 : it->second;
    };

    const auto order_corners_by_point =
        [&mesh](const dof_id_type c0, const dof_id_type c1) -> std::array<dof_id_type, 2>
    {
        const Node *n0 = mesh.query_node_ptr(c0);
        const Node *n1 = mesh.query_node_ptr(c1);
        if (!n0 || !n1)
        {
            std::ostringstream msg;
            msg << "Coupling constraint: missing corner node while ordering edge corners by coordinates ("
                << c0 << ", " << c1 << ").";
            apec::logging::error_and_abort(msg.str());
        }

        // Match libMesh positive-edge-orientation style: compare Point coordinates.
        const Point &p0 = *n0;
        const Point &p1 = *n1;
        return (p0 > p1)
                   ? std::array<dof_id_type, 2>{c0, c1}
                   : std::array<dof_id_type, 2>{c1, c0};
    };

    for (const int cp_id : _az_coupling_cp_ids)
    {
        const auto cp_it = pairs_by_cp.find(cp_id);
        if (cp_it == pairs_by_cp.end())
            continue;

        const auto &cp_pairs = cp_it->second;

        for (const auto &pair : cp_pairs)
        {
            const Node *slave_node = mesh.query_node_ptr(pair.second);
            if (!slave_node || slave_node->processor_id() != proc_id)
                continue;

            const Node *master_node = mesh.query_node_ptr(pair.first);
            if (!master_node)
            {
                std::ostringstream msg;
                msg << "Coupling constraint: slave node " << pair.second
                    << " found without master node " << pair.first
                    << ".";
                apec::logging::error_and_abort(msg.str());
            }

            bool same_edge_orientation = true;
            const auto master_it = edge_info_by_midside_node.find(pair.first);
            const auto slave_it = edge_info_by_midside_node.find(pair.second);
            if (master_it != edge_info_by_midside_node.end() &&
                slave_it != edge_info_by_midside_node.end())
            {
                const std::array<dof_id_type, 2> master_corners = order_corners_by_point(
                    master_it->second.corner0,
                    master_it->second.corner1);
                const std::array<dof_id_type, 2> slave_corners = order_corners_by_point(
                    slave_it->second.corner0,
                    slave_it->second.corner1);

                std::array<int, 3> master_tokens = {
                    cp_of(master_corners[0]),
                    cp_of(pair.first),
                    cp_of(master_corners[1]),
                };
                std::array<int, 3> slave_tokens = {
                    cp_of(slave_corners[0]),
                    cp_of(pair.second),
                    cp_of(slave_corners[1]),
                };

                const int master_corner_cp_count =
                    static_cast<int>(master_tokens[0] != 0) +
                    static_cast<int>(master_tokens[2] != 0);
                const int slave_corner_cp_count =
                    static_cast<int>(slave_tokens[0] != 0) +
                    static_cast<int>(slave_tokens[2] != 0);
                if (master_corner_cp_count == 0 && slave_corner_cp_count == 0)
                {
                    std::ostringstream msg;
                    msg << "Coupling constraint: unable to infer Az edge orientation because "
                        << "no corner CP label was found for either edge in CP_ID="
                        << cp_id << " (master midside node " << pair.first
                        << ", slave midside node " << pair.second << ").";
                    apec::logging::error_and_abort(msg.str());
                }

                std::array<int, 3> master_cp_sorted = master_tokens;
                std::array<int, 3> slave_cp_sorted = slave_tokens;
                std::sort(master_cp_sorted.begin(), master_cp_sorted.end());
                std::sort(slave_cp_sorted.begin(), slave_cp_sorted.end());
                if (master_cp_sorted != slave_cp_sorted)
                {
                    std::ostringstream msg;
                    msg << "Coupling constraint: incompatible edge CP token triplets for CP_ID="
                        << cp_id << " (master midside node " << pair.first
                        << ", slave midside node " << pair.second << ").";
                    apec::logging::error_and_abort(msg.str());
                }

                const bool direct_match =
                    (master_tokens[0] == slave_tokens[0]) &&
                    (master_tokens[1] == slave_tokens[1]) &&
                    (master_tokens[2] == slave_tokens[2]);
                const bool reversed_match =
                    (master_tokens[0] == slave_tokens[2]) &&
                    (master_tokens[1] == slave_tokens[1]) &&
                    (master_tokens[2] == slave_tokens[0]);

                const bool direct_missing_ok =
                    !((master_tokens[0] == 0 && slave_tokens[0] == 0 && master_corners[0] != slave_corners[0]) ||
                      (master_tokens[2] == 0 && slave_tokens[2] == 0 && master_corners[1] != slave_corners[1]));
                const bool reversed_missing_ok =
                    !((master_tokens[0] == 0 && slave_tokens[2] == 0 && master_corners[0] != slave_corners[1]) ||
                      (master_tokens[2] == 0 && slave_tokens[0] == 0 && master_corners[1] != slave_corners[0]));

                const bool direct_valid = direct_match && direct_missing_ok;
                const bool reversed_valid = reversed_match && reversed_missing_ok;

                if (!direct_valid && !reversed_valid)
                {
                    std::ostringstream msg;
                    msg << "Coupling constraint: unable to infer Az edge orientation from CP tokens/coordinates for CP_ID="
                        << cp_id << " (master midside node " << pair.first
                        << ", slave midside node " << pair.second
                        << "). If CP is missing on matched corners, the missing-CP corners must be the same node.";
                    apec::logging::error_and_abort(msg.str());
                }

                same_edge_orientation = direct_valid;
            }

            (*_edge_orientations_map)[std::make_pair(pair.first, pair.second)] =
                same_edge_orientation;
        }
    }
}

template <typename ImplSystem>
void GlobalConstraintsObject<ImplSystem>::constrain()
{
    if (!this->get_parent_object()->is_initialized())
        return;

    if (this->get_parent_object()->has_variable("Az") &&
        (_edge_orientations_map == nullptr || _edge_orientations_map->empty()))
        build_edge_orientations_map();

    ApecMDOFSystem<ImplSystem> *dof_system = this->get_parent_object();
    const auto *time_step = dof_system->current_time_step();
    const bool time_step_updated = !time_step || time_step->updated;

    if (!time_step_updated && _has_constraint_cache)
    {
        apec::logging::info("Reusing cached DOF constraints (time step unchanged)");
        DofMap &dof_map = this->parent_object->get_dof_map();
        for (const auto &entry : _constraint_cache)
        {
            if (!dof_map.is_constrained_dof(entry.first))
                dof_map.add_constraint_row(entry.first, entry.second.row, entry.second.value, true);
        }
        return;
    }

    System &system = *this->get_parent_object();
    DofMap &dof_map = this->parent_object->get_dof_map();
    apec::logging::info("Computing nodal DOF constraints");

    const unsigned int n_vars = system.n_vars();
    const auto dof_to_var_idx = build_dof_to_variable_index_map(system);
    std::vector<std::size_t> total_constraints_per_variable(n_vars, 0);
    std::vector<std::size_t> applied_constraints_per_variable(n_vars, 0);
    std::size_t unknown_total_constraints = 0;
    std::size_t unknown_applied_constraints = 0;

    const auto count_total_for_dof = [&dof_to_var_idx,
                                      &total_constraints_per_variable,
                                      &unknown_total_constraints](const dof_id_type dof)
    {
        const auto it = dof_to_var_idx.find(dof);
        if (it == dof_to_var_idx.end())
        {
            ++unknown_total_constraints;
            return;
        }
        ++total_constraints_per_variable[it->second];
    };

    const auto count_applied_for_dof = [&dof_to_var_idx,
                                        &applied_constraints_per_variable,
                                        &unknown_applied_constraints](const dof_id_type dof)
    {
        const auto it = dof_to_var_idx.find(dof);
        if (it == dof_to_var_idx.end())
        {
            ++unknown_applied_constraints;
            return;
        }
        ++applied_constraints_per_variable[it->second];
    };

    std::unordered_map<dof_id_type, ConstraintEntry> constraints;
    const auto add_or_update_constraint = [&constraints,
                                           &count_total_for_dof](const dof_id_type dof,
                                                                 const ConstraintEntry &entry)
    {
        const auto insert_result = constraints.emplace(dof, entry);
        if (!insert_result.second)
        {
            insert_result.first->second = entry;
            return;
        }
        count_total_for_dof(dof);
    };

    auto coupling_constraints = build_coupling_constraints();
    for (const auto &entry : coupling_constraints)
        add_or_update_constraint(entry.first, entry.second);

    auto point_constraints = build_point_constraints();
    std::size_t overridden_couplings = 0;
    for (const auto &entry : point_constraints)
    {
        if (constraints.find(entry.first) != constraints.end() &&
            coupling_constraints.find(entry.first) != coupling_constraints.end())
        {
            ++overridden_couplings;
            const auto owner_it = this->_point_constraint_owner_cache.find(entry.first);
            if (owner_it != this->_point_constraint_owner_cache.end())
            {
                apec::logging::warn(
                    "Point constraint overrides coupling constraint (node={}, variable={}, dof={})",
                    owner_it->second.node_id,
                    owner_it->second.variable,
                    entry.first);
            }
            else
            {
                apec::logging::warn(
                    "Point constraint overrides coupling constraint (node=<unknown>, variable=<unknown>, dof={})",
                    entry.first);
            }
        }
        add_or_update_constraint(entry.first, entry.second);
    }

    auto time_integration_constraints = build_time_integration_constraints();
    for (const auto &entry : time_integration_constraints)
    {
        ConstraintEntry forced_zero = entry.second;
        forced_zero.row.clear();
        forced_zero.value = 0.0;
        add_or_update_constraint(entry.first, forced_zero);
    }

    std::size_t applied_constraints = 0;
    for (const auto &entry : constraints)
    {
        if (dof_map.is_constrained_dof(entry.first))
            continue;
        dof_map.add_constraint_row(entry.first, entry.second.row, entry.second.value, true);
        ++applied_constraints;
        count_applied_for_dof(entry.first);
    }

    const auto &comm = dof_map.comm();
    std::size_t global_time_integration_constraints = time_integration_constraints.size();
    std::size_t global_point_constraints = point_constraints.size();
    std::size_t global_coupling_constraints = coupling_constraints.size();
    std::size_t global_total_constraints = constraints.size();
    std::size_t global_applied_constraints = applied_constraints;
    for (unsigned int var_idx = 0; var_idx < n_vars; ++var_idx)
    {
        comm.sum(total_constraints_per_variable[var_idx]);
        comm.sum(applied_constraints_per_variable[var_idx]);
    }
    comm.sum(unknown_total_constraints);
    comm.sum(unknown_applied_constraints);
    comm.sum(global_time_integration_constraints);
    comm.sum(global_point_constraints);
    comm.sum(global_coupling_constraints);
    comm.sum(global_total_constraints);
    comm.sum(global_applied_constraints);

    std::vector<apec::logging::ConstraintVariableStat> variable_stats;
    variable_stats.reserve(n_vars);
    for (unsigned int var_idx = 0; var_idx < n_vars; ++var_idx)
    {
        const std::size_t total_count = total_constraints_per_variable[var_idx];
        const std::size_t applied_count = applied_constraints_per_variable[var_idx];
        if (total_count == 0 && applied_count == 0)
            continue;
        variable_stats.push_back({system.variable_name(var_idx), total_count, applied_count});
    }
    std::sort(variable_stats.begin(), variable_stats.end(),
              [](const apec::logging::ConstraintVariableStat &lhs,
                 const apec::logging::ConstraintVariableStat &rhs)
              {
                  if (lhs.total_count != rhs.total_count)
                      return lhs.total_count > rhs.total_count;
                  return lhs.variable_name < rhs.variable_name;
              });

    // Cache constraints for reuse in subsequent Picard iterations at the same time step.
    _constraint_cache.clear();
    _constraint_cache.insert(constraints.begin(), constraints.end());
    _has_constraint_cache = true;

    apec::logging::table_root(apec::logging::format_constraints_summary_table(global_time_integration_constraints,
                                                                              global_point_constraints,
                                                                              global_coupling_constraints,
                                                                              global_total_constraints,
                                                                              global_applied_constraints,
                                                                              variable_stats,
                                                                              unknown_total_constraints,
                                                                              unknown_applied_constraints));
}

template <typename ImplSystem>
bool GlobalConstraintsObject<ImplSystem>::check_update_boundary_conditions()
{
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
        apec::logging::info("Boundary conditions update required: constrained DOF set changed");
    }

    comm.max(update_boundary_conditions);

    if (update_boundary_conditions)
        this->get_parent_object()->delete_rows_cols_constrained_dofs();

    return update_boundary_conditions;
}

template class GlobalConstraintsObject<libMesh::LinearImplicitSystem>;
template class GlobalConstraintsObject<libMesh::NonlinearImplicitSystem>;

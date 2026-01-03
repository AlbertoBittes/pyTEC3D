#include "apec_m_system.h"
#include "system_matrix_reduction.h"
#include "element_types_factory.h"
#include "global_constraints_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/numeric_vector.h"
#include "variables_map.h"

const std::map<ApecMSystemBase::SystemType, std::string> ApecMSystemBase::SYSTEM_INFO_MAP = {
    {ApecMSystemBase::DOF_SOLUTION_SYSTEM, "ApecMLinearSystem"},
    {ApecMSystemBase::NODAL_LOADS, "NodalLoadsSystem"},
    {ApecMSystemBase::ELEMENT_LOADS, "ElementLoadsSystem"},
    {ApecMSystemBase::NODAL_SOLUTION, "NodalSolutionSystem"},
    {ApecMSystemBase::ELEMENT_SOLUTION, "ElementSolutionSystem"},
};

bool ApecMSystemBase::is_nodal_system(ApecMSystemBase::SystemType system_type)
{
    return system_type == ApecMSystemBase::DOF_SOLUTION_SYSTEM ||
           system_type == ApecMSystemBase::NODAL_LOADS ||
           system_type == ApecMSystemBase::NODAL_SOLUTION;
}

const std::string &ApecMSystemBase::system_name_from_type(ApecMSystemBase::SystemType system_type)
{
    auto it = ApecMSystemBase::SYSTEM_INFO_MAP.find(system_type);
    if (it == ApecMSystemBase::SYSTEM_INFO_MAP.end())
        libmesh_error_msg("System type not found.");
    return it->second;
}

template <>
void ApecMSystem<System>::ramp_vectors(const Substep &time_substep, NumericVector<Number> &old_vector, NumericVector<Number> &current_vector, NumericVector<Number> &ramped_vector)
{
    if (time_substep.get_ramped_load())
    {
        Real substep_number = static_cast<Real>(time_substep.get_index()) + 1;
        Real number_of_substeps = static_cast<Real>(time_substep.get_parent_object()->get_n_substeps());
        for (numeric_index_type k = old_vector.first_local_index(); k < old_vector.last_local_index(); ++k)
        {
            Real old_val = old_vector(k);
            Real curr_val = current_vector(k);

            if (std::isnan(curr_val))
            {
                ramped_vector.set(k, std::numeric_limits<Real>::quiet_NaN());
            }
            else
            {
                Real start_val = std::isnan(old_val) ? 0.0 : old_val;
                Real ramped_val = start_val + (curr_val - start_val) / number_of_substeps * substep_number;
                ramped_vector.set(k, ramped_val);
            }
        }
    }
    else
    {
        for (numeric_index_type k = current_vector.first_local_index(); k < current_vector.last_local_index(); ++k)
        {
            Real curr_val = current_vector(k);
            ramped_vector.set(k, curr_val);
        }
    }
}

void ApecMLinearSystem::release_subset_solve()
{
    if (this->_subset != nullptr)
    {
        delete this->_subset;
        this->_subset = nullptr;
    }
}

void ApecMLinearSystem::compute_solution_time_derivative(const Substep &time_substep)
{
    bool transient_effect = time_substep.get_transient_effect();
    NumericVector<Number> *time_derivative_vector = this->request_vector("solution_time_derivative");
    if (!transient_effect)
    {
        // If transient effects are not considered, the solution time derivative is zero
        this->older_local_solution->zero();
        time_derivative_vector->zero();
    }
    else
    {
        Real theta = time_substep.get_theta();
        Real dt = time_substep.get_dt();
        Real theta_x_dt = theta * dt;
        Real theta_minus_1_dt = (theta - 1) * dt;

        this->older_local_solution->scale(theta_minus_1_dt);
        *this->older_local_solution += *this->current_local_solution;
        *this->older_local_solution -= *this->old_local_solution;
        this->older_local_solution->scale(1 / theta_x_dt);

        time_derivative_vector->zero();
        *time_derivative_vector += *this->current_local_solution;
        *time_derivative_vector -= *this->old_local_solution;
        time_derivative_vector->scale(1 / dt);
    }
    this->older_local_solution->close();
    time_derivative_vector->close();
}

void ApecMLinearSystem::delete_rows_cols_constrained_dofs()
{
    libMesh::out << "Removing constrained DOFs from the system..." << std::endl;
    SystemMatrixReduction *subset_ptr = new SystemMatrixReduction(*this);
    this->release_subset_solve();
    this->restrict_solve_to(subset_ptr, SUBSET_ZERO);
}

void ApecMLinearSystem::zero_all()
{
    ApecMSystem<TransientSystem<LinearImplicitSystem>>::zero_all();
    this->get_system_matrix().zero();
}

template <typename SysType>
void ApecMSystem<SysType>::init_variables()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;

    auto element_types_collector = context<ElementTypesFactory>();
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    // Use new helper functions to map system name to SystemType
    ApecMSystemBase::SystemType type = ApecMSystemBase::NUM_SYSTEMS;
    const std::string sys_name = this->name();
    for (const auto &pair : ApecMSystemBase::SYSTEM_INFO_MAP)
    {
        if (pair.second == sys_name)
        {
            type = pair.first;
            break;
        }
    }
    if (type == ApecMSystemBase::NUM_SYSTEMS)
        libmesh_error_msg("System type not found for system name: " + sys_name);

    auto vars_to_subdomain_ids = element_types_collector->associate_variables_with_element_IDs(type);

    for (const auto &pair : *vars_to_subdomain_ids)
    {
        const auto &variable_string = VariablesMap::get_info_by_key(pair.first).global_label;
        const std::set<subdomain_id_type> &subdomain_ids = pair.second;
        auto elem_type = element_types_collector->get_element_type(*subdomain_ids.begin());
        auto fe_family = VariablesMap::get_variable_fe_family(*elem_type, pair.first);
        auto order = VariablesMap::get_variable_order(*elem_type, pair.first);

        // Add the variable to the system, specifying subdomains
        this->add_variable(
            variable_string,
            order,
            fe_family,
            &subdomain_ids);
    }
    // Add more conditions as needed for other combinations
}

template class ApecMSystem<System>;
template class ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class ApecMSystem<libMesh::ExplicitSystem>;

#include "apec_m_system.h"

#include <array>
#include <set>
#include <unordered_set>
#include <vector>

#include "system_matrix_reduction.h"
#include "element_types_factory.h"
#include "global_constraints_object.h"
#include "apec_m_mesh.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/numeric_vector.h"
#include "variables_map.h"

const std::map<ApecMSystemBase::SystemType, std::string> ApecMSystemBase::SYSTEM_INFO_MAP = {
    {ApecMSystemBase::DOF_SOLUTION_SYSTEM, "ApecMLinearSystem"},
    {ApecMSystemBase::NODAL_LOADS, "NodalLoadsSystem"},
    {ApecMSystemBase::ELEMENT_LOADS, "ElementLoadsSystem"},
    {ApecMSystemBase::NODAL_SOLUTION, "NodalSolutionSystem"},
    {ApecMSystemBase::ELEMENT_SOLUTION, "ElementSolutionSystem"},
};

const std::string &ApecMSystemBase::system_name_from_type(ApecMSystemBase::SystemType system_type)
{
    auto it = ApecMSystemBase::SYSTEM_INFO_MAP.find(system_type);
    if (it == ApecMSystemBase::SYSTEM_INFO_MAP.end())
        libmesh_error_msg("System type not found.");
    return it->second;
}

void ApecMLinearSystem::release_subset_solve()
{
    if (this->_subset != nullptr)
    {
        delete this->_subset;
        this->_subset = nullptr;
    }
}

void ApecMLinearSystem::compute_solution_time_derivative(const TimeStep &time_step)
{
    bool transient_effect = time_step.get_data_map();
    NumericVector<Number> *time_derivative_vector = this->request_vector("solution_time_derivative");
    if (!transient_effect)
    {
        // If transient effects are not considered, the solution time derivative is zero
        this->older_local_solution->zero();
        time_derivative_vector->zero();
    }
    else
    {
        Real theta = time_step.get_theta();
        Real dt = time_step.get_dt();
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

    if (!element_types_collector || !element_types_collector->get_element_types_map())
        libmesh_error_msg("Element types map not available for system variables.");
    auto vars_to_subdomain_ids = VariablesMap::associate_variables_with_element_IDs(
        *element_types_collector->get_element_types_map(),
        type);

    std::set<subdomain_id_type> mesh_subdomain_ids;
    this->get_mesh().subdomain_ids(mesh_subdomain_ids);

    for (auto &pair : vars_to_subdomain_ids)
    {
        const auto &variable_string = VariablesMap::get_info_by_key(pair.first).global_label;
        auto &subdomain_ids = pair.second;
        for (auto it = subdomain_ids.begin(); it != subdomain_ids.end();)
        {
            if (!mesh_subdomain_ids.count(*it))
                it = subdomain_ids.erase(it);
            else
                ++it;
        }
        if (subdomain_ids.empty())
            continue;
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

template class ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class ApecMSystem<libMesh::ExplicitSystem>;

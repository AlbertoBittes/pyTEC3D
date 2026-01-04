#include "apec_m_system.h"

#include <array>
#include <vector>

#include "system_matrix_reduction.h"
#include "element_types_factory.h"
#include "global_constraints_object.h"
#include "libmesh/libmesh.h"
#include "libmesh/numeric_vector.h"
#include "variables_map.h"

namespace
{
    VectorValue<Real> read_vector(const NumericVector<Number> &vec,
                                  const std::array<dof_id_type, 3> &dofs)
    {
        return VectorValue<Real>(vec(dofs[0]), vec(dofs[1]), vec(dofs[2]));
    }

    void write_vector(NumericVector<Number> &vec,
                      const std::array<dof_id_type, 3> &dofs,
                      const VectorValue<Real> &values)
    {
        vec.set(dofs[0], values(0));
        vec.set(dofs[1], values(1));
        vec.set(dofs[2], values(2));
    }

    template <typename TransformFn>
    void apply_vector_transform(System &system,
                                bool is_nodal,
                                const std::array<unsigned int, 3> &var_numbers,
                                TransformFn &&transform)
    {
        const auto &dof_map = system.get_dof_map();
        std::vector<dof_id_type> dof_indices;
        const auto &mesh = system.get_mesh();

        if (is_nodal)
        {
            for (const auto *node : mesh.local_node_ptr_range())
            {
                std::array<dof_id_type, 3> dofs{};
                bool complete = true;
                for (unsigned int i = 0; i < 3; ++i)
                {
                    dof_map.dof_indices(node, dof_indices, var_numbers[i]);
                    if (dof_indices.empty())
                    {
                        complete = false;
                        break;
                    }
                    dofs[i] = dof_indices.front();
                }
                if (!complete)
                    continue;

                transform(*node, dofs);
            }
        }
        else
        {
            for (const auto &elem : mesh.active_local_element_ptr_range())
            {
                if (elem->infinite() || elem->type() == NODEELEM)
                    continue;
                std::array<dof_id_type, 3> dofs{};
                bool complete = true;
                for (unsigned int i = 0; i < 3; ++i)
                {
                    dof_map.dof_indices(elem, dof_indices, var_numbers[i]);
                    if (dof_indices.empty())
                    {
                        complete = false;
                        break;
                    }
                    dofs[i] = dof_indices.front();
                }
                if (!complete)
                    continue;

                transform(elem->centroid(), dofs);
            }
        }
    }
} // namespace

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

template <typename SysLike>
bool ApecMSystem<SysLike>::global_to_local_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                                      const std::string &variable_label)
{
    std::array<unsigned int, 3> var_numbers{};

    NumericVector<Number> &solution = *this->solution;
    const bool is_nodal = ApecMSystemBase::is_nodal_system(system_type);

    if (!VariablesMap::resolve_vector_var_numbers(*this, variable_label, var_numbers))
        return true;

    auto transform_values = [&](const Point &point,
                                const std::array<dof_id_type, 3> &dofs)
    {
        VectorValue<Real> global = read_vector(solution, dofs);
        VectorValue<Real> local = coordinate_system.global_to_local_rotation(point) * global;
        write_vector(solution, dofs, local);
    };

    apply_vector_transform(*this, is_nodal, var_numbers, transform_values);

    solution.close();
    this->update();
    return true;
}

template <typename SysLike>
bool ApecMSystem<SysLike>::local_to_global_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                                      const std::string &variable_label)
{
    std::array<unsigned int, 3> var_numbers{};

    NumericVector<Number> &solution = *this->solution;
    const bool is_nodal = ApecMSystemBase::is_nodal_system(system_type);

    if (!VariablesMap::resolve_vector_var_numbers(*this, variable_label, var_numbers))
        return true;

    auto transform_values = [&](const Point &point,
                                const std::array<dof_id_type, 3> &dofs)
    {
        VectorValue<Real> local = read_vector(solution, dofs);
        VectorValue<Real> global = coordinate_system.local_to_global_rotation(point) * local;
        write_vector(solution, dofs, global);
    };

    apply_vector_transform(*this, is_nodal, var_numbers, transform_values);

    solution.close();
    this->update();
    return true;
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

template class ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class ApecMSystem<libMesh::ExplicitSystem>;

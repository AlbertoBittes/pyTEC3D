#include "apec_m_system_objects.h"
#include "libmesh/numeric_vector.h"

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
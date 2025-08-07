#include "apec_m_system_objects.h"
#include "libmesh/numeric_vector.h"

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
    NumericVector<Number> *time_derivative_vector = this->request_vector("solution_derivative");
    if (!transient_effect)
    {
        // If transient effects are not considered, the solution time derivative is zero
        time_derivative_vector->zero();
    }
    else
    {
        Real theta = time_substep.get_theta();
        Real theta_x_dt = theta * time_substep.get_dt();
        *time_derivative_vector += *this->current_local_solution;
        *time_derivative_vector -= *this->old_local_solution;
        time_derivative_vector->scale(1 / theta_x_dt);
        if (theta != 1.0)
        {
            auto time_derivative_vector_old = this->request_vector("solution_derivative_old")->clone();
            Real theta_minus_1_over_theta = (theta - 1) / theta;
            time_derivative_vector_old->scale(theta_minus_1_over_theta);
            *time_derivative_vector += *time_derivative_vector_old;
        }
    }
}
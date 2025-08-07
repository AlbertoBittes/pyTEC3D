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
    // NumericVector<Number> *time_derivative_vector = this->request_vector("solution_derivative");
    if (!transient_effect)
    {
        // If transient effects are not considered, the solution time derivative is zero
        this->older_local_solution->zero();
    }
    else
    {
        Real theta = time_substep.get_theta();
        Real theta_x_dt = theta * time_substep.get_dt();
        Real theta_minus_1_dt = (theta - 1) * time_substep.get_dt();
        this->older_local_solution->scale(theta_minus_1_dt);
        *this->older_local_solution += *this->current_local_solution;
        *this->older_local_solution -= *this->old_local_solution;
        this->older_local_solution->scale(1 / theta_x_dt);
    }
}
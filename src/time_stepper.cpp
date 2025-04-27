#include "../include/time_stepper.h"
#include <fstream>

// Method to update the old solution
// This should go before a solve, but in solve
// void LoadStepper::update_old_solution(TransientLinearImplicitSystem &system)
// {
//     *system.old_local_solution = *system.current_local_solution;
// }

void TimeStep::set_substeps(unsigned int n)
{
    Real dt;
    substeps.clear();
    for (size_t i = 0; i < n; ++i)
    {
        if (n <= 0)
        {
            libmesh_error_msg("Error: Number of substeps " + std::to_string(i) + " is not positive.");
        }
        dt = this->get_delta_t() / n;
        for (unsigned int j = 0; j < n; ++j)
        {
            this->add_substep(Substep(j, dt, this));
        }
    }
}

TimeStep::Iterator TimeStep::begin()
{
    return TimeStep::Iterator(this, 0);
}

TimeStep::Iterator TimeStep::end()
{
    return TimeStep::Iterator(this, this->get_n_substeps());
}

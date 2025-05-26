#include "../include/time_step.h"
#include <fstream>
#include "../include/load_step.h" // Include the header for LoadStep

// Method to update the old solution
// This should go before a solve, but in solve
// void LoadStepper::update_old_solution(TransientLinearImplicitSystem &system)
// {
//     *system.old_local_solution = *system.current_local_solution;
// }
void TimeStep::set_beginning_time(std::shared_ptr<LoadStep> prev_load_step)
{
    if (prev_load_step == nullptr)
        time_beginning = 0.0;
    else
        time_beginning = prev_load_step->get_time_step()->get_time();
}
void TimeStep::set_dt()
{
    delta_t = time - time_beginning;
    if (time <= 0)
    {
        libmesh_error_msg("Error: Time step " + std::to_string(time) + " is not positive.");
    }
}
void TimeStep::populate_time_step(std::shared_ptr<LoadStep> prev_load_step)
{
    // Populate attributes from the data_map
    int n_substep;

    step = this->get_setting<int>("time_step_number");
    time = this->get_setting<Real>("time");
    n_substep = this->get_setting<int>("number_of_substeps");
    transient_effect = this->get_setting<int>("transient_effects");
    ramped_load = this->get_setting<int>("ramped_loads");
    write_results = this->get_setting<int>("write_substep_results");

    this->set_beginning_time(prev_load_step);
    // Populate substeps if they exist in the data_map
    this->set_dt();

    this->set_substeps(n_substep);
}
void TimeStep::set_substeps(int n)
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
        this->add_substep(Substep(i, dt, this));
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

#include "../include/time_stepper.h"
#include <fstream>

void TimeStepper::set_cumulative_substep(unsigned int cum_substep)
{
    cumulative_substep = cum_substep;
}

unsigned int TimeStepper::get_cumulative_substep() const
{
    return cumulative_substep;
}
// Method to update the old solution
// This should go before a solve, but in solve
// void TimeStepper::update_old_solution(TransientLinearImplicitSystem &system)
// {
//     *system.old_local_solution = *system.current_local_solution;
// }

TimeDataVector *TimeStepper::read_time_data_from_csv(const std::string &filename)
{
    TimeDataVector *time_data = new TimeDataVector();
    std::ifstream file(filename);

    if (!file.is_open())
    {
        libmesh_error_msg("Error opening time file");
    }

    std::string line;
    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        unsigned int time_step_num;
        Real time;
        unsigned int number_of_substeps;
        bool transient_effects;
        bool ramped_loads;
        bool write_substep_results;
        char comma; // to consume the commas

        ss >> time_step_num >> comma >> time >> comma >> number_of_substeps >> comma >> transient_effects >> comma >> ramped_loads >> comma >> write_substep_results;

        time_data->emplace_back(time_step_num, time, number_of_substeps, transient_effects, ramped_loads, write_substep_results);
    }

    file.close();
    return time_data;
}

void TimeStepper::allocate_time_steps(TimeDataVector *time_data)
{
    for (const auto &[time_step_num, time, number_of_substeps, transient_effects, ramped_loads, write_substep_results] : *time_data)
    {
        TimeStep time_step(time_step_num, time, *this);
        time_step.set_transient_effect(transient_effects);
        time_step.set_ramped_load(ramped_loads);
        time_step.set_write_results(write_substep_results);
        time_step.set_substeps(number_of_substeps);
        add_time_step(time_step);
    }
}

TimeStep::TimeStep(unsigned int step, Real time, TimeStepper &time_stepper)
    : step(step),
      time(time),
      time_beginning(0.0), // Initialize time_beginning before delta_t
      delta_t(0.0),
      write_results(false),
      transient_effect(true),
      ramped_load(true),
      solved(false),             // Initialize solved before time_stepper
      time_stepper(time_stepper) // Initialize time_stepper last
{
    if (step == 0)
        time_beginning = 0.0;
    else
        time_beginning = time_stepper.get_time_steps()[step - 1].get_time();

    if (time <= 0)
    {
        libmesh_error_msg("Error: Time step " + std::to_string(time) + " is not positive.");
    }
    delta_t = time - time_beginning;
    if (delta_t < 0)
    {
        libmesh_error_msg("Error: Time step " + std::to_string(delta_t) + " is not greater than the previous time step.");
    }
}

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
            this->add_substep(Substep(j, dt, *this, this->get_time_stepper()));
        }
    }
}
// Iteration methods
TimeStepper::Iterator TimeStepper::begin()
{
    return TimeStepper::Iterator(*this, 0);
}

TimeStepper::Iterator TimeStepper::end()
{
    return TimeStepper::Iterator(*this, time_steps.size());
}

TimeStep::Iterator TimeStep::begin()
{
    return TimeStep::Iterator(*this, 0);
}

TimeStep::Iterator TimeStep::end()
{
    return TimeStep::Iterator(*this, this->get_n_substeps());
}

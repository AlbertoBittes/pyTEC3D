#include "time_stepper.h"

// Constructor
TimeStepper::TimeStepper(const nlohmann::json &config, System &system)
    : current_time_step(0), current_substep(0), current_time(0.0), system(system)
{
    read_config(config);
    validate_config();
}

// Setters and getters
void TimeStepper::set_current_time_step(unsigned int step)
{
    current_time_step = step;
}

unsigned int TimeStepper::get_current_time_step() const
{
    return current_time_step;
}

void TimeStepper::set_current_substep(unsigned int substep)
{
    current_substep = substep;
}

unsigned int TimeStepper::get_current_substep() const
{
    return current_substep;
}

void TimeStepper::set_current_time(Real time)
{
    current_time = time;
    system.time = time;
}

Real TimeStepper::get_current_time() const
{
    return current_time;
}

void TimeStepper::set_time_step_dt(Real dt)
{
    time_step_dt = dt;
}

Real TimeStepper::get_time_step_dt() const
{
    return time_step_dt;
}

void TimeStepper::set_time_substep_dt(Real dt)
{
    time_substep_dt = dt;
}

Real TimeStepper::get_time_substep_dt() const
{
    return time_substep_dt;
}

// Method to update the old solution
void TimeStepper::update_old_solution()
{
    *system.old_local_solution = *system.current_local_solution;
}

// Iteration methods
TimeStepper::Iterator TimeStepper::begin()
{
    return Iterator(*this, 0);
}

TimeStepper::Iterator TimeStepper::end()
{
    return Iterator(*this, time_steps.size());
}

TimeStepper::SubstepIterator TimeStepper::substep_begin(unsigned int step_index)
{
    return SubstepIterator(*this, step_index, 0);
}

TimeStepper::SubstepIterator TimeStepper::substep_end(unsigned int step_index)
{
    return SubstepIterator(*this, step_index, time_steps[step_index].get_substeps().size());
}

// Helper methods
void TimeStepper::read_config(const nlohmann::json &config)
{
    if (config.contains("Transient Analysis Settings"))
    {
        const auto &settings = config["Transient Analysis Settings"];

        if (settings.contains("time_steps"))
        {
            auto time_steps_json = settings["time_steps"].get<std::vector<Real>>();
            Real previous_time = 0.0;
            for (size_t i = 0; i < time_steps_json.size(); ++i)
            {
                Real delta_t = time_steps_json[i] - previous_time;
                time_steps.emplace_back(i, time_steps_json[i], delta_t, 0, false, false, false);
                previous_time = time_steps_json[i];
            }
        }
        else
        {
            libmesh_error_msg("time_steps not found in config file");
        }

        if (settings.contains("number_of_substeps"))
        {
            auto substeps_json = settings["number_of_substeps"].get<std::vector<unsigned int>>();
            for (size_t i = 0; i < substeps_json.size(); ++i)
            {
                time_steps[i].set_n_substeps(substeps_json[i]);
                Real dt = time_steps[i].get_delta_t() / substeps_json[i];
                for (unsigned int j = 0; j < substeps_json[i]; ++j)
                {
                    time_steps[i].add_substep(Substep(j, time_steps[i].get_time() + j * dt, dt, &time_steps[i]));
                }
            }
        }
        else
        {
            libmesh_error_msg("number_of_substeps not found in config file");
        }

        if (settings.contains("write_substeps_results"))
        {
            auto write_results_json = settings["write_substeps_results"].get<std::vector<bool>>();
            for (size_t i = 0; i < write_results_json.size(); ++i)
            {
                time_steps[i].set_write_results(write_results_json[i]);
            }
        }
        else
        {
            libmesh_error_msg("write_substeps_results not found in config file");
        }

        if (settings.contains("transient_effects"))
        {
            auto transient_effects_json = settings["transient_effects"].get<std::vector<bool>>();
            for (size_t i = 0; i < transient_effects_json.size(); ++i)
            {
                time_steps[i].set_transient_effect(transient_effects_json[i]);
            }
        }
        else
        {
            libmesh_error_msg("transient_effects not found in config file");
        }

        if (settings.contains("ramped_load"))
        {
            auto ramped_load_json = settings["ramped_load"].get<std::vector<bool>>();
            for (size_t i = 0; i < ramped_load_json.size(); ++i)
            {
                time_steps[i].set_ramped_load(ramped_load_json[i]);
            }
        }
        else
        {
            libmesh_error_msg("ramped_load not found in config file");
        }
    }
    else
    {
        libmesh_error_msg("Transient Analysis Settings not found in config file");
    }
}

void TimeStepper::validate_config()
{
    for (const auto &step : time_steps)
    {
        if (step.get_time() <= 0)
        {
            libmesh_error_msg("Error: Time step " + std::to_string(step.get_step()) + " is not positive.");
        }
        if (step.get_step() > 0 && step.get_time() <= time_steps[step.get_step() - 1].get_time())
        {
            libmesh_error_msg("Error: Time step " + std::to_string(step.get_step()) + " is not greater than the previous time step.");
        }
        for (const auto &substep : step.get_substeps())
        {
            if (substep.get_time() <= 0)
            {
                libmesh_error_msg("Error: Substep " + std::to_string(substep.get_index()) + " in time step " + std::to_string(step.get_step()) + " is not positive.");
            }
        }
    }
}
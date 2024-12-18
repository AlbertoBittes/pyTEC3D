#include "../include/time_stepper.h"

bool Substep::get_write_results() const
{
    return parent->get_write_results();
}

bool Substep::get_transient_effect() const
{
    return parent->get_transient_effect();
}

bool Substep::get_ramped_load() const
{
    return parent->get_ramped_load();
}
// Constructor
TimeStepper::TimeStepper(const nlohmann::json &config, TransientLinearImplicitSystem &system)
    : current_time_step(0), current_substep(0), current_time(0.0), cumulative_substep(0), system(system)
{
    read_config(config);
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

void TimeStepper::set_cumulative_substep(unsigned int cum_substep)
{
    cumulative_substep = cum_substep;
}

unsigned int TimeStepper::get_cumulative_substep() const
{
    return cumulative_substep;
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
    return TimeStep::Iterator(*this, n_substeps);
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
            auto write_results_json = settings["write_substeps_results"].get<std::vector<bool>>();
            auto transient_effects_json = settings["transient_effects"].get<std::vector<bool>>();
            auto ramped_load_json = settings["ramped_load"].get<std::vector<bool>>();

            Real previous_time = 0.0;
            for (unsigned int i = 0; i < time_steps_json.size(); ++i)
            {
                if (time_steps_json[i] <= 0)
                {
                    libmesh_error_msg("Error: Time step " + std::to_string(time_steps_json[i]) + " is not positive.");
                }
                Real delta_t = time_steps_json[i] - previous_time;
                if (delta_t < 0)
                {
                    libmesh_error_msg("Error: Time step " + std::to_string(time_steps_json[i]) + " is not greater than the previous time step.");
                }
                time_steps.emplace_back(i, time_steps_json[i], delta_t, 0, write_results_json[i], transient_effects_json[i], ramped_load_json[i], *this);
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
                if (substeps_json[i] <= 0)
                {
                    libmesh_error_msg("Error: Number of substeps " + std::to_string(i) + " is not positive.");
                }
                time_steps[i].set_n_substeps(substeps_json[i]);
                Real dt = time_steps[i].get_delta_t() / substeps_json[i];
                for (unsigned int j = 0; j < substeps_json[i]; ++j)
                {
                    time_steps[i].add_substep(Substep(j, time_steps[i].get_time() + j * dt, dt, &time_steps[i], *this));
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

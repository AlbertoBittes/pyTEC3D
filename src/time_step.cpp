#include <fstream>
#include "time_step.h"
#include "fem_data_repository.h"
#include "fe_parent_object.h"
Substep::Substep(unsigned int index, Real dt, TimeStep *parent)
    : FEParentObjectBase<TimeStep>(parent),
      index(index),
      time(0.0),
      dt(dt),
      solved(false)
{
    time = parent->get_time_beginning() + (index + 1) * dt;
}

bool Substep::get_write_results() const { return this->get_parent_object()->get_write_results(); };
bool Substep::get_transient_effect() const { return this->get_parent_object()->get_transient_effect(); };
bool Substep::get_ramped_load() const { return this->get_parent_object()->get_ramped_load(); };
Real Substep::get_theta() const { return this->get_parent_object()->get_theta(); };

void TimeStep::set_beginning_time(const std::shared_ptr<TimeStep> &previous_time_step)
{
    if (!previous_time_step)
        time_beginning = 0.0;
    else
    {
        time_beginning = previous_time_step->get_time();
    }
}
void TimeStep::set_dt()
{
    delta_t = time - time_beginning;
    if (time <= 0)
    {
        libmesh_error_msg("Error: Time step " + std::to_string(time) + " is not positive.");
    }
}
void TimeStep::populate_time_step(const std::shared_ptr<TimeStep> &previous_time_step)
{
    // Populate attributes from the data_map
    int n_substep;

    step = this->get_int("time_step_number");
    time = this->get_real("time");
    n_substep = this->get_int("number_of_substeps");
    transient_effect = this->get_int("transient_effects");
    ramped_load = this->get_int("ramped_loads");
    theta = this->get_real("theta");
    write_results = this->get_int("write_substep_results");

    this->set_beginning_time(previous_time_step);
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

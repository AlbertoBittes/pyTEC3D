#include <fstream>
#include "time_step.h"
#include "fem_data_repository.h"
#include "fe_parent_object.h"

Substep::Substep(unsigned int index, Real dt, TimeStep *parent)
    : FEParentObjectBase<TimeStep>(parent),
      index(index),
      time(0.0),
      dt(dt),
      solved(false),
      number_picard_iterations(0)
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

void TimeStep::print_data_info() const
{
    if (empty() || this->get_data_map()->columns.empty())
    {
        libMesh::out << "  Time step: no data loaded." << std::endl;
        return;
    }

    libMesh::out << "  Step " << step << " at time " << time << " (start " << time_beginning << ", delta_t=" << delta_t << ")" << std::endl;
    libMesh::out << "  Transient effect: " << (transient_effect ? "on" : "off")
                 << ", ramped loads: " << (ramped_load ? "on" : "off")
                 << ", theta=" << theta
                 << ", write substep results: " << (write_results ? "yes" : "no") << std::endl;

    libMesh::out << "  Substeps: " << substeps.size();
    if (!substeps.empty())
    {
        libMesh::out << " (dt per substep ~ ";
        Real dt_per = substeps.front().get_dt();
        libMesh::out << dt_per << ")";
        libMesh::out << "\n    Times:";
        for (const auto &substep : substeps)
            libMesh::out << " " << substep.get_time();
    }
    libMesh::out << std::endl;
}
void TimeStep::populate_time_step(const std::shared_ptr<TimeStep> &previous_time_step)
{
    substeps.clear();
    active_substep = nullptr;
    previous_substep = nullptr;
    previous_time_step_ptr = previous_time_step;

    if (empty() || this->get_data_map()->columns.empty())
    {
        step = 0;
        time = 0.0;
        time_beginning = previous_time_step ? previous_time_step->get_time() : 0.0;
        delta_t = 0.0;
        write_results = 0;
        transient_effect = 0;
        ramped_load = 0;
        theta = 1.0;
        return;
    }

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

std::size_t TimeStep::count_steps_including_previous() const
{
    std::size_t count = 0;
    const TimeStep *ptr = this;
    while (ptr)
    {
        ++count;
        ptr = ptr->previous_time_step_ptr.get();
    }
    return count;
}

std::size_t TimeStep::count_substeps_before_current() const
{
    std::size_t total = 0;
    const TimeStep *ptr = previous_time_step_ptr.get();
    while (ptr)
    {
        total += ptr->get_n_substeps();
        ptr = ptr->previous_time_step_ptr.get();
    }
    return total;
}

std::size_t TimeStep::count_total_substeps_including_current() const
{
    return count_substeps_before_current() + static_cast<std::size_t>(get_n_substeps());
}

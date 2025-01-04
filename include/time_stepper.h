#ifndef TIME_STEPPER_H
#define TIME_STEPPER_H

#include <stdexcept>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include "../include/json.hpp"
#include <libmesh/libmesh.h>
#include <libmesh/system.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include "libmesh/numeric_vector.h"
#include <libmesh/sparse_matrix.h> // Add this line to include the SparseMatrix definition

using namespace libMesh;

class TimeStep;
class TimeStepper;

class Substep
{
public:
    Substep(unsigned int index, Real time, Real dt, TimeStep *parent, TimeStepper &time_stepper)
        : index(index), time(time), dt(dt), parent(parent), time_stepper(time_stepper) {}

    unsigned int get_index() const { return index; }
    Real get_time() const { return time; }
    Real get_dt() const { return dt; }
    bool get_write_results() const;
    bool get_transient_effect() const;

    void set_index(unsigned int i) { index = i; }
    void set_time(Real t) { time = t; }
    void set_dt(Real d) { dt = d; }

private:
    unsigned int index;
    Real time;
    Real dt;
    TimeStep *parent;
    TimeStepper &time_stepper;
};

class TimeStep
{
public:
    TimeStep(unsigned int step, Real time, Real delta_t, unsigned int n_substeps, bool write_results, bool transient_effect, TimeStepper &time_stepper)
        : step(step), time(time), delta_t(delta_t), n_substeps(n_substeps), write_results(write_results), transient_effect(transient_effect), time_stepper(time_stepper) {}

    unsigned int get_step() const { return step; }
    Real get_time() const { return time; }
    Real get_delta_t() const { return delta_t; }
    unsigned int get_n_substeps() const { return n_substeps; }
    bool get_write_results() const { return write_results; }
    bool get_transient_effect() const { return transient_effect; }
    const std::vector<Substep> &get_substeps() const { return substeps; }
    TimeStepper &get_time_stepper() { return time_stepper; }

    void set_step(unsigned int s) { step = s; }
    void set_time(Real t) { time = t; }
    void set_delta_t(Real dt) { delta_t = dt; }
    void set_n_substeps(unsigned int n) { n_substeps = n; }
    void set_write_results(bool w) { write_results = w; }
    void set_transient_effect(bool t) { transient_effect = t; }
    void add_substep(const Substep &substep) { substeps.push_back(substep); }

    // Custom substep iterator for TimeStep

    class Iterator;
    Iterator begin();
    Iterator end();

private:
    unsigned int step;
    Real time;
    Real delta_t;
    unsigned int n_substeps;
    bool write_results;
    bool transient_effect;
    std::vector<Substep> substeps;
    TimeStepper &time_stepper;
};

class TimeStepper
{
public:
    // Constructor
    TimeStepper(const nlohmann::json &config, TransientLinearImplicitSystem &system);

    // Setters and getters
    void set_current_time_step(unsigned int step);
    unsigned int get_current_time_step() const;

    void set_current_substep(unsigned int substep);
    unsigned int get_current_substep() const;

    void set_cumulative_substep(unsigned int cum_substep);
    unsigned int get_cumulative_substep() const;

    void set_current_time(Real time);
    Real get_current_time() const;

    void set_time_step_dt(Real dt);
    Real get_time_step_dt() const;

    void set_time_substep_dt(Real dt);
    Real get_time_substep_dt() const;

    void set_prev_time_substep_dt(Real dt);
    Real get_prev_time_substep_dt() const;

    const std::vector<TimeStep> &get_time_steps() const { return time_steps; }
    TransientLinearImplicitSystem &get_system() { return system; }
    // Method to update the old solution
    void update_old_solution();

    // Iteration methods
    class Iterator;
    Iterator begin();
    Iterator end();

private:
    unsigned int current_time_step;
    unsigned int current_substep;
    unsigned int cumulative_substep;
    Real current_time;
    Real time_step_dt;
    Real time_substep_dt;
    Real prev_time_substep_dt;
    std::vector<TimeStep> time_steps;
    TransientLinearImplicitSystem &system;

    // Helper methods
    void read_config(const nlohmann::json &config);
};

// Custom iterator for TimeStepper
class TimeStepper::Iterator
{
public:
    Iterator(TimeStepper &ts, unsigned int index) : ts(ts), index(index), first(true) {}

    bool operator!=(const Iterator &other) const { return index != other.index; }
    TimeStep &operator*()
    {
        if (first)
        {
            ts.set_current_time_step(index);
            ts.set_time_step_dt(ts.time_steps[index].get_delta_t());
            ts.get_system().get_system_matrix().zero();
            first = false;
        }
        return ts.time_steps[index];
    }
    Iterator &operator++()
    {
        first = true;
        index++;
        return *this;
    }

private:
    TimeStepper &ts;
    unsigned int index;
    bool first;
};

class TimeStep::Iterator
{
public:
    Iterator(TimeStep &ts, unsigned int substep_index)
        : ts(ts), time_stepper(ts.get_time_stepper()), substep_index(substep_index), first(true) {}

    bool operator!=(const Iterator &other) const { return substep_index != other.substep_index; }
    const Substep &operator*()
    {
        if (first)
        {
            time_stepper.update_old_solution();
            time_stepper.set_current_time(time_stepper.get_current_time() + ts.get_substeps()[substep_index].get_dt());
            time_stepper.get_system().time = time_stepper.get_current_time();
            time_stepper.set_current_substep(substep_index);
            time_stepper.set_time_substep_dt(ts.get_substeps()[substep_index].get_dt());
            if (substep_index > 0)
            {
                time_stepper.set_prev_time_substep_dt(ts.get_substeps()[substep_index - 1].get_dt());
            }
            else
            {
                if (ts.get_step() > 0)
                {
                    time_stepper.set_prev_time_substep_dt(time_stepper.get_time_steps()[ts.get_step() - 1].get_substeps().back().get_dt());
                }
                else
                {
                    time_stepper.set_prev_time_substep_dt(0);
                }
            }
            time_stepper.set_cumulative_substep(time_stepper.get_cumulative_substep() + 1);
            time_stepper.get_system().rhs->zero();
            first = false;
        }
        return ts.get_substeps()[substep_index];
    }
    Iterator &operator++()
    {
        first = true;
        substep_index++;
        return *this;
    }

private:
    TimeStep &ts;
    TimeStepper &time_stepper;
    unsigned int substep_index;
    bool first;
};

#endif // TIME_STEPPER_H
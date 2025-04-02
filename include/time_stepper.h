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
    Substep(unsigned int index, Real dt, TimeStep &parent, TimeStepper &time_stepper)
        : index(index), dt(dt), parent(parent), time_stepper(time_stepper)
    {
        time = parent.get_time_beginning() + (index + 1) * dt;
    }

    unsigned int get_index() const { return index; }
    Real get_time() const { return time; }
    Real get_dt() const { return dt; }
    bool get_write_results() const { return parent.get_write_results(); };
    bool get_transient_effect() const { return parent.get_transient_effect(); };
    bool get_ramped_load() const { return parent.get_ramped_load(); };
    TimeStep &get_parent() { return parent; }

    void set_index(unsigned int i) { index = i; }
    void set_time(Real t) { time = t; }
    void set_dt(Real d) { dt = d; }

private:
    unsigned int index;
    Real time;
    Real dt;
    TimeStep &parent;
    TimeStepper &time_stepper;
};

class TimeStep
{
public:
    TimeStep(unsigned int step, Real time, TimeStepper &time_stepper)
        : time_stepper(time_stepper) {}

    unsigned int get_step() const { return step; }
    Real get_time() const { return time; }
    Real get_time_beginning() const { return time_beginning; }
    Real get_delta_t() const { return delta_t; }
    unsigned int get_n_substeps() const { return substeps.size(); }
    bool get_write_results() const { return write_results; }
    bool get_transient_effect() const { return transient_effect; }
    bool get_ramped_load() const { return ramped_load; }
    std::vector<Substep> &get_substeps() { return substeps; }
    const std::vector<Substep> &get_substeps() const { return substeps; }
    TimeStepper &get_time_stepper() { return time_stepper; }

    void set_step(unsigned int s) { step = s; }
    void set_time(Real t) { time = t; }
    void set_time_beginning(Real t) { time_beginning = t; }
    void set_delta_t(Real dt) { delta_t = dt; }
    // void set_n_substeps(unsigned int n) { n_substeps = n; }
    void set_write_results(bool w) { write_results = w; }
    void set_transient_effect(bool t) { transient_effect = t; }
    void set_ramped_load(bool r) { ramped_load = r; }
    void add_substep(const Substep &substep) { substeps.push_back(substep); }
    void set_substeps(unsigned int);
    // Custom substep iterator for TimeStep

    class Iterator;
    Iterator begin();
    Iterator end();

private:
    unsigned int step;
    Real time;
    Real time_beginning;
    Real delta_t;
    bool write_results;
    bool transient_effect;
    bool ramped_load;
    std::vector<Substep> substeps;
    TimeStepper &time_stepper;
};

typedef std::vector<std::tuple<unsigned int, Real, unsigned int, bool, bool, bool>> TimeDataVector;

class TimeStepper
{
public:
    // Constructor
    TimeStepper() {};

    void set_cumulative_substep(unsigned int cum_substep);
    unsigned int get_cumulative_substep() const;

    const std::vector<TimeStep> &get_time_steps() const { return time_steps; }
    TimeStep &get_current_time_step() { return time_steps.back(); }
    // Method to update the old solution
    // Has to go in solve
    void update_old_solution(TransientLinearImplicitSystem &system);

    void add_time_step(const TimeStep &time_step) { time_steps.push_back(time_step); }
    TimeDataVector &read_time_data_from_csv(const std::string &filename);
    void allocate_time_steps(TimeDataVector &time_data);

    // Iteration methods
    class Iterator;
    Iterator begin();
    Iterator end();

private:
    std::vector<TimeStep> time_steps;
    unsigned int active_time_step;
    unsigned int current_substep;
    unsigned int cumulative_substep;
    // Helper methods
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
            first = false;
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
    Iterator(TimeStep &time_step, unsigned int substep_index)
        : time_step(time_step), time_stepper(time_step.get_time_stepper()), substep_index(substep_index), first(true) {}

    bool operator!=(const Iterator &other) const { return substep_index != other.substep_index; }
    const Substep &operator*()
    {
        if (first)
        {
            time_stepper.set_cumulative_substep(time_stepper.get_cumulative_substep() + 1);
            // This should go after a solve
            // time_stepper.get_system().rhs->zero();
            first = false;
        }
        return time_step.get_substeps()[substep_index];
    }
    Iterator &operator++()
    {
        first = true;
        substep_index++;
        return *this;
    }

private:
    TimeStep &time_step;
    TimeStepper &time_stepper;
    unsigned int substep_index;
    bool first;
};

#endif // TIME_STEPPER_H
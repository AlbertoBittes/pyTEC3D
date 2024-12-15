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

using namespace libMesh;

class TimeStep;

class Substep
{
public:
    Substep(unsigned int index, Real time, Real dt, TimeStep *parent)
        : index(index), time(time), dt(dt), parent(parent) {}

    unsigned int get_index() const { return index; }
    Real get_time() const { return time; }
    Real get_dt() const { return dt; }
    bool get_write_results() const { return parent->get_write_results(); }
    bool get_transient_effect() const { return parent->get_transient_effect(); }
    bool get_ramped_load() const { return parent->get_ramped_load(); }

    void set_index(unsigned int i) { index = i; }
    void set_time(Real t) { time = t; }
    void set_dt(Real d) { dt = d; }

private:
    unsigned int index;
    Real time;
    Real dt;
    TimeStep *parent;
};

class TimeStep
{
public:
    TimeStep(unsigned int step, Real time, Real delta_t, unsigned int n_substeps, bool write_results, bool transient_effect, bool ramped_load)
        : step(step), time(time), delta_t(delta_t), n_substeps(n_substeps), write_results(write_results), transient_effect(transient_effect), ramped_load(ramped_load) {}

    unsigned int get_step() const { return step; }
    Real get_time() const { return time; }
    Real get_delta_t() const { return delta_t; }
    unsigned int get_n_substeps() const { return n_substeps; }
    bool get_write_results() const { return write_results; }
    bool get_transient_effect() const { return transient_effect; }
    bool get_ramped_load() const { return ramped_load; }
    const std::vector<Substep> &get_substeps() const { return substeps; }

    void set_step(unsigned int s) { step = s; }
    void set_time(Real t) { time = t; }
    void set_delta_t(Real dt) { delta_t = dt; }
    void set_n_substeps(unsigned int n) { n_substeps = n; }
    void set_write_results(bool w) { write_results = w; }
    void set_transient_effect(bool t) { transient_effect = t; }
    void set_ramped_load(bool r) { ramped_load = r; }
    void add_substep(const Substep &substep) { substeps.push_back(substep); }

    // Iteration methods for substeps
    std::vector<Substep>::iterator begin() { return substeps.begin(); }
    std::vector<Substep>::iterator end() { return substeps.end(); }

private:
    unsigned int step;
    Real time;
    Real delta_t;
    unsigned int n_substeps;
    bool write_results;
    bool transient_effect;
    bool ramped_load;
    std::vector<Substep> substeps;
};

class TimeStepper
{
public:
    // Constructor
    TimeStepper(const nlohmann::json &config, System &system);

    // Setters and getters
    void set_current_time_step(unsigned int step);
    unsigned int get_current_time_step() const;

    void set_current_substep(unsigned int substep);
    unsigned int get_current_substep() const;

    void set_current_time(Real time);
    Real get_current_time() const;

    void set_time_step_dt(Real dt);
    Real get_time_step_dt() const;

    void set_time_substep_dt(Real dt);
    Real get_time_substep_dt() const;

    // Method to update the old solution
    void update_old_solution();

    // Iteration methods
    class Iterator;
    Iterator begin();
    Iterator end();

    class SubstepIterator;
    SubstepIterator substep_begin(unsigned int step_index);
    SubstepIterator substep_end(unsigned int step_index);

private:
    unsigned int current_time_step;
    unsigned int current_substep;
    Real current_time;
    Real time_step_dt;
    Real time_substep_dt;
    std::vector<TimeStep> time_steps;
    System &system;

    // Helper methods
    void read_config(const nlohmann::json &config);
    void validate_config();
};

// Custom iterator for TimeStepper
class TimeStepper::Iterator
{
public:
    Iterator(TimeStepper &ts, size_t index) : ts(ts), index(index) {}

    bool operator!=(const Iterator &other) const { return index != other.index; }
    const TimeStep &operator*() { return ts.time_steps[index]; }
    Iterator &operator++()
    {
        index++;
        return *this;
    }

private:
    TimeStepper &ts;
    size_t index;
};

// Custom substep iterator for TimeStepper
class TimeStepper::SubstepIterator
{
public:
    SubstepIterator(TimeStepper &ts, size_t step_index, size_t substep_index)
        : ts(ts), step_index(step_index), substep_index(substep_index), first(true) {}

    bool operator!=(const SubstepIterator &other) const { return substep_index != other.substep_index; }
    const Substep &operator*()
    {
        if (first)
        {
            ts.update_old_solution();
            ts.current_time += ts.time_steps[step_index].get_substeps()[substep_index].get_dt();
            ts.system.time = ts.current_time;
            ts.current_substep++;
            first = false;
        }
        return ts.time_steps[step_index].get_substeps()[substep_index];
    }
    SubstepIterator &operator++()
    {
        first = true;
        substep_index++;
        return *this;
    }

private:
    TimeStepper &ts;
    size_t step_index;
    size_t substep_index;
    bool first;
};

#endif // TIME_STEPPER_H
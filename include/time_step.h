#ifndef TIME_STEP_H
#define TIME_STEP_H

#include <stdexcept>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <libmesh/libmesh.h>
#include "../include/load_step.h"

using namespace libMesh;

class Substep;

class TimeStep
{
public:
    // Constructor with binary file stream
    TimeStep(std::istream &file, LoadStep *load_step, LoadStep *prev_load_step = nullptr)
        : load_step(load_step)
    {
        // Read step
        file.read(reinterpret_cast<char *>(&step), sizeof(step));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read 'step' from the binary stream.");
        }

        // Read time
        file.read(reinterpret_cast<char *>(&time), sizeof(time));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read 'time' from the binary stream.");
        }

        // Read write_results as unsigned int
        unsigned int write_results_uint;
        file.read(reinterpret_cast<char *>(&write_results_uint), sizeof(write_results_uint));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read 'write_results' from the binary stream.");
        }
        bool write_results = (write_results_uint != 0);

        // Read transient_effect as unsigned int
        unsigned int transient_effect_uint;
        file.read(reinterpret_cast<char *>(&transient_effect_uint), sizeof(transient_effect_uint));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read 'transient_effect' from the binary stream.");
        }
        bool transient_effect = (transient_effect_uint != 0);

        // Read ramped_load as unsigned int
        unsigned int ramped_load_uint;
        file.read(reinterpret_cast<char *>(&ramped_load_uint), sizeof(ramped_load_uint));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read 'ramped_load' from the binary stream.");
        }
        bool ramped_load = (ramped_load_uint != 0);

        // Read the number of substeps
        unsigned int num_substeps;
        file.read(reinterpret_cast<char *>(&num_substeps), sizeof(num_substeps));
        if (!file)
        {
            libmesh_error_msg("Error: Failed to read the number of substeps from the binary stream.");
        }

        *this = TimeStep(step, time, load_step, write_results, num_substeps,
                         transient_effect, ramped_load, prev_load_step);
    }

    TimeStep(unsigned int step, Real time, LoadStep *load_step, unsigned int num_substeps = 1, bool write_results = true,
             bool transient_effect = true, bool ramped_load = true, LoadStep *prev_load_step = nullptr)
        : step(step),
          time(time),
          write_results(write_results),
          transient_effect(transient_effect),
          ramped_load(ramped_load),
          load_step(load_step) // Initialize time_stepper last
    {
        if (prev_load_step == nullptr)
            time_beginning = 0.0;
        else
            time_beginning = prev_load_step->get_time_step()->get_time();

        if (time <= 0)
        {
            libmesh_error_msg("Error: Time step " + std::to_string(time) + " is not positive.");
        }

        delta_t = time - time_beginning;

        if (delta_t < 0)
        {
            libmesh_error_msg("Error: Time step " + std::to_string(delta_t) + " is not greater than the previous time step.");
        }

        this->set_substeps(num_substeps);
    }

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
    LoadStep &get_load_step() { return *load_step; }

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

    void set_beginning_time(LoadStep *load_step)
    {
        if (load_step == nullptr)
            time_beginning = 0.0;
        else
            time_beginning = load_step->get_time_step()->get_time();
    }
    void set_dt(LoadStep *load_step)
    {
        if (load_step == nullptr)
            delta_t = 0.0;
        else
            delta_t = load_step->get_time_step()->get_time() - time_beginning;
    }

private:
    unsigned int step;
    Real time;
    Real time_beginning;
    Real delta_t;
    bool write_results;
    bool transient_effect;
    bool ramped_load;

    std::vector<Substep> substeps;
    LoadStep *load_step;
};

class Substep
{
public:
    Substep(unsigned int index, Real dt, TimeStep *parent)
        : index(index), dt(dt), parent(parent), solved(false)
    {
        time = parent->get_time_beginning() + (index + 1) * dt;
    }

    unsigned int get_index() const { return index; }
    Real get_time() const { return time; }
    Real get_dt() const { return dt; }
    bool get_solved() const { return solved; }

    bool get_write_results() const { return parent->get_write_results(); };
    bool get_transient_effect() const { return parent->get_transient_effect(); };
    bool get_ramped_load() const { return parent->get_ramped_load(); };

    const TimeStep &get_parent() const { return *parent; }
    TimeStep &get_parent() { return *parent; }

    void set_index(unsigned int i) { index = i; }
    void set_time(Real t) { time = t; }
    void set_dt(Real d) { dt = d; }
    void set_solved(bool s) { solved = s; }

private:
    unsigned int index;
    Real time;
    Real dt;
    TimeStep *parent;
    bool solved;
};

class TimeStep::Iterator
{
public:
    Iterator(TimeStep *time_step, unsigned int substep_index)
        : time_step(time_step), substep_index(substep_index), first(true) {}

    bool operator!=(const Iterator &other) const { return substep_index != other.substep_index; }
    const Substep &operator*()
    {
        if (first)
        {
            // This should go after a solve
            // time_stepper.get_system().rhs->zero();
            first = false;
        }
        return time_step->get_substeps()[substep_index];
    }
    Iterator &operator++()
    {
        first = true;
        substep_index++;
        return *this;
    }

private:
    TimeStep *time_step;
    unsigned int substep_index;
    bool first;
};

#endif // TIME_STEP_H
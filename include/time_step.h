#ifndef TIME_STEP_H
#define TIME_STEP_H

#include <memory>
#include <stdexcept>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <libmesh/libmesh.h>
#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "globals.h"

using namespace libMesh;

// The TimeStep object collects all information regarding the transient settings of the analysis:
// Time step size
// Time substeps
// Integration scheme parameter theta
// Ramping of loads, to be used together with time substeps
// Transient effects: when OFF, a static analysis load step is made

class TimeStep; // Forward declaration

class Substep : public FEParentObjectBase<TimeStep>
{
public:
    Substep(unsigned int index, Real dt, TimeStep *parent);

    unsigned int get_index() const { return index; }
    Real get_time() const { return time; }
    Real get_dt() const { return dt; }
    bool get_solved() const { return solved; }

    bool get_write_results() const;
    bool get_transient_effect() const;
    bool get_ramped_load() const;
    Real get_theta() const; // Not used in this context, but kept for compatibility

    void set_index(unsigned int i) { index = i; }
    void set_time(Real t) { time = t; }
    void set_dt(Real d) { dt = d; }
    void set_solved(bool s) { solved = s; }

    void print_solve_msg() const
    {
        libMesh::out << "Substep: " << get_index() + 1 << std::endl
                     << "Time: " << get_time() << std::endl
                     << "dt: " << get_dt() << std::endl
                     << "Write Results: " << get_write_results() << std::endl
                     << "Transient Effect: " << get_transient_effect() << std::endl
                     << "Ramped Load: " << get_ramped_load() << std::endl
                     << "Theta: " << get_theta() << std::endl;
    }

private:
    unsigned int index;
    Real time;
    Real dt;
    bool solved;
};

class TimeStep : public BinaryDataFEDataObject
{
public:
    TimeStep(std::ifstream &file, std::shared_ptr<TimeStep> previous_time_step = nullptr)
        : BinaryDataFEDataObject(file)
    {
        if (previous_time_step)
            this->populate_time_step(previous_time_step);
        else
            this->populate_time_step(nullptr);
    }
    // Constructor with binary file stream
    // TimeStep(std::ifstream &file, std::shared_ptr<TimeStep> previous_time_step = nullptr)
    //     : BinaryDataFEDataObject(file)
    // {
    //     this->populate_time_step(previous_time_step);
    // }

    int get_step() const { return step; }
    Real get_time() const { return time; }
    Real get_time_beginning() const { return time_beginning; }
    Real get_delta_t() const { return delta_t; }
    unsigned int get_n_substeps() const { return substeps.size(); }
    int get_write_results() const { return write_results; }
    int get_transient_effect() const { return transient_effect; }
    int get_ramped_load() const { return ramped_load; }
    Real get_theta() const { return theta; } // Not used in this context, but kept for compatibility
    Substep *get_active_substep() { return active_substep; }
    Substep *get_previous_substep() { return previous_substep; }
    std::vector<Substep> &get_substeps() { return substeps; }
    const std::vector<Substep> &get_substeps() const { return substeps; }

    void set_step(int s) { step = s; }
    void set_time(Real t) { time = t; }
    void set_time_beginning(Real t) { time_beginning = t; }
    void set_delta_t(Real dt) { delta_t = dt; }
    // void set_n_substeps(unsigned int n) { n_substeps = n; }
    void set_write_results(int w) { write_results = w; }
    void set_transient_effect(int t) { transient_effect = t; }
    void set_ramped_load(int r) { ramped_load = r; }
    void set_theta(Real th) { theta = th; } // Not used in this context, but kept for compatibility
    void set_active_substep(Substep *substep) { active_substep = substep; }
    void set_previous_substep(Substep *substep) { previous_substep = substep; }
    void add_substep(const Substep &substep) { substeps.push_back(substep); }
    void set_substeps(int);
    // Custom substep iterator for TimeStep

    void populate_time_step(const std::shared_ptr<TimeStep> &previous_time_step);

    void print_solve_msg() const
    {
        libMesh::out << "Solving time step " << get_step() << std::endl;
        libMesh::out << "Time increment at the end of time step: " << get_delta_t() << std::endl;
        libMesh::out << "Number of time substeps: " << get_n_substeps() << std::endl;
    }
    class Iterator;
    Iterator begin();
    Iterator end();

    void set_beginning_time(const std::shared_ptr<TimeStep> &previous_time_step);

    void set_dt();

private:
    int step;
    Real time;
    Real time_beginning;
    Real delta_t;
    int write_results;
    int transient_effect;
    int ramped_load;
    Real theta; // Not used in this context, but kept for compatibility

    std::vector<Substep> substeps;
    Substep *active_substep = nullptr;
    Substep *previous_substep = nullptr;
};

// The iterator takes care of solving the time substeps sequentially
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
            // Update active_substep and previous_substep
            time_step->set_previous_substep(time_step->get_active_substep());
            time_step->set_active_substep(&time_step->get_substeps()[substep_index]);

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

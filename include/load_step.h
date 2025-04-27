#ifndef LOAD_STEP_H
#define LOAD_STEP_H

#include <memory>
#include "../include/materials_object.h"
#include "../include/loads_object.h"
#include "../include/node_bc_object.h"
#include "../include/time_step.h"
#include "../include/solver_configuration.h"

class LoadStep
{
public:
    // Default constructor: initializes everything to nullptr
    LoadStep()
        : materials(nullptr), loads(nullptr), node_bc(nullptr), time_step(nullptr), solver_config(nullptr) {}

    // Constructor that calls a reading function
    LoadStep(const std::string &filename)
        : materials(nullptr), loads(nullptr), node_bc(nullptr), time_step(nullptr), solver_config(nullptr)
    {
        read_from_file(filename);
    }

    // Dummy reading function
    void read_from_file(const std::string &filename)
    {
        // Placeholder for actual file reading logic
        std::cout << "Reading LoadStep data from file: " << filename << std::endl;

        // Example: Initialize objects (replace with actual logic)
        materials = std::make_shared<MaterialsObject>();
        loads = std::make_shared<LoadsObject>();
        node_bc = std::make_shared<NodeBCObject>();
        time_step = std::make_shared<TimeStep>(0, 0.0, *this); // Example initialization
        solver_config = std::make_shared<SolverConfiguration>();
    }

    // Getters
    std::shared_ptr<MaterialsObject> get_materials() const { return materials; }
    std::shared_ptr<LoadsObject> get_loads() const { return loads; }
    std::shared_ptr<NodeBCObject> get_node_bc() const { return node_bc; }
    std::shared_ptr<TimeStep> get_time_step() const { return time_step; }
    std::shared_ptr<SolverConfiguration> get_solver_config() const { return solver_config; }

    // Setters
    void set_materials(std::shared_ptr<MaterialsObject> materials) { this->materials = materials; }
    void set_loads(std::shared_ptr<LoadsObject> loads) { this->loads = loads; }
    void set_node_bc(std::shared_ptr<NodeBCObject> node_bc) { this->node_bc = node_bc; }
    void set_time_step(std::shared_ptr<TimeStep> time_step) { this->time_step = time_step; }
    void set_solver_config(std::shared_ptr<SolverConfiguration> solver_config) { this->solver_config = solver_config; }

private:
    std::shared_ptr<MaterialsObject> materials;
    std::shared_ptr<LoadsObject> loads;
    std::shared_ptr<NodeBCObject> node_bc;
    std::shared_ptr<TimeStep> time_step;
    std::shared_ptr<SolverConfiguration> solver_config;

    LoadStepper *load_stepper; // Pointer to the LoadStepper (if needed for TimeStep initialization)
};

class LoadStepper
{
public:
    // Constructor
    LoadStepper()
        : active_load_step(nullptr), old_load_step(nullptr) {}

    // Getters
    std::shared_ptr<LoadStep> get_active_load_step() const { return active_load_step; }
    std::shared_ptr<LoadStep> get_old_load_step() const { return old_load_step; }

    // Setters
    void set_active_load_step(std::shared_ptr<LoadStep> load_step) { active_load_step = load_step; }
    void set_old_load_step(std::shared_ptr<LoadStep> load_step) { old_load_step = load_step; }

    // Switch pointers and read a new LoadStep
    void read_new_load_step(const std::string &filename)
    {
        // Move the current active load step to the old load step
        old_load_step = active_load_step;

        // Create a new LoadStep and assign it to the active load step
        active_load_step = std::make_shared<LoadStep>(filename);
    }

    // Method to update the old solution (if needed)
    // void update_old_solution(TransientLinearImplicitSystem &system)
    // {
    //     if (old_load_step)
    //     {
    //         // Placeholder for logic to update the old solution
    //         std::cout << "Updating old solution for the system." << std::endl;
    //     }
    // }

private:
    std::shared_ptr<LoadStep> active_load_step;
    std::shared_ptr<LoadStep> old_load_step;
};

#endif // LOAD_STEP_H
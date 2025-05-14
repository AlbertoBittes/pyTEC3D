#ifndef LOAD_STEP_H
#define LOAD_STEP_H

#include <memory>
#include "../include/fe_data_object.h"
#include "../include/materials_object.h"
#include "../include/loads_object.h"
#include "../include/real_constants_object.h"
#include "../include/node_bc_object.h"
#include "../include/time_step.h"
#include "../include/solver_configuration.h"
#include "libmesh/equation_systems.h" // Include the header for EquationSystems
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/mesh.h"

class LoadStepper;

class LoadStep
{
public:
    // Default constructor: initializes everything to nullptr
    LoadStep(LoadStepper *load_stepper = nullptr, EquationSystems *equation_systems = nullptr)
        : load_stepper(load_stepper),
          time_step(nullptr),
          solver_config(nullptr),
          materials(nullptr),
          real_constants(nullptr),
          loads(nullptr),
          node_bc(nullptr),
          equation_systems(equation_systems)
    {
    }

    void read_from_file(const std::string &filename);

    // Getters
    std::shared_ptr<TimeStep> get_time_step() const { return time_step; }
    std::shared_ptr<PetscSolverConfiguration> get_solver_config() const { return solver_config; }
    std::shared_ptr<MaterialsObject> get_materials() const { return materials; }
    std::shared_ptr<RealConstantsObject> get_real_constants() const { return real_constants; }
    std::shared_ptr<LoadsObject> get_loads() const { return loads; }
    std::shared_ptr<NodeBCObject> get_node_bc() const { return node_bc; }
    LoadStepper *get_load_stepper() const { return load_stepper; }
    EquationSystems *get_equation_systems() const { return equation_systems; }

    // Setters
    void set_time_step(std::shared_ptr<TimeStep> time_step) { this->time_step = time_step; }
    void set_solver_config(std::shared_ptr<PetscSolverConfiguration> solver_config) { this->solver_config = solver_config; }
    void set_materials(std::shared_ptr<MaterialsObject> materials) { this->materials = materials; }
    void set_real_constants(std::shared_ptr<RealConstantsObject> real_constants) { this->real_constants = real_constants; }
    void set_loads(std::shared_ptr<LoadsObject> loads) { this->loads = loads; }
    void set_node_bc(std::shared_ptr<NodeBCObject> node_bc) { this->node_bc = node_bc; }
    void set_load_stepper(LoadStepper *load_stepper) { this->load_stepper = load_stepper; }
    void set_equation_systems(EquationSystems *equation_systems) { this->equation_systems = equation_systems; }

private:
    std::shared_ptr<TimeStep> time_step;                     // Pointer to the TimeStep
    std::shared_ptr<PetscSolverConfiguration> solver_config; // Pointer to the SolverConfiguration
    std::shared_ptr<MaterialsObject> materials;              // Pointer to the MaterialsObject
    std::shared_ptr<RealConstantsObject> real_constants;     // Pointer to the RealConstantsObject
    std::shared_ptr<LoadsObject> loads;                      // Pointer to the LoadsObject
    std::shared_ptr<NodeBCObject> node_bc;                   // Pointer to the NodeBCObject

    LoadStepper *load_stepper;         // Pointer to the LoadStepper (if needed for TimeStep initialization)
    EquationSystems *equation_systems; // Pointer to the EquationSystems (if needed)
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
    void read_new_load_step(LoadStepper *load_stepper, EquationSystems *equation_systems, const std::string &filename)
    {
        // Move the current active load step to the old load step
        old_load_step = active_load_step;

        // Create a new LoadStep and assign it to the active load step
        active_load_step = std::make_shared<LoadStep>(load_stepper, equation_systems);
        active_load_step->read_from_file(filename);

        // Update the load step index
        load_step_index++;
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
    unsigned int load_step_index = 0;
};

#endif // LOAD_STEP_H
#ifndef LOAD_STEP_H
#define LOAD_STEP_H

#include <memory>
#include "fe_data_object.h"
#include "materials_object.h"
#include "loads_object.h"
#include "real_constants_object.h"
#include "node_bc_object.h"
#include "time_step.h"
#include "solver_configuration.h"
#include "libmesh/equation_systems.h" // Include the header for EquationSystems
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/mesh.h"

class LoadStepper;

enum LoadStepObjectType : uint32_t
{
    TIME_STEP = 0,
    PETSC_SOLVER_CONFIG,
    MATERIALS_OBJECT,
    REAL_CONSTANTS_OBJECT,
    LOADS_OBJECT,
    NODE_BC_OBJECT
};

class LoadStep : public FEDataObjectBase<LoadStepper>, public FEMDataContext
{
public:
    // Default constructor
    LoadStep(LoadStepper *load_stepper = nullptr)
        : FEDataObjectBase<LoadStepper>(load_stepper)
    {
    }

    void read_from_file(const std::string &filename);
    // // Add a FEDataObject to the array
    // void add_data_object(const std::shared_ptr<FEDataObjectBase<LoadStepper>> &obj)
    // {
    //     data_objects.push_back(obj);
    // }
};

class LoadStepper : public FEDataObjectBase<EquationSystems>
{
public:
    // Constructor
    LoadStepper(EquationSystems *equation_systems = nullptr)
        : FEDataObjectBase(equation_systems),
          active_load_step(nullptr),
          old_load_step(nullptr),
          load_step_index(0) {}

    // Getters
    std::shared_ptr<LoadStep> get_active_load_step() const { return active_load_step; }
    std::shared_ptr<LoadStep> get_old_load_step() const { return old_load_step; }
    unsigned int get_load_step_index() const { return load_step_index; }

    // // Setters
    // void set_active_load_step(std::shared_ptr<LoadStep> load_step) { active_load_step = load_step; }
    // void set_old_load_step(std::shared_ptr<LoadStep> load_step) { old_load_step = load_step; }
    // void set_load_step_index(unsigned int index) { load_step_index = index; }

    // Switch pointers and read a new LoadStep
    void read_new_load_step(const std::string &filename)
    {
        // Move the current active load step to the old load step
        std::swap(old_load_step, active_load_step);

        // Create a new LoadStep and assign it to the active load step
        active_load_step = std::make_shared<LoadStep>(this);
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
#ifndef SIMULATION_CONTROLLER_H
#define SIMULATION_CONTROLLER_H

#include <memory>
#include "session.h"
#include "element_builder.h"
#include "system_builder.h"
#include "apec_m_linear_system.h"
#include "load_step.h"
// Add other service headers as needed

class SimulationController
{
public:
    SimulationController()
        : session(session)
    {
        // Initialize and collect all main services here
        // Example:
        // element_builder = std::make_shared<ElementBuilder>();
        // system_builder = std::make_shared<SystemBuilder>();
        // apec_m_linear_system = std::make_shared<ApecMLinearSystem>();
        // load_stepper = std::make_shared<LoadStepper>();
    }

    void run_batch_mode(const std::string &mesh_filename, libMesh::Parallel::Communicator &comm)
    {
        // Place your batch mode logic here, using the collected services
        // For example:
        // 1. Build mesh
        // 2. Build equation systems
        // 3. Loop over load steps and solve
    }

    // Add other simulation control methods as needed

    // Service getters (optional)
    // std::shared_ptr<ElementBuilder> get_element_builder() { return element_builder; }
    // ...

private:
    std::shared_ptr<Session> session;
    // std::shared_ptr<ElementBuilder> element_builder;
    // std::shared_ptr<SystemBuilder> system_builder;
    std::shared_ptr<TransientSystem<LinearImplicitSystem>> apec_m_linear_system;
    std::shared_ptr<LoadStepper> load_stepper;
    // Add other services as members
};

#endif // SIMULATION_CONTROLLER_H
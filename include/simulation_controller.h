#ifndef SIMULATION_CONTROLLER_H
#define SIMULATION_CONTROLLER_H

#include <memory>
#include <string>
#include "apec_m_system_objects.h"
#include "apec_m_eq_system_object.h"
#include "load_step.h"
// #include "apec_m_linear_system.h"
// Add other service headers as needed

using namespace libMesh;

// The simulation controller object does the following:
// Reads input model/mesh
// Distributes degrees of freedom based on model information and setups the FEM systems
// Solves the system
// given user input model and load steps, it manages in a logical order all steps of a EMAG simulation.
class SimulationController
{
public:
    explicit SimulationController(libMesh::Parallel::Communicator &comm, bool batch_mode = true, const std::string &mesh_filename = "");

    // CLI command methods
    void read_mesh();
    void clear();
    void physics();
    void solve();
    void finish();
    void exit();

    void set_transient_effects();
    void check_constrained_dofs();
    void set_load_step_data();
    void set_solver_settings();
    void set_boundary_info(const std::set<std::string> &bc_arrays);
    void generate_infinite_elements();

    // Available CLI commands
    enum class Command
    {
        READ_MESH,
        CLEAR,
        PHYSICS,
        SOLVE,
        FINISH,
        EXIT,
        INVALID
    };

    // Execute in batch mode
    void run_batch_mode();

    // Main CLI loop
    void run_cli(libMesh::Parallel::Communicator &comm);

    // Other simulation control methods as needed
    // Helper to parse command from input string
    inline Command parse_command(const std::string &input) const
    {
        if (input.find("-read_mesh") == 0)
            return Command::READ_MESH;
        if (input == "-clear")
            return Command::CLEAR;
        if (input == "-physics")
            return Command::PHYSICS;
        if (input.find("-solve") == 0)
            return Command::SOLVE;
        if (input == "-finish")
            return Command::FINISH;
        if (input == "-exit")
            return Command::EXIT;
        return Command::INVALID;
    }

    // Helper to extract argument after command (e.g., filepath)
    inline std::string extract_argument(const std::string &input) const
    {
        auto pos = input.find(' ');
        if (pos != std::string::npos)
            return input.substr(pos + 1);
        return "";
    }

    // The main attributes of the Simulation contrller, which represent all user input:
    // Mesh/Element Types/etc. in apec_m_model_data
    // The equation systems enerated from the input model
    // The load steps of the simulation to be executed
    libMesh::Parallel::Communicator &comm;
    bool batch_mode;

    std::shared_ptr<FEMDataRepository> apec_m_model_data;

    std::shared_ptr<ApecMEquationSystems> equation_systems;

    std::shared_ptr<LoadStepper> load_stepper;
};

#endif // SIMULATION_CONTROLLER_H
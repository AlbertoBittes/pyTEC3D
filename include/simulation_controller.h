#ifndef SIMULATION_CONTROLLER_H
#define SIMULATION_CONTROLLER_H

#include <memory>
#include <string>
#include "apec_m_linear_system.h"
#include "load_step.h"
// #include "apec_m_linear_system.h"
// Add other service headers as needed

using namespace libMesh;

class SimulationController
{
public:
    explicit SimulationController(libMesh::Parallel::Communicator &comm, bool batch_mode = true);

    // CLI command methods
    void read_mesh();
    void clear();
    void physics();
    void solve();
    void finish();
    void exit();

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

    void run_batch_mode();

    // Main CLI loop
    void run_cli(libMesh::Parallel::Communicator &comm);
    // Add other simulation control methods as needed
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

    // Service getters (optional)
    // std::shared_ptr<ElementBuilder> get_element_builder() { return element_builder; }
    // ...
    libMesh::Parallel::Communicator &comm;

    std::shared_ptr<Session> session;
    std::shared_ptr<LoadStepper> load_stepper;

    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<FEMDataRepository> mesh_input_data;

    std::shared_ptr<ApecMEquationSystems> equation_systems;
    bool batch_mode;
    // Add other services as members
};

#endif // SIMULATION_CONTROLLER_H
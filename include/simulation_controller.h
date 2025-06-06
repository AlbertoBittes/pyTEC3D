#ifndef SIMULATION_CONTROLLER_H
#define SIMULATION_CONTROLLER_H

#include <memory>
#include <string>
#include "session.h"
#include "element_builder.h"
#include "system_builder.h"
#include "apec_m_linear_system.h"
#include "load_step.h"
// Add other service headers as needed

using namespace libMesh;

class SimulationController
{
public:
    explicit SimulationController(libMesh::Parallel::Communicator &comm, bool batch_mode = true)
        : comm(comm), batch_mode(batch_mode)
    {
        // Initialize and collect all main services here
        // Example:
        // element_builder = std::make_shared<ElementBuilder>();
        // system_builder = std::make_shared<SystemBuilder>();
        // apec_m_linear_system = std::make_shared<ApecMLinearSystem>();
        // load_stepper = std::make_shared<LoadStepper>();
        session = std::make_shared<Session>();
        mesh = std::make_shared<Mesh>(comm);
        equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
        equation_systems->init_element_builder();
        load_stepper = std::make_shared<LoadStepper>();
    }

    // CLI command methods
    void read_mesh(const std::string &filepath);
    void clear();
    void physics();
    void solve(int load_step);
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

    void run_batch_mode()
    {
        // Place your batch mode logic here, using the collected services
        // For example:
        // 1. Build mesh
        // 2. Build equation systems
        // 3. Loop over load steps and solve
    }

    // Main CLI loop
    void run_cli(libMesh::Parallel::Communicator &comm)
    {
        std::string input;
        while (true)
        {
            std::cout << "Enter command: ";
            std::getline(std::cin, input);

            Command cmd = parse_command(input);
            switch (cmd)
            {
            case Command::READ_MESH:
            {
                std::string filepath = extract_argument(input);
                session->set_input_mesh_filepath(filepath);
                read_mesh();
                break;
            }
            case Command::CLEAR:
                clear();
                break;
            case Command::PHYSICS:
                physics();
                break;
            case Command::SOLVE:
            {
                std::string step = extract_argument(input);
                load_stepper->read_new_load_step(step);
                solve();
                break;
            }
            case Command::FINISH:
                finish();
                break;
            case Command::EXIT:
                exit();
                return;
            default:
                libMesh::out << "Unknown command.\n";
            }
        }
    }
    // Add other simulation control methods as needed
    // Helper to parse command from input string
    Command parse_command(const std::string &input) const
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
    std::string extract_argument(const std::string &input) const
    {
        auto pos = input.find(' ');
        if (pos != std::string::npos)
            return input.substr(pos + 1);
        return "";
    }

    void read_mesh();

    void physics();

    void solve();
    // Service getters (optional)
    // std::shared_ptr<ElementBuilder> get_element_builder() { return element_builder; }
    // ...
    libMesh::Parallel::Communicator &comm;

    std::shared_ptr<Session> session;
    std::shared_ptr<LoadStepper> load_stepper;

    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<ElementTypeFactory> element_types;

    std::shared_ptr<ApecMEquationSystems> equation_systems;
    bool batch_mode;
    // Add other services as members
};

#endif // SIMULATION_CONTROLLER_H
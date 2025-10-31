#include "apec_m_cli.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>
#include <unordered_map>

#include <libmesh/libmesh.h>

namespace
{
    std::string trim(const std::string &input)
    {
        const auto first = input.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            return "";
        const auto last = input.find_last_not_of(" \t\r\n");
        return input.substr(first, last - first + 1);
    }

    void report_session_paths(const Session &session)
    {
        libMesh::out << "Current mode: " << Session::mode_to_string(session.get_mode()) << std::endl;
        libMesh::out << "Working folder: " << session.get_path("working_folderpath") << std::endl;

        if (!session.get_model_data_filepath().empty())
            libMesh::out << "Model data file: " << session.get_model_data_filepath() << std::endl;

        if (session.has_active_mesh())
        {
            libMesh::out << "Active mesh: " << session.get_active_mesh() << std::endl;
            const auto &mesh_folder = session.get_path("mesh_folderpath");
            if (!mesh_folder.empty())
                libMesh::out << "Mesh folder: " << mesh_folder << std::endl;

            const auto &input_folder = session.get_path("input_folderpath");
            if (!input_folder.empty())
                libMesh::out << "Input folder: " << input_folder << std::endl;

            const auto &output_folder = session.get_path("output_folderpath");
            if (!output_folder.empty())
                libMesh::out << "Output folder: " << output_folder << std::endl;

            const auto &load_steps_folder = session.get_path("load_steps_folderpath");
            if (!load_steps_folder.empty())
                libMesh::out << "Load steps folder: " << load_steps_folder << std::endl;

            const auto &final_output_path = session.get_path("final_output_path");
            if (!final_output_path.empty())
                libMesh::out << "Final output path: " << final_output_path << std::endl;

            const auto &emat_folder = session.get_path("emat_folderpath");
            if (!emat_folder.empty())
                libMesh::out << "Emat folder: " << emat_folder << std::endl;

            const auto &temp_folder = session.get_path("temp_folderpath");
            if (!temp_folder.empty())
                libMesh::out << "Temp folder: " << temp_folder << std::endl;
        }
        else
        {
            libMesh::out << "Active mesh: <none>" << std::endl;
        }
    }

    bool require_mode(const Session &session, Session::Mode expected, const std::string &command_name)
    {
        if (session.get_mode() != expected)
        {
            libMesh::out << command_name << " is only available in "
                         << Session::mode_to_string(expected) << " mode (current: "
                         << Session::mode_to_string(session.get_mode()) << ")." << std::endl;
            return false;
        }
        return true;
    }
}

typedef BasicCommandLineInterface::Command Command;

void BasicCommandLineInterface::run_interactive(APEC_M &program)
{
    std::string input;
    while (true)
    {
        std::cout << "apec_m> " << std::flush;
        if (!std::getline(std::cin, input))
            break;

        auto cleaned = trim(input);
        if (cleaned.empty())
            continue;

        if (!execute_command(program, cleaned, true))
            break;
    }
}

void BasicCommandLineInterface::run_batch(APEC_M &program, const std::string &batch_file)
{
    std::ifstream file(batch_file);
    if (!file)
    {
        libMesh::out << "Unable to open batch file: " << batch_file << std::endl;
        return;
    }

    std::string line;
    unsigned int line_number = 0;
    while (std::getline(file, line))
    {
        ++line_number;
        auto cleaned = trim(line);
        if (cleaned.empty() || cleaned[0] == '#')
            continue;

        if (!execute_command(program, cleaned, false))
        {
            libMesh::out << "Batch execution stopped at line " << line_number << std::endl;
            break;
        }
    }
}

Command BasicCommandLineInterface::parse_command(const std::string &input) const
{
    if (input.rfind("-read_mesh", 0) == 0)
        return Command::READ_MESH;
    if (input == "-clear")
        return Command::CLEAR_MESH;
    if (input.rfind("-set_workdir", 0) == 0)
        return Command::SET_WORKDIR;
    if (input.rfind("-active_mesh", 0) == 0)
        return Command::SET_ACTIVE_MESH;
    if (input.rfind("-read_model_data", 0) == 0)
        return Command::READ_MODEL_DATA;
    if (input == "-distribute_boundary_ids")
        return Command::DISTRIBUTE_BOUNDARIES;
    if (input == "-generate_infinite_elements")
        return Command::GENERATE_INFINITE_ELEMENTS;
    if (input == "-generate_nodal_couplings")
        return Command::GENERATE_NODAL_COUPLINGS;
    if (input == "-generate_periodic_bcs")
        return Command::GENERATE_PERIODIC_BCS;
    if (input.rfind("-mode", 0) == 0)
        return Command::SET_MODE;
    if (input == "-help" || input == "help")
        return Command::HELP;
    if (input == "-exit" || input == "exit" || input == "quit")
        return Command::EXIT;
    return Command::INVALID;
}

std::string BasicCommandLineInterface::extract_argument(const std::string &input) const
{
    auto pos = input.find(' ');
    if (pos == std::string::npos)
        return "";
    return trim(input.substr(pos + 1));
}

bool BasicCommandLineInterface::execute_command(APEC_M &program, const std::string &input, bool interactive)
{
    const auto command = parse_command(input);
    switch (command)
    {
    case Command::READ_MESH:
    {
        auto args = extract_argument(input);
        if (args.empty())
        {
            libMesh::out << "Usage: -read_mesh <mesh_file>" << std::endl;
            return true;
        }

        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-read_mesh"))
            return true;

        auto mesh_name = default_mesh_name();
        auto filepath = args;

        if (interactive)
            libMesh::out << "Reading mesh '" << mesh_name << "' from " << filepath << std::endl;

        program.read_mesh(mesh_name, filepath);
        report_session_paths(session);
        return true;
    }
    case Command::CLEAR_MESH:
    {
        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-clear"))
            return true;

        auto mesh_name = default_mesh_name();
        program.clear_mesh(mesh_name);
        libMesh::out << "Cleared mesh context '" << mesh_name << "'." << std::endl;
        return true;
    }
    case Command::SET_WORKDIR:
    {
        auto arg = extract_argument(input);
        if (arg.empty())
        {
            libMesh::out << "Usage: -set_workdir <folder_path>" << std::endl;
            return true;
        }

        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-set_workdir"))
            return true;

        session.set_working_directory(std::filesystem::path(arg));

        libMesh::out << "Working directory set to: " << session.get_working_directory() << std::endl;
        report_session_paths(session);
        return true;
    }
    case Command::SET_ACTIVE_MESH:
    {
        auto mesh_name = extract_argument(input);
        if (mesh_name.empty())
        {
            libMesh::out << "Usage: -active_mesh <mesh_name>" << std::endl;
            return true;
        }

        if (!program.has_mesh(mesh_name))
        {
            libMesh::out << "Mesh '" << mesh_name << "' is not loaded. Use -read_mesh first." << std::endl;
            return true;
        }

        if (!program.set_active_mesh(mesh_name))
        {
            libMesh::out << "Unable to activate mesh '" << mesh_name << "'." << std::endl;
            return true;
        }

        auto &session = program.session();
        libMesh::out << "Active mesh set to '" << session.get_active_mesh() << "'." << std::endl;
        report_session_paths(session);
        return true;
    }
    case Command::READ_MODEL_DATA:
    {
        auto filepath = extract_argument(input);
        if (filepath.empty())
        {
            libMesh::out << "Usage: -read_model_data <model_file>" << std::endl;
            return true;
        }

        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-read_model_data"))
            return true;

        std::filesystem::path candidate = session.get_working_directory() / filepath;
        if (!std::filesystem::exists(candidate))
        {
            libMesh::out << "Model data file must exist in the working directory: " << candidate << std::endl;
            return true;
        }

        if (!program.load_model_data(filepath))
        {
            libMesh::out << "Failed to load model data from '" << filepath << "'." << std::endl;
            return true;
        }

        libMesh::out << "Model data loaded from '" << candidate << "'." << std::endl;
        report_session_paths(session);
        return true;
    }
    case Command::DISTRIBUTE_BOUNDARIES:
    {
        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-distribute_boundary_ids"))
            return true;

        if (!program.distribute_boundary_ids())
            libMesh::out << "Boundary distribution failed." << std::endl;
        return true;
    }
    case Command::GENERATE_INFINITE_ELEMENTS:
    {
        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-generate_infinite_elements"))
            return true;

        program.generate_infinite_elements();
        return true;
    }
    case Command::GENERATE_NODAL_COUPLINGS:
    {
        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-generate_nodal_couplings"))
            return true;

        program.generate_nodal_couplings();
        return true;
    }
    case Command::GENERATE_PERIODIC_BCS:
    {
        auto &session = program.session();
        if (!require_mode(session, Session::Mode::PreProcessing, "-generate_periodic_bcs"))
            return true;

        if (!program.generate_periodic_boundaries())
            libMesh::out << "Periodic boundary generation failed." << std::endl;
        return true;
    }
    case Command::SET_MODE:
    {
        auto mode_arg = extract_argument(input);
        if (mode_arg.empty())
        {
            libMesh::out << "Usage: -mode <preprocessing|analysis|other>" << std::endl;
            return true;
        }

        static const std::unordered_map<std::string, Session::Mode> mode_map = {
            {"pre", Session::Mode::PreProcessing},
            {"preprocessing", Session::Mode::PreProcessing},
            {"analysis", Session::Mode::Analysis},
            {"other", Session::Mode::Other}};

        std::string key = mode_arg;
        std::transform(key.begin(), key.end(), key.begin(), [](unsigned char c)
                       { return static_cast<char>(std::tolower(c)); });

        auto it = mode_map.find(key);
        if (it == mode_map.end())
        {
            libMesh::out << "Unknown mode '" << mode_arg << "'." << std::endl;
            return true;
        }

        if (!program.set_mode(it->second))
            return true;

        report_session_paths(program.session());
        return true;
    }
    case Command::HELP:
        libMesh::out << "Available commands:\n"
                     << "  -read_mesh <mesh_file>\n"
                     << "  -set_workdir <folder_path>\n"
                     << "  -active_mesh <mesh_name>\n"
                     << "  -distribute_boundary_ids\n"
                     << "  -generate_infinite_elements\n"
                     << "  -generate_nodal_couplings\n"
                     << "  -generate_periodic_bcs\n"
                     << "  -mode <preprocessing|analysis|other>\n"
                     << "  -read_model_data <model_file>\n"
                     << "  -clear\n"
                     << "  -exit" << std::endl;
        return true;
    case Command::EXIT:
        libMesh::out << "Exiting CLI." << std::endl;
        return false;
    case Command::INVALID:
    default:
        libMesh::out << "Unknown command. Type -help for a list of commands." << std::endl;
        return true;
    }
}

std::string BasicCommandLineInterface::default_mesh_name()
{
    return "default";
}

#include "apec_m_cli.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>

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
    libMesh::out << "Working folder: " << session.get_path("working_folderpath") << std::endl;

    if (!session.get_model_data_filepath().empty())
        libMesh::out << "Model data file: " << session.get_model_data_filepath() << std::endl;

    if (session.has_active_mesh())
    {
        libMesh::out << "Active mesh: " << session.get_active_mesh() << std::endl;
        const auto &mesh_folder = session.get_path("mesh_folderpath");
        if (!mesh_folder.empty())
            libMesh::out << "Mesh folder: " << mesh_folder << std::endl;
    }
    else
    {
        libMesh::out << "Active mesh: <none>" << std::endl;
    }

    if (!session.has_active_analysis())
    {
        libMesh::out << "Active analysis: <none>" << std::endl;
        libMesh::out << "Analysis-specific folders will be created once an analysis is instantiated." << std::endl;
        return;
    }

    libMesh::out << "Active analysis: " << session.get_active_analysis() << std::endl;

    const auto &analysis_folder = session.get_path("analysis_folderpath");
    if (!analysis_folder.empty())
        libMesh::out << "Analysis folder: " << analysis_folder << std::endl;

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
    if (input.rfind("-active_analysis", 0) == 0)
        return Command::SET_ACTIVE_ANALYSIS;
    if (input.rfind("-read_model_data", 0) == 0)
        return Command::READ_MODEL_DATA;
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

        auto mesh_name = default_mesh_name();
        auto filepath = args;

        if (interactive)
            libMesh::out << "Reading mesh '" << mesh_name << "' from " << filepath << std::endl;

        program.read_mesh(mesh_name, filepath);
        report_session_paths(program.session());
        return true;
    }
    case Command::CLEAR_MESH:
    {
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
        session.set_working_directory(std::filesystem::path(arg));

        libMesh::out << "Working directory set to: " << session.get_working_directory() << std::endl;
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
        if (session.has_active_analysis())
        {
            libMesh::out << "Active mesh set to '" << session.get_active_mesh() << "' with analysis '"
                         << session.get_active_analysis() << "'." << std::endl;
        }
        else
        {
            libMesh::out << "Active mesh set to '" << session.get_active_mesh() << "'. No active analysis yet." << std::endl;
        }
        report_session_paths(session);
        return true;
    }
    case Command::SET_ACTIVE_ANALYSIS:
    {
        auto analysis_name = extract_argument(input);
        if (analysis_name.empty())
        {
            libMesh::out << "Usage: -active_analysis <analysis_name>" << std::endl;
            return true;
        }

        auto &session = program.session();
        if (!session.has_active_mesh())
        {
            libMesh::out << "No active mesh. Use -active_mesh to select one first." << std::endl;
            return true;
        }

        const auto mesh_name = session.get_active_mesh();
        auto activation = program.activate_analysis(mesh_name, analysis_name);
        if (!activation.has_value())
        {
            libMesh::out << "Mesh '" << mesh_name << "' has no equation system. Read the mesh first." << std::endl;
            return true;
        }

        using ActivationState = ApecMEquationSystems::AnalysisActivationState;
        switch (activation.value())
        {
        case ActivationState::Created:
            libMesh::out << "Analysis '" << analysis_name << "' created and activated for mesh '" << mesh_name << "'." << std::endl;
            break;
        case ActivationState::ActivatedExisting:
            libMesh::out << "Switched to analysis '" << analysis_name << "' for mesh '" << mesh_name << "'." << std::endl;
            break;
        case ActivationState::AlreadyActive:
            libMesh::out << "Analysis '" << analysis_name << "' is already active for mesh '" << mesh_name << "'." << std::endl;
            break;
        }

        report_session_paths(session);
        return true;
    }
    case Command::HELP:
        libMesh::out << "Available commands:\n"
                     << "  -read_mesh <mesh_file>\n"
                     << "  -set_workdir <folder_path>\n"
                     << "  -active_mesh <mesh_name>\n"
                     << "  -active_analysis <analysis_name>\n"
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
*** End of File

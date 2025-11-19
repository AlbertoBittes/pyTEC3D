#include "apec_m_cli.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>
#include <system_error>
#include <unordered_map>

#include <libmesh/libmesh.h>
#include <libmesh/parallel.h>

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
        const auto normalize_path = [](const std::filesystem::path &path) -> std::filesystem::path
        {
            if (path.empty())
                return path;

            std::error_code ec;
            auto absolute_path = std::filesystem::absolute(path, ec);
            if (!ec)
                return absolute_path.lexically_normal();

            return path.lexically_normal();
        };

        const auto print_path = [&](const std::string &label, const std::filesystem::path &path)
        {
            if (path.empty())
                return;
            libMesh::out << label << normalize_path(path) << std::endl;
        };

        libMesh::out << "Current mode: " << Session::mode_to_string(session.get_mode()) << std::endl;
        print_path("Working folder: ", session.get_path("working_folderpath"));

        if (!session.get_model_data_filepath().empty())
            print_path("Model data file: ", session.get_model_data_filepath());

        if (session.has_active_mesh())
        {
            libMesh::out << "Active mesh: " << session.get_active_mesh() << std::endl;
            print_path("Mesh folder: ", session.get_path("mesh_folderpath"));
            print_path("Input folder: ", session.get_path("input_folderpath"));
            print_path("Output folder: ", session.get_path("output_folderpath"));
            print_path("Load steps folder: ", session.get_path("load_steps_folderpath"));
            print_path("Final output path: ", session.get_path("final_output_path"));
            print_path("Emat folder: ", session.get_path("emat_folderpath"));
            print_path("Temp folder: ", session.get_path("temp_folderpath"));
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
    auto &comm = program.communicator();
    const unsigned int root = 0;

    bool keep_running = true;
    while (keep_running)
    {
        std::string input;

        if (comm.rank() == static_cast<int>(root))
        {
            std::cout << "apec_m> " << std::flush;
            if (!std::getline(std::cin, input))
                input = "-exit";
            input = trim(input);
        }

        comm.broadcast(input, root);

        if (input.empty())
            continue;

        bool continue_loop = execute_command(program, input, comm.rank() == static_cast<int>(root));
        comm.broadcast(continue_loop, root);
        keep_running = continue_loop;
    }
}

void BasicCommandLineInterface::run_batch(APEC_M &program, const std::string &batch_file)
{
    auto &comm = program.communicator();
    const unsigned int root = 0;

    std::ifstream file;
    bool file_ok = true;
    if (comm.rank() == static_cast<int>(root))
    {
        file.open(batch_file);
        if (!file)
        {
            libMesh::out << "Unable to open batch file: " << batch_file << std::endl;
            file_ok = false;
        }
    }

    comm.broadcast(file_ok, root);
    if (!file_ok)
        return;

    unsigned int line_number = 0;
    bool keep_running = true;
    while (keep_running)
    {
        std::string command;

        if (comm.rank() == static_cast<int>(root))
        {
            std::string line;
            bool found_command = false;
            while (std::getline(file, line))
            {
                ++line_number;
                auto cleaned = trim(line);
                if (cleaned.empty() || cleaned[0] == '#')
                    continue;
                command = cleaned;
                found_command = true;
                break;
            }

            if (!found_command)
                command = "-exit";
        }

        comm.broadcast(command, root);

        if (command.empty())
            continue;

        bool continue_loop = execute_command(program, command, comm.rank() == static_cast<int>(root));
        comm.broadcast(continue_loop, root);
        keep_running = continue_loop;

        if (!keep_running && comm.rank() == static_cast<int>(root) && command != "-exit")
            libMesh::out << "Batch execution stopped at line " << line_number << std::endl;
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
    if (input == "-print_model_data")
        return Command::PRINT_MODEL_DATA;
    if (input == "-list_load_steps")
        return Command::LIST_LOAD_STEPS;
    if (input.rfind("-solve_transient", 0) == 0)
        return Command::SOLVE_TRANSIENT;
    if (input.rfind("-solve", 0) == 0)
        return Command::SOLVE;
    if (input == "-finish")
        return Command::FINISH;
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

        const Session &session = *program.model_data()->get_object<Session>();
        if (!require_mode(session, Session::Mode::PreProcessing, "-read_mesh"))
            return true;

        auto filepath = args;
        std::filesystem::path mesh_path(filepath);
        if (mesh_path.is_relative())
            mesh_path = std::filesystem::path(session.get_working_directory()) / mesh_path;
        auto mesh_name = mesh_path.stem().string();
        if (mesh_name.empty())
            mesh_name = default_mesh_name();

        if (interactive)
            libMesh::out << "Reading mesh '" << mesh_name << "' from " << filepath << std::endl;

        program.read_mesh(mesh_name, filepath);
        report_session_paths(session);
        return true;
    }
    case Command::CLEAR_MESH:
    {
        const Session &session = *program.model_data()->get_object<Session>();
        if (!require_mode(session, Session::Mode::PreProcessing, "-clear"))
            return true;

        std::string mesh_name = session.has_active_mesh() ? session.get_active_mesh() : default_mesh_name();
        if (mesh_name.empty())
            mesh_name = default_mesh_name();

        if (!program.has_mesh(mesh_name))
        {
            libMesh::out << "Mesh '" << mesh_name << "' is not loaded. Use -read_mesh first." << std::endl;
            return true;
        }

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

        Session &session = *program.model_data()->get_object<Session>();
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

        const Session &session = *program.model_data()->get_object<Session>();
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

        const Session &session = *program.model_data()->get_object<Session>();
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
        program.print_model_data();
        return true;
    }
    case Command::PRINT_MODEL_DATA:
    {
        const Session &session = *program.model_data()->get_object<Session>();
        // if (!require_mode(session, Session::Mode::PreProcessing, "-print_model_data"))
        //     return true;

        program.print_model_data();
        return true;
    }
    case Command::LIST_LOAD_STEPS:
    {
        const Session &session = *program.model_data()->get_object<Session>();
        if (!require_mode(session, Session::Mode::Analysis, "-list_load_steps"))
            return true;

        const auto files = program.list_load_step_files();
        if (files.empty())
        {
            libMesh::out << "No load step binaries found for the active mesh." << std::endl;
            return true;
        }

        libMesh::out << files.size() << " load step file(s) available:" << std::endl;
        for (const auto &path : files)
            libMesh::out << "  - " << path.filename().string() << std::endl;
        return true;
    }
    case Command::SOLVE:
    {
        const Session &session = *program.model_data()->get_object<Session>();
        if (!require_mode(session, Session::Mode::Analysis, "-solve"))
            return true;

        auto argument = extract_argument(input);
        if (!program.solve_load_step(argument))
            libMesh::out << "Load step solve failed." << std::endl;
        return true;
    }
    case Command::SOLVE_TRANSIENT:
    {
        const Session &session = *program.model_data()->get_object<Session>();
        if (!require_mode(session, Session::Mode::Analysis, "-solve_transient"))
            return true;

        auto argument = extract_argument(input);
        if (!program.solve_transient_load_steps(argument))
            libMesh::out << "Transient solve failed." << std::endl;
        return true;
    }
    case Command::FINISH:
    {
        const Session &session = *program.model_data()->get_object<Session>();
        if (!require_mode(session, Session::Mode::Analysis, "-finish"))
            return true;

        if (!program.finish_analysis())
            libMesh::out << "Failed to reset analysis state." << std::endl;
        else
            libMesh::out << "Analysis state cleared." << std::endl;
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

        report_session_paths(*program.model_data()->get_object<Session>());
        return true;
    }
    case Command::HELP:
        libMesh::out << "Available commands:\n"
                     << "  -read_mesh <mesh_file>\n"
                     << "  -set_workdir <folder_path>\n"
                     << "  -active_mesh <mesh_name>\n"
                     << "  -list_load_steps\n"
                     << "  -solve [load_step_filename]\n"
                     << "  -solve_transient [base_name]\n"
                     << "  -finish\n"
                     << "  -mode <preprocessing|analysis|other>\n"
                     << "  -read_model_data <model_file>\n"
                     << "  -print_model_data\n"
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

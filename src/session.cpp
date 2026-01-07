#include "session.h"

#include <cctype>
#include <iostream>
#include <string>

#include <libmesh/libmesh.h>

Session::Session()
    : input_mesh_filepath(""),
      mesh_filename_without_extension(""),
      output_filename(""),
      active_mesh(""),
      working_directory(std::filesystem::current_path()),
      program_mode(Mode::PreProcessing)
{
    initialize_folder_paths();
}

void Session::initialize_folder_paths()
{
    update_dynamic_paths();

    auto &paths_map = *this->get_data_map();

    const auto &mesh_folder = paths_map["mesh_folderpath"];
    if (mesh_folder.empty())
        return;

    std::filesystem::create_directories(mesh_folder);

    const auto &input_folder = paths_map["input_folderpath"];
    if (!input_folder.empty())
        std::filesystem::create_directories(input_folder);

    const auto &load_steps_folder = paths_map["load_steps_folderpath"];
    if (!load_steps_folder.empty())
        std::filesystem::create_directories(load_steps_folder);

    const auto &time_steps_folder = paths_map["time_steps_folderpath"];
    if (!time_steps_folder.empty())
        std::filesystem::create_directories(time_steps_folder);

    const auto &model_data_folder = paths_map["model_data_folderpath"];
    if (!model_data_folder.empty())
        std::filesystem::create_directories(model_data_folder);

    create_temp_folder();
    create_output_folder();
    create_emat_folder();
    check_load_steps_folder();
}

void Session::set_input_mesh_filepath(const std::string &filepath)
{
    input_mesh_filepath = filepath;
    if (filepath.empty())
    {
        update_dynamic_paths();
        return;
    }

    std::filesystem::path input_path(filepath);
    mesh_filename_without_extension = input_path.stem().string();

    update_dynamic_paths();
}

void Session::set_active_mesh(const std::string &mesh_name)
{
    active_mesh = mesh_name;
    if (!mesh_name.empty())
        mesh_filename_without_extension = mesh_name;
    else
        mesh_filename_without_extension.clear();

    output_filename.clear();
    update_dynamic_paths();
}

void Session::set_mode(Mode mode)
{
    program_mode = mode;
}

std::string Session::mode_to_string(Mode mode)
{
    switch (mode)
    {
    case Mode::PreProcessing:
        return "PreProcessing";
    case Mode::Analysis:
        return "Analysis";
    case Mode::Other:
        return "Other";
    }
    return "Unknown";
}

void Session::set_path(const std::string &key, const std::filesystem::path &path)
{
    (*data_map)[key] = path;
}

const std::filesystem::path &Session::get_path(const std::string &key) const
{
    return data_map->at(key);
}

std::filesystem::path Session::resolve_step_path(const std::string &folder_key,
                                                 const std::string &filename) const
{
    const auto &folder = get_path(folder_key);
    if (folder.empty())
        return {};

    if (filename.empty())
        return {};

    auto candidate = std::filesystem::path(filename);
    if (!candidate.is_absolute())
        candidate = folder / candidate;

    if (!candidate.has_extension())
        candidate.replace_extension(".bin");

    return candidate;
}

std::filesystem::path Session::resolve_default_step_base(const std::string &folder_key,
                                                         const std::string &default_token,
                                                         const std::string &suffix) const
{
    const auto &folder = get_path(folder_key);
    if (folder.empty())
        return {};

    auto token = get_mesh_filename_without_extension();
    if (token.empty())
        token = default_token;

    std::filesystem::path base = folder / (token + suffix);
    return base;
}

void Session::set_output_filename(const std::string &filename)
{
    output_filename = filename;
    update_dynamic_paths();
}

void Session::set_output_filepath(unsigned int output_label, std::optional<unsigned int> picard_iteration)
{
    if (!has_active_mesh())
    {
        libMesh::out << "No active mesh. Skipping output filepath assignment." << std::endl;
        output_filename.clear();
        update_dynamic_paths();
        return;
    }

    std::string base_name = mesh_filename_without_extension;
    if (base_name.empty())
    {
        base_name = active_mesh.empty() ? std::string("mesh") : sanitize_token(active_mesh);
    }

    output_filename = base_name + "_out_" + std::to_string(output_label);
    if (picard_iteration)
        output_filename += "_" + std::to_string(*picard_iteration);
    output_filename += ".pvtu";
    update_dynamic_paths();
}

void Session::create_temp_folder()
{
    auto &paths_map = *this->get_data_map();
    const auto &temp_folderpath = paths_map["temp_folderpath"];
    if (temp_folderpath.empty())
        return;

    if (!std::filesystem::exists(temp_folderpath))
    {
        std::filesystem::create_directories(temp_folderpath);
        libMesh::out << "Temporary folder created at " << temp_folderpath << std::endl;
    }
    else
    {
        libMesh::out << "Temporary folder already exists at " << temp_folderpath << std::endl;
    }
}

void Session::create_output_folder()
{
    auto &paths_map = *this->get_data_map();
    const auto &output_folderpath = paths_map["output_folderpath"];
    if (output_folderpath.empty())
        return;

    if (!std::filesystem::exists(output_folderpath))
    {
        std::filesystem::create_directories(output_folderpath);
        libMesh::out << "Output folder created at " << output_folderpath << std::endl;
    }
    else
    {
        libMesh::out << "Output folder already exists at " << output_folderpath << std::endl;
    }
}

void Session::create_emat_folder()
{
    auto &paths_map = *this->get_data_map();
    const auto &emat_folderpath = paths_map["emat_folderpath"];
    if (emat_folderpath.empty())
        return;

    if (!std::filesystem::exists(emat_folderpath))
    {
        std::filesystem::create_directories(emat_folderpath);
        libMesh::out << "Emat folder created at " << emat_folderpath << std::endl;
        emat_folders_cleared.insert(emat_folderpath);
    }
    else
    {
        if (emat_folders_cleared.insert(emat_folderpath).second)
        {
            for (const auto &entry : std::filesystem::directory_iterator(emat_folderpath))
            {
                if (entry.is_regular_file())
                    std::filesystem::remove(entry.path());
            }
            libMesh::out << "Emat folder already exists at " << emat_folderpath << ". All files deleted." << std::endl;
        }
        else
        {
            libMesh::out << "Emat folder already exists at " << emat_folderpath << "." << std::endl;
        }
    }
}

void Session::check_load_steps_folder()
{
    auto &paths_map = *this->get_data_map();
    const auto &load_steps_folderpath = paths_map["load_steps_folderpath"];
    if (load_steps_folderpath.empty())
        return;
    if (!std::filesystem::exists(load_steps_folderpath))
    {
        std::filesystem::create_directories(load_steps_folderpath);
        libMesh::out << "Input folder created at " << load_steps_folderpath << std::endl;
    }
}

void Session::print_data_info() const
{
    libMesh::out << "  Mode: " << mode_to_string(program_mode) << std::endl;
    libMesh::out << "  Working directory: " << (working_directory.empty() ? "<unset>" : working_directory.string()) << std::endl;
    libMesh::out << "  Active mesh: " << (active_mesh.empty() ? "<none>" : active_mesh) << std::endl;
    libMesh::out << "  Output filename: " << (output_filename.empty() ? "<unset>" : output_filename) << std::endl;

    const auto *paths = this->get_data_map();
    if (paths && !paths->empty())
    {
        libMesh::out << "  Registered paths:" << std::endl;
        for (const auto &entry : *paths)
            libMesh::out << "    " << entry.first << " -> " << entry.second << std::endl;
    }
    else
        libMesh::out << "  Registered paths: none." << std::endl;
}

void Session::update_dynamic_paths()
{
    auto &paths_map = *this->get_data_map();

    const auto effective_working_directory = working_directory.empty() ? std::filesystem::current_path() : working_directory;
    const auto mesh_token = sanitize_token(active_mesh.empty() ? mesh_filename_without_extension : active_mesh);

    std::filesystem::path mesh_root;
    if (!mesh_token.empty())
        mesh_root = effective_working_directory / mesh_token;

    paths_map["working_folderpath"] = effective_working_directory;
    paths_map["time_steps_folderpath"] = effective_working_directory / "_time_";
    paths_map["model_data_folderpath"] = effective_working_directory / "_model_data_";
    paths_map["mesh_folderpath"] = mesh_root;

    if (!mesh_root.empty())
    {
        paths_map["input_folderpath"] = mesh_root / "input";
        paths_map["output_folderpath"] = mesh_root / "output";
        paths_map["emat_folderpath"] = mesh_root / "emat";
        paths_map["temp_folderpath"] = mesh_root / "temp";
        paths_map["load_steps_folderpath"] = paths_map["input_folderpath"];

        if (!output_filename.empty())
            paths_map["final_output_path"] = paths_map["output_folderpath"] / output_filename;
        else
            paths_map["final_output_path"] = std::filesystem::path();
    }
    else
    {
        paths_map["input_folderpath"] = std::filesystem::path();
        paths_map["output_folderpath"] = std::filesystem::path();
        paths_map["emat_folderpath"] = std::filesystem::path();
        paths_map["temp_folderpath"] = std::filesystem::path();
        paths_map["load_steps_folderpath"] = std::filesystem::path();
        paths_map["final_output_path"] = std::filesystem::path();
    }

    paths_map.erase("analysis_folderpath");
}

std::string Session::sanitize_token(const std::string &value)
{
    std::string result;
    result.reserve(value.size());
    for (char ch : value)
    {
        if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' || ch == '-')
            result.push_back(ch);
        else
            result.push_back('_');
    }
    return result;
}

#include "../include/session.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <regex>
#include <string>
#include "../include/json.hpp"
#include <libmesh/libmesh.h>

Session::Session(const std::string &configFile)
{
    read_config(configFile);
    input_mesh_filepath = config["Mesh Input Section"]["mesh_filepath"].get<std::string>();
    libMesh::out << "Setting output path..." << std::endl;
    mesh_folderpath = std::filesystem::path(input_mesh_filepath).parent_path();

    if (config["Output Section"]["output_folderpath"].is_null())
    {
        // Create an output folder path based on the input mesh file path
        output_folderpath = (mesh_folderpath / "output");
        create_output_folder();
    }
    else
        output_folderpath = std::filesystem::path(config["Output Section"]["output_filepath"].get<std::string>());

    // Use filesystem to modify the basename
    mesh_path = std::filesystem::path(input_mesh_filepath);
    mesh_filename_without_extension = mesh_path.stem().string();

    create_temp_folder();
}

void Session::read_config(const std::string &configFile)
{
    std::ifstream i(configFile);
    if (!i.is_open())
    {
        libmesh_error_msg("Could not open config file: " + configFile);
    }
    config = nlohmann::json::parse(i);

    // Perform checks on the configuration
    // checkTransientAnalysisSettings();
    // checkFilePaths();
    // checkStringValues();
    // checkBooleanValues();
    // checkMandatoryVariables();
}

void Session::set_output_filepath(unsigned int output_label)
{
    output_filename = mesh_filename_without_extension + "_out_" + std::to_string(output_label) + ".pvtu";
    final_output_path = output_folderpath / output_filename;
    // std::filesystem::path output_path(output_folderpath);
    // std::filesystem::path mesh_path(mesh_filepath);
    // std::filesystem::path final_output_path;
    // std::string mesh_filename_without_extension = mesh_path.stem().string();
    // std::string output_filename;
}

void Session::create_temp_folder()
{
    temp_folder_path = std::filesystem::current_path() / "temp_files";
    if (!std::filesystem::exists(temp_folder_path))
    {
        std::filesystem::create_directory(temp_folder_path);
        libMesh::out << "Temporary folder 'temp_files' created at " << temp_folder_path << std::endl;
    }
    else
    {
        libMesh::out << "Temporary folder 'temp_files' already exists at " << temp_folder_path << std::endl;
    }
}

void Session::create_output_folder()
{
    if (!std::filesystem::exists(output_folderpath))
    {
        std::filesystem::create_directory(output_folderpath);
        libMesh::out << "Output folder created at " << output_path << std::endl;
    }
    else
    {
        libMesh::out << "Output folder already exists at " << output_path << std::endl;
    }
}
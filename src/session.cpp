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
    mesh_filepath = config["Mesh Input Section"]["mesh_filepath"].get<std::string>();
    if (config["Output Section"]["output_filepath"].is_null())
        output_folderpath = std::filesystem::path(mesh_filepath).parent_path().string();
    else
        output_folderpath = config["Output Section"]["output_filepath"].get<std::string>();
    libMesh::out << "Setting output path..." << std::endl;
    // Use filesystem to modify the basename
    output_path = std::filesystem::path(output_folderpath);
    mesh_path = std::filesystem::path(mesh_filepath);
    mesh_filename_without_extension = mesh_path.stem().string();
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
    final_output_path = output_path / output_filename;
    // std::filesystem::path output_path(output_folderpath);
    // std::filesystem::path mesh_path(mesh_filepath);
    // std::filesystem::path final_output_path;
    // std::string mesh_filename_without_extension = mesh_path.stem().string();
    // std::string output_filename;
}
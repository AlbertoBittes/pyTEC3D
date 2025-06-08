#include "session.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <string>
#include <filesystem>
#include <libmesh/libmesh.h>

std::unique_ptr<std::ios> Session::open_file_stream(const std::filesystem::path &folder,
                                                    const std::string &filename,
                                                    FileMode mode)
{
    std::filesystem::path filepath = folder / filename;
    std::unique_ptr<std::ios> stream;

    if (mode == FileMode::Read)
    {
        auto ifs = std::make_unique<std::ifstream>(filepath, std::ios::binary);
        if (!ifs->is_open())
            libmesh_error_msg("Failed to open file for reading: " + filepath.string());
        stream = std::move(ifs);
    }
    else // FileMode::Write
    {
        auto ofs = std::make_unique<std::ofstream>(filepath, std::ios::binary);
        if (!ofs->is_open())
            libmesh_error_msg("Failed to open file for writing: " + filepath.string());
        stream = std::move(ofs);
    }
    return stream;
}

std::vector<std::unique_ptr<std::ios>> Session::open_file_streams(
    const std::filesystem::path &folder,
    const std::vector<std::string> &filenames,
    FileMode mode)
{
    std::vector<std::unique_ptr<std::ios>> streams;
    streams.reserve(filenames.size());
    for (const auto &fname : filenames)
    {
        streams.push_back(open_file_stream(folder, fname, mode));
    }
    return streams;
}

void Session::close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams)
{
    for (auto &stream : streams)
    {
        if (stream)
        {
            // Try to cast to ifstream
            if (auto *ifs = dynamic_cast<std::ifstream *>(stream.get()))
            {
                if (ifs->is_open())
                    ifs->close();
            }
            // Try to cast to ofstream
            else if (auto *ofs = dynamic_cast<std::ofstream *>(stream.get()))
            {
                if (ofs->is_open())
                    ofs->close();
            }
        }
    }
}

void Session::set_output_filepath(unsigned int output_label)
{
    auto paths_map = *this->get_data_map();
    auto &output_folderpath = paths_map["output_folderpath"];
    auto &final_output_path = paths_map["final_output_path"];
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
    auto paths_map = *this->get_data_map();
    auto &temp_folderpath = paths_map["temp_folderpath"];
    if (!std::filesystem::exists(temp_folderpath))
    {
        std::filesystem::create_directory(temp_folderpath);
        libMesh::out << "Temporary folder 'temp_files' created at " << temp_folderpath << std::endl;
    }
    else
    {
        libMesh::out << "Temporary folder 'temp_files' already exists at " << temp_folderpath << std::endl;
    }
}

void Session::create_output_folder()
{
    auto paths_map = *this->get_data_map();
    auto &output_folderpath = paths_map["output_folderpath"];
    if (!std::filesystem::exists(output_folderpath))
    {
        std::filesystem::create_directory(output_folderpath);
        libMesh::out << "Output folder created at " << output_folderpath << std::endl;
    }
    else
    {
        libMesh::out << "Output folder already exists at " << output_folderpath << std::endl;
    }
}

void Session::create_emat_folder()
{
    auto paths_map = *this->get_data_map();
    auto &emat_folderpath = paths_map["emat_folderpath"];
    if (!std::filesystem::exists(emat_folderpath))
    {
        std::filesystem::create_directory(emat_folderpath);
        libMesh::out << "Emat folder created at " << emat_folderpath << std::endl;
    }
    else
    {
        // Delete all files in the folder
        for (const auto &entry : std::filesystem::directory_iterator(emat_folderpath))
        {
            if (entry.is_regular_file())
                std::filesystem::remove(entry.path());
        }
        libMesh::out << "Emat folder already exists at " << emat_folderpath << ". All files deleted." << std::endl;
    }
}

void Session::check_load_steps_folder()
{
    auto paths_map = *this->get_data_map();
    auto &load_steps_folderpath = paths_map["load_steps_folderpath"];
    if (!std::filesystem::exists(load_steps_folderpath))
    {
        libmesh_error_msg("Load steps folder does not exist");
    }
    else
    {
        libMesh::out << "Load steps folder already exists at " << load_steps_folderpath << std::endl;

        // Count the number of load step files
        unsigned int load_step_count = 0;
        for (unsigned int i = 1;; ++i)
        {
            std::filesystem::path load_step_file = load_steps_folderpath / ("load_step_" + std::to_string(i) + ".bin");
            if (std::filesystem::exists(load_step_file))
            {
                ++load_step_count;
            }
            else
            {
                // Stop counting if a file is missing in the sequence
                break;
            }
        }

        // Set the number of load steps
        // if (load_step_count > 0)
        // {
        n_load_steps = load_step_count;
        libMesh::out << "Number of load steps: " << n_load_steps << std::endl;
        // }
        // else
        // {
        //     libmesh_error_msg("No valid load step files found in the load steps folder.");
        // }
    }
}
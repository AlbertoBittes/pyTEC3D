#ifndef SESSION_H
#define SESSION_H

#include <filesystem>
#include <string>
#include <libmesh/libmesh.h>
#include "simulation_controller.h"

using namespace libMesh;

enum FileMode
{
    Read,
    Write
};

class Session : public FEDataObjectBase<SimulationController>
{
public:
    // Default constructor
    Session()
        : input_mesh_filepath(""),
          input_directory(""),
          output_folderpath("output"),
          temp_folderpath("temp"),
          load_steps_folderpath("load_steps"),
          emat_folderpath("emat_files")
    {
    }

    // Method to initialize folder paths based on input_mesh_filepath
    void initialize_folder_paths()
    {
        // Convert input_mesh_filepath to a filesystem path and extract the directory
        std::filesystem::path input_path = std::filesystem::path(input_mesh_filepath);
        input_directory = input_path.parent_path(); // Store the directory as an attribute

        // Prepend the extracted directory to the folder paths
        output_folderpath = input_directory / "output";
        temp_folderpath = input_directory / "temp";
        load_steps_folderpath = input_directory / "load_steps";
        emat_folderpath = input_directory / "emat_files";

        this->create_temp_folder();
        this->create_output_folder();
        this->create_emat_folder();

        this->check_load_steps_folder();
    }

    // Setters
    void set_input_mesh_filepath(const std::string &filepath)
    {
        input_mesh_filepath = filepath;
        initialize_folder_paths(); // Reinitialize folder paths when input filepath changes
    }
    void set_output_folderpath(const std::filesystem::path &path) { output_folderpath = path; }
    void set_output_path(const std::filesystem::path &path) { output_path = path; }
    void set_mesh_filename_without_extension(const std::string &filename) { mesh_filename_without_extension = filename; }
    void set_output_filename(const std::string &filename) { output_filename = filename; }
    void set_final_output_path(const std::filesystem::path &path) { final_output_path = path; }
    void set_temp_folder_path(const std::filesystem::path &path) { temp_folderpath = path; }
    void set_load_steps_folderpath(const std::filesystem::path &path) { load_steps_folderpath = path; }
    void set_emat_folderpath(const std::filesystem::path &path) { emat_folderpath = path; }
    void set_n_load_steps(unsigned int n) { n_load_steps = n; }

    // Getters
    const std::string &get_input_mesh_filepath() const { return input_mesh_filepath; }
    const std::filesystem::path &get_input_directory() const { return input_directory; }
    const std::filesystem::path &get_output_folderpath() const { return output_folderpath; }
    const std::filesystem::path &get_output_path() const { return output_path; }
    const std::string &get_mesh_filename_without_extension() const { return mesh_filename_without_extension; }
    const std::string &get_output_filename() const { return output_filename; }
    const std::filesystem::path &get_final_output_path() const { return final_output_path; }
    const std::filesystem::path &get_temp_folder_path() const { return temp_folderpath; }
    const std::filesystem::path &get_load_steps_folderpath() const { return load_steps_folderpath; }
    const std::filesystem::path &get_emat_folderpath() const { return emat_folderpath; }
    unsigned int get_n_load_steps() const { return n_load_steps; }

    // Methods
    void set_output_filepath(unsigned int output_label);
    void create_temp_folder();
    void create_output_folder();
    void create_emat_folder();
    void check_load_steps_folder();

    std::unique_ptr<std::ios> open_file_stream(const std::filesystem::path &folder,
                                               const std::string &filename,
                                               FileMode mode);

    std::vector<std::unique_ptr<std::ios>> open_file_streams(
        const std::filesystem::path &folder,
        const std::vector<std::string> &filenames,
        FileMode mode);

    void close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams);

private:
    std::string input_mesh_filepath;
    std::filesystem::path input_directory; // New attribute to store the extracted directory
    std::filesystem::path output_folderpath;
    std::filesystem::path output_path;
    std::string mesh_filename_without_extension;
    std::string output_filename;
    std::filesystem::path final_output_path;
    std::filesystem::path temp_folderpath;
    std::filesystem::path load_steps_folderpath;
    std::filesystem::path emat_folderpath;
    unsigned int n_load_steps;
};

#endif // SESSION_H
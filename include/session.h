#ifndef SESSION_H
#define SESSION_H

#include <filesystem>
#include <string>
#include <libmesh/libmesh.h>

using namespace libMesh;

class Session
{
public:
    // Default constructor
    Session()
        : batch_mode(false),
          input_mesh_filepath(""),
          output_folderpath("output"),
          temp_folderpath("temp"),
          load_steps_folderpath("load_steps"),
          emat_folderpath("emat_files") {}

    // Constructor with mesh filename
    explicit Session(const std::string &mesh_filename)
        : batch_mode(false),
          input_mesh_filepath(mesh_filename),
          mesh_filename_without_extension(std::filesystem::path(mesh_filename).stem().string()),
          output_folderpath("output"),
          temp_folderpath("temp"),
          load_steps_folderpath("load_steps"),
          emat_folderpath("emat_files")
    {
    }

    // Constructor with mesh filename and batch mode
    Session(const std::string &mesh_filename, bool batch_mode)
        : batch_mode(batch_mode),
          input_mesh_filepath(mesh_filename),
          mesh_filename_without_extension(std::filesystem::path(mesh_filename).stem().string()),
          output_folderpath("output"),
          temp_folderpath("temp"),
          load_steps_folderpath("load_steps"),
          emat_folderpath("emat_files")
    {
    }

    // Setters
    void set_batch_mode(bool mode) { batch_mode = mode; }
    void set_input_mesh_filepath(const std::string &filepath) { input_mesh_filepath = filepath; }
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
    bool get_batch_mode() const { return batch_mode; }
    const std::string &get_input_mesh_filepath() const { return input_mesh_filepath; }
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

private:
    bool batch_mode;
    std::string input_mesh_filepath;
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
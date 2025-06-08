#ifndef SESSION_H
#define SESSION_H

#include <filesystem>
#include <string>
#include <unordered_map>

#include "fe_data_object.h" // or wherever FEDataObjectMap is defined
#include "globals.h"        // or wherever PathMap is defined

typedef std::unordered_map<std::string, std::filesystem::path> PathMap;

enum FileMode
{
    Read,
    Write
};

class Session : public FEDataObjectMap<PathMap>
{
public:
    Session()
        : input_mesh_filepath(""),
          mesh_filename_without_extension(""),
          output_filename(""),
          n_load_steps(0)
    {
        auto paths_map = *this->get_data_map();
        // Initialize default paths in the map
        paths_map["output_folderpath"] = "output";
        paths_map["temp_folderpath"] = "temp";
        paths_map["load_steps_folderpath"] = "load_steps";
        paths_map["emat_folderpath"] = "emat_files";
    }

    void initialize_folder_paths()
    {
        auto paths_map = *this->get_data_map();

        std::filesystem::path input_path = std::filesystem::path(input_mesh_filepath);
        paths_map["input_directory"] = input_path.parent_path();

        // Prepend the extracted directory to the folder paths
        paths_map["output_folderpath"] = paths_map["input_directory"] / "output";
        paths_map["temp_folderpath"] = paths_map["input_directory"] / "temp";
        paths_map["load_steps_folderpath"] = paths_map["input_directory"] / "load_steps";
        paths_map["emat_folderpath"] = paths_map["input_directory"] / "emat_files";

        create_temp_folder();
        create_output_folder();
        create_emat_folder();
        check_load_steps_folder();
    }

    // Setters
    void set_input_mesh_filepath(const std::string &filepath)
    {
        input_mesh_filepath = filepath;
        initialize_folder_paths();
    }
    void set_path(const std::string &key, const std::filesystem::path &path) { paths[key] = path; }
    void set_mesh_filename_without_extension(const std::string &filename) { mesh_filename_without_extension = filename; }
    void set_output_filename(const std::string &filename) { output_filename = filename; }
    void set_n_load_steps(unsigned int n) { n_load_steps = n; }

    // Getters
    const std::string &get_input_mesh_filepath() const { return input_mesh_filepath; }
    const std::filesystem::path &get_path(const std::string &key) const { return paths.at(key); }
    const std::string &get_mesh_filename_without_extension() const { return mesh_filename_without_extension; }
    const std::string &get_output_filename() const { return output_filename; }
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
    std::string mesh_filename_without_extension;
    std::string output_filename;
    unsigned int n_load_steps;
    std::unordered_map<std::string, std::filesystem::path> paths;
};

#endif // SESSION_H
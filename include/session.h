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

// The session class takes care of interacting with data files, by properly storing the filepaths
// of important files, such as the element matrices file, the output folder, the temporary files folder, etc.
class Session : public FEDataObjectMap<PathMap>
{
public:
    enum class Mode
    {
        PreProcessing,
        Analysis,
        Other
    };

    Session();

    void initialize_folder_paths();

    void set_input_mesh_filepath(const std::string &filepath);
    const std::string &get_input_mesh_filepath() const { return input_mesh_filepath; }

    void set_active_mesh(const std::string &mesh_name);
    const std::string &get_active_mesh() const { return active_mesh; }
    bool has_active_mesh() const { return !active_mesh.empty(); }

    void set_working_directory(const std::filesystem::path &directory);
    const std::filesystem::path &get_working_directory() const { return working_directory; }

    void set_model_data_filepath(const std::filesystem::path &filepath);
    const std::filesystem::path &get_model_data_filepath() const { return model_data_filepath; }

    void set_mode(Mode mode);
    Mode get_mode() const { return program_mode; }
    static std::string mode_to_string(Mode mode);

    void set_path(const std::string &key, const std::filesystem::path &path);
    const std::filesystem::path &get_path(const std::string &key) const;

    const std::string &get_mesh_filename_without_extension() const { return mesh_filename_without_extension; }

    void set_output_filename(const std::string &filename);
    const std::string &get_output_filename() const { return output_filename; }

    void set_output_filepath(unsigned int output_label);
    void create_temp_folder();
    void create_output_folder();
    void create_emat_folder();
    void check_load_steps_folder();

private:
    void update_dynamic_paths();
    static std::string sanitize_token(const std::string &value);

    std::string input_mesh_filepath;
    std::string mesh_filename_without_extension;
    std::string output_filename;

    std::string active_mesh;
    std::filesystem::path working_directory;
    std::filesystem::path model_data_filepath;
    bool custom_working_directory;
    Mode program_mode;
};

#endif // SESSION_H

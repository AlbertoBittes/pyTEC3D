#ifndef SESSION_H
#define SESSION_H

#include "../include/session.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <regex>
#include <string>
#include "../include/json.hpp"
#include <libmesh/libmesh.h>

using namespace libMesh;

class Session
{
public:
    bool interactive;
    nlohmann::json config;
    std::string output_folderpath;
    std::string mesh_filepath;
    std::filesystem::path output_path;
    std::filesystem::path mesh_path;
    std::filesystem::path final_output_path;
    std::string mesh_filename_without_extension;
    std::string output_filename;
    Session(const std::string &configFile);

    void read_config(const std::string &configFile);

    void set_output_filepath(unsigned int output_label);
    // private:
    //     void checkTransientAnalysisSettings()
    //     {
    //         if (config["Transient Analysis Settings"]["transient_analysis"].get<bool>())
    //         {
    //             Real analysis_time = config["Transient Analysis Settings"]["analysis time"].get<Real>();
    //             auto time_step_size = config["Transient Analysis Settings"]["time_step_size"];
    //             auto number_time_steps = config["Transient Analysis Settings"]["number_time_steps"];

    //             if (analysis_time <= 0)
    //             {
    //                 libmesh_error_msg("analysis_time must be > 0.");
    //             }
    //             if (time_step_size.is_null() && number_time_steps.is_null())
    //             {
    //                 libmesh_error_msg("At least one between time_step_size and number_time_steps must be defined.");
    //             }
    //             if (!time_step_size.is_null() && !number_time_steps.is_null())
    //             {
    //                 libmesh_error_msg("One between time_step_size and number_time_steps must be defined.");
    //             }

    //             if (time_step_size.is_number() && analysis_time <= time_step_size.get<Real>())
    //             {
    //                 libmesh_error_msg("analysis_time must be greater than time_step_size.");
    //             }

    //             if (time_step_size.is_number() && static_cast<int>(analysis_time) % static_cast<int>(time_step_size.get<Real>()) != 0)
    //             {
    //                 libmesh_error_msg("The remainder of the division analysis_time/time_step_size must be 0.");
    //             }
    //             if (!number_time_steps.is_null() && !number_time_steps.is_number_integer())
    //             {
    //                 libmesh_error_msg("number_time_steps must be an integer.");
    //             }
    //         }
    //     void checkFilePaths()
    //     {
    //         // Implement checks for filepath/folderpath variables to be strings and path-like
    //     }

    //     void checkStringValues()
    //     {
    //         // Implement checks for string values like approx_order, input_type, etc., to have proper values
    //     }

    //     void checkIntegerValues()
    //     {
    //         // List of integer variables to check
    //         std::vector<std::string> integerVariables = {
    //             "Finite Elements Settings/APEC_iterations",
    //             "Finite Elements Settings/quadrature_order"};

    //         for (const auto &varPath : integerVariables)
    //         {
    //             auto keys = nlohmann::json::parse("[\"" + std::regex_replace(varPath, std::regex("/"), "\",\"") + "\"]");
    //             nlohmann::json temp = config;
    //             bool found = true;
    //             for (const auto &key : keys)
    //             {
    //                 if (temp.find(key) == temp.end())
    //                 {
    //                     found = false;
    //                     break;
    //                 }
    //                 temp = temp[key];
    //             }
    //             if (!found || !temp.is_number_integer())
    //             {
    //                 libmesh_error_msg("Configuration error: " + varPath + " must be an integer value.");
    //             }
    //         }
    //     }
    //     void checkBooleanValues()
    //     {
    //         // List of boolean variables to check
    //         std::vector<std::string> booleanVariables = {
    //             "APEC_algorithm",
    //             "Transient Analysis Settings/transient_analysis",
    //             "nonlinear_analysis_ferrom"};

    //         for (const auto &varPath : booleanVariables)
    //         {
    //             auto keys = nlohmann::json::parse("[\"" + std::regex_replace(varPath, std::regex("/"), "\",\"") + "\"]");
    //             nlohmann::json temp = config;
    //             bool found = true;
    //             for (const auto &key : keys)
    //             {
    //                 if (temp.find(key) == temp.end())
    //                 {
    //                     found = false;
    //                     break;
    //                 }
    //                 temp = temp[key];
    //             }
    //             if (!found || !temp.is_boolean())
    //             {
    //                 libmesh_error_msg("Configuration error: " + varPath + " must be a boolean value.");
    //             }
    //         }
    //     }

    //     void checkMandatoryVariables()
    //     {
    //         std::vector<std::string> mandatoryPaths = {"Mesh Input Section/mesh_filepath", "Loads and BC Files Section/loads_folderpath", "Loads and BC Files Section/bc_folderpath", "Output Section/output_filepath"};
    //         for (const auto &path : mandatoryPaths)
    //         {
    //             auto keys = nlohmann::json::parse("[\"" + std::regex_replace(path, std::regex("/"), "\",\"") + "\"]");
    //             nlohmann::json temp = config;
    //             for (const auto &key : keys)
    //             {
    //                 if (temp.find(key) == temp.end() || !temp[key].is_string() || temp[key].get<std::string>().empty())
    //                 {
    //                     libmesh_error_msg("Mandatory variable " + path + " is missing or not a valid path.");
    //                 }
    //                 temp = temp[key];
    //             }
    //         }
    //     }
    // };
};

#endif // SESSION_H
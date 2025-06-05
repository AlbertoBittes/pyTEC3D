#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

namespace Globals
{
    // Enum for global variable names
    enum class GlobalVariableName
    {
        A_x,
        A_y,
        A_z,
        VOLT,
        Az,
        EMF,
        CURR,
        Js_x,
        Js_y,
        Js_z,
        B_x,
        B_y,
        B_z,
        J_x,
        J_y,
        J_z,
        UNKNOWN // Fallback for invalid strings
    };

    // Enum for API variable names
    enum class APIVariableName
    {
        A,
        VOLT,
        EMF,
        CURR,
        Js,
        B,
        J,
        UNKNOWN // Fallback for invalid strings
    };

    enum class VariableCategory
    {
        DOFVariables,
        LoadsVariables,
        NodalSolutionVariables,
        ElementSolutionVariables
        // NodalBCsVariables,
    };
    // Declare mappings
    extern const std::unordered_map<std::string, GlobalVariableName> STRING_TO_ENUM;
    extern const std::unordered_map<GlobalVariableName, std::string> ENUM_TO_STRING;
    extern const std::unordered_map<GlobalVariableName, std::pair<std::string, unsigned int>> GLOBAL_TO_API_VARIABLE_MAP;
    extern const std::unordered_map<APIVariableName, std::vector<std::string>> API_TO_GLOBAL_VARIABLE_MAP;
    extern const std::unordered_map<VariableCategory, std::vector<GlobalVariableName>> CATEGORY_TO_GLOBAL_VARIABLES_MAP;
    extern const std::unordered_map<GlobalVariableName, VariableCategory> GLOBAL_VARIABLE_TO_CATEGORY_MAP;

    // Declare global constants
    extern const double MU_0;

    // Utility functions
    GlobalVariableName string_to_global_variable_name(const std::string &name);
    std::string global_variable_name_to_string(GlobalVariableName name);
}

#endif // GLOBALS_H
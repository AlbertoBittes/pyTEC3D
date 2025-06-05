#include "globals.h"
#include <map>
#include <unordered_map>

namespace Globals
{
    // const std::vector<std::string> GLOBAL_VARIABLES_LABELS = {"A_x", "A_y", "A_z", "VOLT", "Az", "EMF", "CURR"};
    // const std::vector<std::string> API_VARIABLES_LABELS = {"A", "VOLT", "EMF", "CURR"};
    // const std::vector<unsigned int> API_VARIABLES_COMPONENTS = {3, 1, 1, 1};
    // const std::vector<std::string> GLOBAL_LOADS_LABELS = {"Js_x", "Js_y", "Js_z"};
    // const std::vector<std::string> API_LOADS_LABELS = {"Js"};
    // const std::vector<unsigned int> API_LOADS_COMPONENTS = {3};
    // const std::vector<std::string> GLOBAL_NODAL_DERIVED_VARIABLES_LABELS = {"B_x", "B_y", "B_z"};
    const std::unordered_map<std::string, GlobalVariableName> STRING_TO_ENUM = {
        {"A_x", GlobalVariableName::A_x},
        {"A_y", GlobalVariableName::A_y},
        {"A_z", GlobalVariableName::A_z},
        {"VOLT", GlobalVariableName::VOLT},
        {"Az", GlobalVariableName::Az},
        {"EMF", GlobalVariableName::EMF},
        {"CURR", GlobalVariableName::CURR},
        {"Js_x", GlobalVariableName::Js_x},
        {"Js_y", GlobalVariableName::Js_y},
        {"Js_z", GlobalVariableName::Js_z},
        {"B_x", GlobalVariableName::B_x},
        {"B_y", GlobalVariableName::B_y},
        {"B_z", GlobalVariableName::B_z},
        {"J_x", GlobalVariableName::J_x},
        {"J_y", GlobalVariableName::J_y},
        {"J_z", GlobalVariableName::J_z},
    };

    // Enum-to-string mapping
    const std::unordered_map<GlobalVariableName, std::string> ENUM_TO_STRING = {
        {GlobalVariableName::A_x, "A_x"},
        {GlobalVariableName::A_y, "A_y"},
        {GlobalVariableName::A_z, "A_z"},
        {GlobalVariableName::VOLT, "VOLT"},
        {GlobalVariableName::Az, "Az"},
        {GlobalVariableName::EMF, "EMF"},
        {GlobalVariableName::CURR, "CURR"},
        {GlobalVariableName::Js_x, "Js_x"},
        {GlobalVariableName::Js_y, "Js_y"},
        {GlobalVariableName::Js_z, "Js_z"},
        {GlobalVariableName::B_x, "B_x"},
        {GlobalVariableName::B_y, "B_y"},
        {GlobalVariableName::B_z, "B_z"},
        {GlobalVariableName::J_x, "J_x"},
        {GlobalVariableName::J_y, "J_y"},
        {GlobalVariableName::J_z, "J_z"},
    };
    // Use unordered_map for fast lookup
    const std::unordered_map<GlobalVariableName, std::pair<std::string, unsigned int>> GLOBAL_TO_API_VARIABLE_MAP = {
        {GlobalVariableName::A_x, {"A", 0}},
        {GlobalVariableName::A_y, {"A", 1}},
        {GlobalVariableName::A_z, {"A", 2}},
        {GlobalVariableName::VOLT, {"VOLT", 0}},
        {GlobalVariableName::Az, {"A", 2}}, // Az is the same as A_z
        {GlobalVariableName::EMF, {"EMF", 0}},
        {GlobalVariableName::CURR, {"CURR", 0}},
        {GlobalVariableName::Js_x, {"Js", 0}},
        {GlobalVariableName::Js_y, {"Js", 1}},
        {GlobalVariableName::Js_z, {"Js", 2}},
        {GlobalVariableName::B_x, {"B", 0}},
        {GlobalVariableName::B_y, {"B", 1}},
        {GlobalVariableName::B_z, {"B", 2}},
        {GlobalVariableName::J_x, {"J", 0}},
        {GlobalVariableName::J_y, {"J", 1}},
        {GlobalVariableName::J_z, {"J", 2}}};

    // Reverse mapping: API variable -> vector of global variable names (ordered by component)
    const std::unordered_map<APIVariableName, std::vector<std::string>> API_TO_GLOBAL_VARIABLE_MAP = {
        {APIVariableName::A, {"A_x", "A_y", "A_z"}},
        {APIVariableName::VOLT, {"VOLT"}},
        {APIVariableName::EMF, {"EMF"}},
        {APIVariableName::CURR, {"CURR"}},
        {APIVariableName::Js, {"Js_x", "Js_y", "Js_z"}},
        {APIVariableName::B, {"B_x", "B_y", "B_z"}},
        {APIVariableName::J, {"J_x", "J_y", "J_z"}}};

    const std::unordered_map<VariableCategory, std::vector<GlobalVariableName>> CATEGORY_TO_GLOBAL_VARIABLES_MAP = {
        {VariableCategory::DOFVariables, {GlobalVariableName::A_x, GlobalVariableName::A_y, GlobalVariableName::A_z, GlobalVariableName::VOLT, GlobalVariableName::Az, GlobalVariableName::EMF, GlobalVariableName::CURR}},
        {VariableCategory::LoadsVariables, {GlobalVariableName::Js_x, GlobalVariableName::Js_y, GlobalVariableName::Js_z}},
        {VariableCategory::NodalSolutionVariables, {GlobalVariableName::B_x, GlobalVariableName::B_y, GlobalVariableName::B_z}},
        {VariableCategory::ElementSolutionVariables, {GlobalVariableName::J_x, GlobalVariableName::J_y, GlobalVariableName::J_z}},
    };

    const std::unordered_map<GlobalVariableName, VariableCategory> GLOBAL_VARIABLE_TO_CATEGORY_MAP = {
        {GlobalVariableName::A_x, VariableCategory::DOFVariables},
        {GlobalVariableName::A_y, VariableCategory::DOFVariables},
        {GlobalVariableName::A_z, VariableCategory::DOFVariables},
        {GlobalVariableName::VOLT, VariableCategory::DOFVariables},
        {GlobalVariableName::Az, VariableCategory::DOFVariables},
        {GlobalVariableName::EMF, VariableCategory::DOFVariables},
        {GlobalVariableName::CURR, VariableCategory::DOFVariables},
        {GlobalVariableName::Js_x, VariableCategory::LoadsVariables},
        {GlobalVariableName::Js_y, VariableCategory::LoadsVariables},
        {GlobalVariableName::Js_z, VariableCategory::LoadsVariables},
        {GlobalVariableName::B_x, VariableCategory::NodalSolutionVariables},
        {GlobalVariableName::B_y, VariableCategory::NodalSolutionVariables},
        {GlobalVariableName::B_z, VariableCategory::NodalSolutionVariables},
        {GlobalVariableName::J_x, VariableCategory::ElementSolutionVariables},
        {GlobalVariableName::J_y, VariableCategory::ElementSolutionVariables},
        {GlobalVariableName::J_z, VariableCategory::ElementSolutionVariables},
    };
    // Convert string to GlobalVariableName
    inline GlobalVariableName string_to_global_variable_name(const std::string &name)
    {
        auto it = STRING_TO_ENUM.find(name);
        if (it != STRING_TO_ENUM.end())
        {
            return it->second;
        }
        return GlobalVariableName::UNKNOWN; // Fallback for invalid strings
    }

    // Convert GlobalVariableName to string
    inline std::string global_variable_name_to_string(GlobalVariableName name)
    {
        auto it = ENUM_TO_STRING.find(name);
        if (it != ENUM_TO_STRING.end())
        {
            return it->second;
        }
        return "UNKNOWN"; // Fallback for invalid enum values
    }

    const double MU_0 = 3.14159265358979323846 * 4.0 * 1e-7; // 4 * pi constant
}
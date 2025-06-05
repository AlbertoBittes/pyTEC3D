#include "globals.h"
#include <unordered_map>

namespace Globals
{
    const std::vector<std::string> GLOBAL_VARIABLES_LABELS = {"A_x", "A_y", "A_z", "VOLT", "Az", "EMF", "CURR"};
    const std::vector<std::string> API_VARIABLES_LABELS = {"A", "VOLT", "EMF", "CURR"};
    const std::vector<unsigned int> API_VARIABLES_COMPONENTS = {3, 1, 1, 1};
    const std::vector<std::string> GLOBAL_LOADS_LABELS = {"Js_x", "Js_y", "Js_z"};
    const std::vector<std::string> API_LOADS_LABELS = {"Js"};
    const std::vector<unsigned int> API_LOADS_COMPONENTS = {3};
    const std::vector<std::string> GLOBAL_NODAL_DERIVED_VARIABLES_LABELS = {"B_x", "B_y", "B_z"};

    // Use unordered_map for fast lookup
    const std::unordered_map<std::string, std::pair<std::string, unsigned int>> GLOBAL_TO_API_VARIABLE_MAP = {
        {"A_x", {"A", 0}},
        {"A_y", {"A", 1}},
        {"A_z", {"A", 2}},
        {"VOLT", {"VOLT", 0}},
        {"Az", {"A", 2}}, // If "Az" is a synonym for "A_z"
        {"EMF", {"EMF", 0}},
        {"CURR", {"CURR", 0}}};

    const std::unordered_map<std::string, std::pair<std::string, unsigned int>> GLOBAL_TO_API_LOAD_MAP = {
        {"Js_x", {"Js", 0}},
        {"Js_y", {"Js", 1}},
        {"Js_z", {"Js", 2}}};

    // Reverse mapping: API variable -> vector of global variable names (ordered by component)
    const std::unordered_map<std::string, std::vector<std::string>> API_TO_GLOBAL_VARIABLE_MAP = {
        {"A", {"A_x", "A_y", "A_z"}},
        {"VOLT", {"VOLT"}},
        {"EMF", {"EMF"}},
        {"CURR", {"CURR"}}};

    // Reverse mapping: API load -> vector of global load names (ordered by component)
    const std::unordered_map<std::string, std::vector<std::string>> API_TO_GLOBAL_LOAD_MAP = {
        {"Js", {"Js_x", "Js_y", "Js_z"}}};

    const double FOUR_PI = 3.14159265358979323846 * 4.0 * 1e-7; // 4 * pi constant
}
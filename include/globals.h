#ifndef GLOBALS_H
#define GLOBALS_H

#include <vector>
#include <string>

namespace Globals
{
    extern const std::vector<std::string> GLOBAL_VARIABLES_LABELS = {"A_x", "A_y", "A_z", "VOLT", "Az", "EMF", "CURR"};
    extern const std::vector<std::string> API_VARIABLES_LABELS = {"A", "VOLT", "EMF", "CURR"};
    extern const std::vector<unsigned int> API_VARIABLES_COMPONENTS = {3, 1, 1, 1};
    extern const std::vector<std::string> GLOBAL_LOADS_LABELS = {"Js"};
    extern const std::vector<std::string> API_LOADS_LABELS = {"Js"};
    extern const std::vector<unsigned int> API_LOADS_COMPONENTS = {3};
    extern const std::vector<std::string> GLOBAL_NODAL_DERIVED_VARIABLES_LABELS = {"B_x", "B_y", "B_z"};
}

#endif // GLOBALS_H
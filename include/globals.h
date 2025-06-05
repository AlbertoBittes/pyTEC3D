#ifndef GLOBALS_H
#define GLOBALS_H

#include <vector>
#include <string>

namespace Globals
{
    extern const std::vector<std::string> GLOBAL_VARIABLES_LABELS;
    extern const std::vector<std::string> API_VARIABLES_LABELS;
    extern const std::vector<unsigned int> API_VARIABLES_COMPONENTS;
    extern const std::vector<std::string> GLOBAL_LOADS_LABELS;
    extern const std::vector<std::string> API_LOADS_LABELS;
    extern const std::vector<unsigned int> API_LOADS_COMPONENTS;
    extern const std::vector<std::string> GLOBAL_NODAL_DERIVED_VARIABLES_LABELS;

    extern const std::unordered_map<std::string, std::pair<std::string, unsigned int>> GLOBAL_TO_API_VARIABLE_MAP;
    extern const std::unordered_map<std::string, std::pair<std::string, unsigned int>> GLOBAL_TO_API_LOAD_MAP;
    extern const std::unordered_map<std::string, std::vector<std::string>> API_TO_GLOBAL_VARIABLE_MAP;
    extern const std::unordered_map<std::string, std::vector<std::string>> API_TO_GLOBAL_LOAD_MAP;

    extern const double MU_0;
}

#endif // GLOBALS_H
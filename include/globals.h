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
}

#endif // GLOBALS_H
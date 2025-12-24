#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <map>
#include <vector>
#include <utility>
#include <variant>
#include <cstdint>
#include <libmesh/libmesh.h>
#include "nd_array.h"
// Some useful global variables
using namespace libMesh;
namespace Globals
{
    //
    // System type globals
    //
    enum SystemType
    {
        DOF_SOLUTION_SYSTEM = 0,
        NODAL_LOADS,
        ELEMENT_LOADS,
        NODAL_SOLUTION,
        ELEMENT_SOLUTION,
        NUM_SYSTEMS
    };

    // Global system info map
    extern const std::map<SystemType, std::string> SYSTEM_INFO_MAP;

    extern const libMesh::boundary_id_type FLUX_PARALLEL_ID_BC;
    extern const libMesh::boundary_id_type FLUX_NORMAL_ID_BC;

    // Declare physical constants
    extern const Real MU_0;

    extern const std::vector<std::string> VTK_DATA_ARRAY_NAMES;
}

#endif // GLOBALS_H

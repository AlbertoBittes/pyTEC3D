#include "globals.h"
#include <libmesh/libmesh.h>
#include <map>

using namespace libMesh;

namespace Globals
{
    const double MU_0 = pi * 4.0 * 1e-7;                // 4 * pi constant
    const boundary_id_type FLUX_PARALLEL_ID_BC = 32760; // Example ID for flux parallel boundary condition
    const boundary_id_type FLUX_NORMAL_ID_BC = 32761;   // Example ID for flux normal boundary condition

    const std::map<SystemType, std::string> SYSTEM_INFO_MAP = {
        {DOF_SOLUTION_SYSTEM, "ApecMLinearSystem"},
        {NODAL_LOADS, "NodalLoadsSystem"},
        {ELEMENT_LOADS, "ElementLoadsSystem"},
        {NODAL_SOLUTION, "NodalSolutionSystem"},
        {ELEMENT_SOLUTION, "ElementSolutionSystem"},
    };

    const std::vector<std::string> VTK_DATA_ARRAY_NAMES = {
        // "ansys_node_num",
        // "ansys_elem_num",
        "ansys_etype",
        "ansys_material_type",
        "ansys_real_constant",
        // "infinite_boundary",
        // "flux_normal",
        // "flux_parallel",
        // "periodic_boundary",
        "angles",
        "boundary_id_",
        "PLASMA",
        "VACUUM_CHAMBER",
        "VACUUM_CHAMBER_NODES",
        "VERTICAL_CONTROL_COILS",
        "RADIAL_CONTROL_COILS",
        "VACUUM_CHAMBER_BOUNDARY",
        "NEW_VACUUM_CHAMBER_BOUNDARY"};
}

#include "globals.h"
#include <libmesh/libmesh.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <map>
#include <unordered_map>

using namespace libMesh;

namespace Globals
{
    const std::map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP = {
        // DOF Variables
        {GlobalVariableName::A_x, {GlobalVariableName::A_x, APIVariableName::A, "A_x", "A", 0, SystemType::DOF_SOLUTION_SYSTEM, LAGRANGE, FIRST}},
        {GlobalVariableName::A_y, {GlobalVariableName::A_y, APIVariableName::A, "A_y", "A", 1, SystemType::DOF_SOLUTION_SYSTEM, LAGRANGE, FIRST}},
        {GlobalVariableName::A_z, {GlobalVariableName::A_z, APIVariableName::A, "A_z", "A", 2, SystemType::DOF_SOLUTION_SYSTEM, LAGRANGE, FIRST}},
        {GlobalVariableName::Az, {GlobalVariableName::Az, APIVariableName::Az, "Az", "Az", 0, SystemType::DOF_SOLUTION_SYSTEM, NEDELEC_ONE, FIRST}},
        {GlobalVariableName::VOLT, {GlobalVariableName::VOLT, APIVariableName::VOLT, "VOLT", "VOLT", 0, SystemType::DOF_SOLUTION_SYSTEM, LAGRANGE, FIRST}},
        {GlobalVariableName::EMF, {GlobalVariableName::EMF, APIVariableName::EMF, "EMF", "EMF", 0, SystemType::DOF_SOLUTION_SYSTEM, LAGRANGE, FIRST}},
        {GlobalVariableName::V, {GlobalVariableName::V, APIVariableName::V, "V", "V", 0, SystemType::DOF_SOLUTION_SYSTEM, LAGRANGE, FIRST}},

        // Nodal Loads Variables
        {GlobalVariableName::CSG_x, {GlobalVariableName::CSG_x, APIVariableName::CSG, "CSG_x", "CSG", 0, SystemType::NODAL_LOADS, LAGRANGE, FIRST}},
        {GlobalVariableName::CSG_y, {GlobalVariableName::CSG_y, APIVariableName::CSG, "CSG_y", "CSG", 1, SystemType::NODAL_LOADS, LAGRANGE, FIRST}},
        {GlobalVariableName::CSG_z, {GlobalVariableName::CSG_z, APIVariableName::CSG, "CSG_z", "CSG", 2, SystemType::NODAL_LOADS, LAGRANGE, FIRST}},
        {GlobalVariableName::AMPS, {GlobalVariableName::AMPS, APIVariableName::AMPS, "AMPS", "AMPS", 0, SystemType::NODAL_LOADS, LAGRANGE, FIRST}},
        {GlobalVariableName::CURT, {GlobalVariableName::CURT, APIVariableName::CURT, "CURT", "CURT", 0, SystemType::NODAL_LOADS, LAGRANGE, FIRST}},

        // Element Loads Variables
        {GlobalVariableName::Js_x, {GlobalVariableName::Js_x, APIVariableName::Js, "Js_x", "Js", 0, SystemType::ELEMENT_LOADS, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Js_y, {GlobalVariableName::Js_y, APIVariableName::Js, "Js_y", "Js", 1, SystemType::ELEMENT_LOADS, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Js_z, {GlobalVariableName::Js_z, APIVariableName::Js, "Js_z", "Js", 2, SystemType::ELEMENT_LOADS, MONOMIAL, CONSTANT}},

        // Nodal Solution Variables
        {GlobalVariableName::B_x, {GlobalVariableName::B_x, APIVariableName::B, "B_x", "B", 0, SystemType::NODAL_SOLUTION, L2_LAGRANGE, FIRST}},
        {GlobalVariableName::B_y, {GlobalVariableName::B_y, APIVariableName::B, "B_y", "B", 1, SystemType::NODAL_SOLUTION, L2_LAGRANGE, FIRST}},
        {GlobalVariableName::B_z, {GlobalVariableName::B_z, APIVariableName::B, "B_z", "B", 2, SystemType::NODAL_SOLUTION, L2_LAGRANGE, FIRST}},
        {GlobalVariableName::H_x, {GlobalVariableName::H_x, APIVariableName::H, "H_x", "H", 0, SystemType::NODAL_SOLUTION, L2_LAGRANGE, FIRST}},
        {GlobalVariableName::H_y, {GlobalVariableName::H_y, APIVariableName::H, "H_y", "H", 1, SystemType::NODAL_SOLUTION, L2_LAGRANGE, FIRST}},
        {GlobalVariableName::H_z, {GlobalVariableName::H_z, APIVariableName::H, "H_z", "H", 2, SystemType::NODAL_SOLUTION, L2_LAGRANGE, FIRST}},
        {GlobalVariableName::PSI, {GlobalVariableName::PSI, APIVariableName::PSI, "PSI", "PSI", 0, SystemType::NODAL_SOLUTION, LAGRANGE, FIRST}},

        // Element Solution Variables
        {GlobalVariableName::Jt_x, {GlobalVariableName::Jt_x, APIVariableName::Jt, "Jt_x", "Jt", 0, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Jt_y, {GlobalVariableName::Jt_y, APIVariableName::Jt, "Jt_y", "Jt", 1, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Jt_z, {GlobalVariableName::Jt_z, APIVariableName::Jt, "Jt_z", "Jt", 2, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Be_x, {GlobalVariableName::Be_x, APIVariableName::Be, "Be_x", "Be", 0, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Be_y, {GlobalVariableName::Be_y, APIVariableName::Be, "Be_y", "Be", 1, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::Be_z, {GlobalVariableName::Be_z, APIVariableName::Be, "Be_z", "Be", 2, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::He_x, {GlobalVariableName::He_x, APIVariableName::He, "He_x", "He", 0, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::He_y, {GlobalVariableName::He_y, APIVariableName::He, "He_y", "He", 1, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::He_z, {GlobalVariableName::He_z, APIVariableName::He, "He_z", "He", 2, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
        {GlobalVariableName::PSIe, {GlobalVariableName::PSIe, APIVariableName::PSIe, "PSIe", "PSIe", 0, SystemType::ELEMENT_SOLUTION, MONOMIAL, CONSTANT}},
    };

    const double MU_0 = pi * 4.0 * 1e-7;                // 4 * pi constant
    const boundary_id_type FLUX_PARALLEL_ID_BC = 32760; // Example ID for flux parallel boundary condition
    const boundary_id_type FLUX_NORMAL_ID_BC = 32761;   // Example ID for flux normal boundary condition

    const std::map<SystemType, SystemInfo> SYSTEM_INFO_MAP = {
        {DOF_SOLUTION_SYSTEM, {DOF_SOLUTION_SYSTEM, "ApecMLinearSystem"}},
        {NODAL_LOADS, {NODAL_LOADS, "NodalLoadsSystem"}},
        {ELEMENT_LOADS, {ELEMENT_LOADS, "ElementLoadsSystem"}},
        {NODAL_SOLUTION, {NODAL_SOLUTION, "NodalSolutionSystem"}},
        {ELEMENT_SOLUTION, {ELEMENT_SOLUTION, "ElementSolutionSystem"}},
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
        // "MidsideNodes",
        "BOUNDARY_ID_",
        "PLASMA",
        "VACUUM_CHAMBER",
        "VACUUM_CHAMBER_NODES",
        "VERTICAL_CONTROL_COILS",
        "RADIAL_CONTROL_COILS",
        "VACUUM_CHAMBER_BOUNDARY",
        "NEW_VACUUM_CHAMBER_BOUNDARY"};
}

#include "globals.h"
#include <libmesh/libmesh.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <map>
#include <unordered_map>

using namespace libMesh;

namespace Globals
{
    const std::unordered_map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP = {
        // DOF Variables
        {GlobalVariableName::A_x, {GlobalVariableName::A_x, APIVariableName::A, VariableCategory::DOFVariables, "A_x", "A", LAGRANGE, FIRST, 0, true, 0}},
        {GlobalVariableName::A_y, {GlobalVariableName::A_y, APIVariableName::A, VariableCategory::DOFVariables, "A_y", "A", LAGRANGE, FIRST, 1, true, 0}},
        {GlobalVariableName::A_z, {GlobalVariableName::A_z, APIVariableName::A, VariableCategory::DOFVariables, "A_z", "A", LAGRANGE, FIRST, 2, true, 0}},
        {GlobalVariableName::VOLT, {GlobalVariableName::VOLT, APIVariableName::VOLT, VariableCategory::DOFVariables, "VOLT", "VOLT", LAGRANGE, FIRST, 0, false, 0}},
        {GlobalVariableName::Az, {GlobalVariableName::Az, APIVariableName::Az, VariableCategory::DOFVariables, "Az", "Az", LAGRANGE, FIRST, 0, true, 0}},
        {GlobalVariableName::EMF, {GlobalVariableName::EMF, APIVariableName::EMF, VariableCategory::DOFVariables, "EMF", "EMF", LAGRANGE, FIRST, 0, false, 0}},
        {GlobalVariableName::CURR, {GlobalVariableName::CURR, APIVariableName::CURR, VariableCategory::DOFVariables, "CURR", "CURR", LAGRANGE, FIRST, 0, false, 0}},

        // Nodal Loads Variables
        {GlobalVariableName::AMPS, {GlobalVariableName::AMPS, APIVariableName::AMPS, VariableCategory::NodalLoadsVariables, "AMPS", "AMPS", LAGRANGE, FIRST, 0, true, 0}},
        {GlobalVariableName::CSGZ, {GlobalVariableName::CSGZ, APIVariableName::CSGZ, VariableCategory::NodalLoadsVariables, "CSGZ", "CSGZ", LAGRANGE, FIRST, 0, true, 0}},

        // Element Loads Variables
        {GlobalVariableName::Js_x, {GlobalVariableName::Js_x, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_x", "Js", MONOMIAL, CONSTANT, 0, true, 0}},
        {GlobalVariableName::Js_y, {GlobalVariableName::Js_y, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_y", "Js", MONOMIAL, CONSTANT, 1, true, 0}},
        {GlobalVariableName::Js_z, {GlobalVariableName::Js_z, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_z", "Js", MONOMIAL, CONSTANT, 2, true, 0}},

        // Nodal Solution Variables
        {GlobalVariableName::B_x, {GlobalVariableName::B_x, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_x", "B", L2_LAGRANGE, FIRST, 0, true, 0}},
        {GlobalVariableName::B_y, {GlobalVariableName::B_y, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_y", "B", L2_LAGRANGE, FIRST, 1, true, 0}},
        {GlobalVariableName::B_z, {GlobalVariableName::B_z, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_z", "B", L2_LAGRANGE, FIRST, 2, true, 0}},

        // Element Solution Variables
        {GlobalVariableName::J_x, {GlobalVariableName::J_x, APIVariableName::J, VariableCategory::ElementSolutionVariables, "J_x", "J", MONOMIAL, CONSTANT, 0, true, 0}},
        {GlobalVariableName::J_y, {GlobalVariableName::J_y, APIVariableName::J, VariableCategory::ElementSolutionVariables, "J_y", "J", MONOMIAL, CONSTANT, 1, true, 0}},
        {GlobalVariableName::J_z, {GlobalVariableName::J_z, APIVariableName::J, VariableCategory::ElementSolutionVariables, "J_z", "J", MONOMIAL, CONSTANT, 2, true, 0}}};

    const double MU_0 = pi * 4.0 * 1e-7;                // 4 * pi constant
    const boundary_id_type FLUX_PARALLEL_ID_BC = 32760; // Example ID for flux parallel boundary condition
    const boundary_id_type FLUX_NORMAL_ID_BC = 32761;   // Example ID for flux normal boundary condition

    const std::unordered_map<SystemType, SystemInfo> SYSTEM_INFO_MAP = {
        {MAIN, {MAIN, "ApecMLinearSystem", Globals::VariableCategory::DOFVariables}},
        {NODAL_LOADS, {NODAL_LOADS, "NodalLoadsSystem", Globals::VariableCategory::NodalLoadsVariables}},
        {ELEMENT_LOADS, {ELEMENT_LOADS, "ElementLoadsSystem", Globals::VariableCategory::ElementLoadsVariables}},
        {NODAL_EXPLICIT, {NODAL_EXPLICIT, "NodalSolutionSystem", Globals::VariableCategory::NodalSolutionVariables}},
        {ELEMENT_EXPLICIT, {ELEMENT_EXPLICIT, "ElementSolutionSystem", Globals::VariableCategory::ElementSolutionVariables}},
    };

    const std::unordered_map<BoundaryConditionType, BoundaryConditionInfo> BOUNDARY_CONDITION_INFO_MAP = {
        {BoundaryConditionType::INFINITE_BOUNDARY, {BoundaryConditionType::INFINITE_BOUNDARY, "infinite_boundary", 0}},
        {BoundaryConditionType::FLUX_NORMAL, {BoundaryConditionType::FLUX_NORMAL, "flux_normal", FLUX_NORMAL_ID_BC}},
        {BoundaryConditionType::FLUX_PARALLEL, {BoundaryConditionType::FLUX_PARALLEL, "flux_parallel", FLUX_PARALLEL_ID_BC}},
        {BoundaryConditionType::UNKNOWN, {BoundaryConditionType::UNKNOWN, "UNKNOWN", 0}}};
}
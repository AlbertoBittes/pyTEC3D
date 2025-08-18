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
        {GlobalVariableName::A_x, {GlobalVariableName::A_x, APIVariableName::A, VariableCategory::DOFVariables, "A_x", "A", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::A_y, {GlobalVariableName::A_y, APIVariableName::A, VariableCategory::DOFVariables, "A_y", "A", LAGRANGE, FIRST, 1, 0}},
        {GlobalVariableName::A_z, {GlobalVariableName::A_z, APIVariableName::A, VariableCategory::DOFVariables, "A_z", "A", LAGRANGE, FIRST, 2, 0}},
        {GlobalVariableName::VOLT, {GlobalVariableName::VOLT, APIVariableName::VOLT, VariableCategory::DOFVariables, "VOLT", "VOLT", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::Az, {GlobalVariableName::Az, APIVariableName::Az, VariableCategory::DOFVariables, "Az", "Az", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::EMF, {GlobalVariableName::EMF, APIVariableName::EMF, VariableCategory::DOFVariables, "EMF", "EMF", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::CURR, {GlobalVariableName::CURR, APIVariableName::CURR, VariableCategory::DOFVariables, "CURR", "CURR", LAGRANGE, FIRST, 0, 0}},

        // Nodal Loads Variables
        {GlobalVariableName::CSG_x, {GlobalVariableName::CSG_x, APIVariableName::CSG, VariableCategory::NodalLoadsVariables, "CSG_x", "CSG", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::CSG_y, {GlobalVariableName::CSG_y, APIVariableName::CSG, VariableCategory::NodalLoadsVariables, "CSG_y", "CSG", LAGRANGE, FIRST, 1, 0}},
        {GlobalVariableName::CSG_z, {GlobalVariableName::CSG_z, APIVariableName::CSG, VariableCategory::NodalLoadsVariables, "CSG_z", "CSG", LAGRANGE, FIRST, 2, 0}},
        {GlobalVariableName::AMPS, {GlobalVariableName::AMPS, APIVariableName::AMPS, VariableCategory::NodalLoadsVariables, "AMPS", "AMPS", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::VLTG, {GlobalVariableName::VLTG, APIVariableName::VLTG, VariableCategory::NodalLoadsVariables, "VLTG", "VLTG", LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::CURT, {GlobalVariableName::CURT, APIVariableName::CURT, VariableCategory::NodalLoadsVariables, "CURT", "CURT", LAGRANGE, FIRST, 0, 0}},

        // Element Loads Variables
        {GlobalVariableName::Js_x, {GlobalVariableName::Js_x, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_x", "Js", MONOMIAL, CONSTANT, 0, 0}},
        {GlobalVariableName::Js_y, {GlobalVariableName::Js_y, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_y", "Js", MONOMIAL, CONSTANT, 1, 0}},
        {GlobalVariableName::Js_z, {GlobalVariableName::Js_z, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_z", "Js", MONOMIAL, CONSTANT, 2, 0}},

        // Nodal Solution Variables
        {GlobalVariableName::B_x, {GlobalVariableName::B_x, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_x", "B", L2_LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::B_y, {GlobalVariableName::B_y, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_y", "B", L2_LAGRANGE, FIRST, 1, 0}},
        {GlobalVariableName::B_z, {GlobalVariableName::B_z, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_z", "B", L2_LAGRANGE, FIRST, 2, 0}},
        {GlobalVariableName::H_x, {GlobalVariableName::H_x, APIVariableName::H, VariableCategory::NodalSolutionVariables, "H_x", "H", L2_LAGRANGE, FIRST, 0, 0}},
        {GlobalVariableName::H_y, {GlobalVariableName::H_y, APIVariableName::H, VariableCategory::NodalSolutionVariables, "H_y", "H", L2_LAGRANGE, FIRST, 1, 0}},
        {GlobalVariableName::H_z, {GlobalVariableName::H_z, APIVariableName::H, VariableCategory::NodalSolutionVariables, "H_z", "H", L2_LAGRANGE, FIRST, 2, 0}},

        // Element Solution Variables
        {GlobalVariableName::Jt_x, {GlobalVariableName::Jt_x, APIVariableName::Jt, VariableCategory::ElementSolutionVariables, "Jt_x", "Jt", MONOMIAL, CONSTANT, 0, 0}},
        {GlobalVariableName::Jt_y, {GlobalVariableName::Jt_y, APIVariableName::Jt, VariableCategory::ElementSolutionVariables, "Jt_y", "Jt", MONOMIAL, CONSTANT, 1, 0}},
        {GlobalVariableName::Jt_z, {GlobalVariableName::Jt_z, APIVariableName::Jt, VariableCategory::ElementSolutionVariables, "Jt_z", "Jt", MONOMIAL, CONSTANT, 2, 0}},
        {GlobalVariableName::Be_x, {GlobalVariableName::Be_x, APIVariableName::Be, VariableCategory::ElementSolutionVariables, "Be_x", "Be", MONOMIAL, CONSTANT, 0, 0}},
        {GlobalVariableName::Be_y, {GlobalVariableName::Be_y, APIVariableName::Be, VariableCategory::ElementSolutionVariables, "Be_y", "Be", MONOMIAL, CONSTANT, 1, 0}},
        {GlobalVariableName::Be_z, {GlobalVariableName::Be_z, APIVariableName::Be, VariableCategory::ElementSolutionVariables, "Be_z", "Be", MONOMIAL, CONSTANT, 2, 0}},
        {GlobalVariableName::He_x, {GlobalVariableName::He_x, APIVariableName::He, VariableCategory::ElementSolutionVariables, "He_x", "He", MONOMIAL, CONSTANT, 0, 0}},
        {GlobalVariableName::He_y, {GlobalVariableName::He_y, APIVariableName::He, VariableCategory::ElementSolutionVariables, "He_y", "He", MONOMIAL, CONSTANT, 1, 0}},
        {GlobalVariableName::He_z, {GlobalVariableName::He_z, APIVariableName::He, VariableCategory::ElementSolutionVariables, "He_z", "He", MONOMIAL, CONSTANT, 2, 0}}};

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
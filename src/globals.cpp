#include "globals.h"
#include <libmesh/libmesh.h>
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
        {GlobalVariableName::Az, {GlobalVariableName::Az, APIVariableName::A, VariableCategory::DOFVariables, "Az", "Az", LAGRANGE, FIRST, 0, true, 0}},
        {GlobalVariableName::EMF, {GlobalVariableName::EMF, APIVariableName::EMF, VariableCategory::DOFVariables, "EMF", "EMF", LAGRANGE, FIRST, 0, false, 0}},
        {GlobalVariableName::CURR, {GlobalVariableName::CURR, APIVariableName::CURR, VariableCategory::DOFVariables, "CURR", "CURR", LAGRANGE, FIRST, 0, false, 0}},

        // Element Loads Variables
        {GlobalVariableName::Js_x, {GlobalVariableName::Js_x, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_x", "Js", MONOMIAL, CONSTANT, 0, true, 0}},
        {GlobalVariableName::Js_y, {GlobalVariableName::Js_y, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_y", "Js", MONOMIAL, CONSTANT, 1, true, 0}},
        {GlobalVariableName::Js_z, {GlobalVariableName::Js_z, APIVariableName::Js, VariableCategory::ElementLoadsVariables, "Js_z", "Js", MONOMIAL, CONSTANT, 2, true, 0}},

        // Nodal Loads Variables
        {GlobalVariableName::AMPS, {GlobalVariableName::AMPS, APIVariableName::AMPS, VariableCategory::NodalLoadsVariables, "AMPS", "AMPS", LAGRANGE, FIRST, 0, true, 0}},

        // Nodal Solution Variables
        {GlobalVariableName::B_x, {GlobalVariableName::B_x, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_x", "B", L2_LAGRANGE, FIRST, 0, true, 0}},
        {GlobalVariableName::B_y, {GlobalVariableName::B_y, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_y", "B", L2_LAGRANGE, FIRST, 1, true, 0}},
        {GlobalVariableName::B_z, {GlobalVariableName::B_z, APIVariableName::B, VariableCategory::NodalSolutionVariables, "B_z", "B", L2_LAGRANGE, FIRST, 2, true, 0}},

        // Element Solution Variables
        {GlobalVariableName::J_x, {GlobalVariableName::J_x, APIVariableName::J, VariableCategory::ElementSolutionVariables, "J_x", "J", MONOMIAL, CONSTANT, 0, true, 0}},
        {GlobalVariableName::J_y, {GlobalVariableName::J_y, APIVariableName::J, VariableCategory::ElementSolutionVariables, "J_y", "J", MONOMIAL, CONSTANT, 1, true, 0}},
        {GlobalVariableName::J_z, {GlobalVariableName::J_z, APIVariableName::J, VariableCategory::ElementSolutionVariables, "J_z", "J", MONOMIAL, CONSTANT, 2, true, 0}},
    };

    const double MU_0 = libMesh::pi * 4.0 * 1e-7;                // 4 * pi constant
    const libMesh::boundary_id_type FLUX_PARALLEL_ID_BC = 32760; // Example ID for flux parallel boundary condition
    const libMesh::boundary_id_type FLUX_NORMAL_ID_BC = 32761;   // Example ID for flux normal boundary condition

    const std::unordered_map<SystemType, SystemInfo> SYSTEM_INFO_MAP = {
        {MAIN, {MAIN, "ApecMLinearSystem", Globals::VariableCategory::DOFVariables}},
        {NODAL_EXPLICIT, {NODAL_EXPLICIT, "NodalSolutionSystem", Globals::VariableCategory::NodalSolutionVariables}},
        {ELEMENT_EXPLICIT, {ELEMENT_EXPLICIT, "ElementSolutionSystem", Globals::VariableCategory::ElementSolutionVariables}},
        {ELEMENT_LOADS, {ELEMENT_LOADS, "ElementLoadsSystem", Globals::VariableCategory::ElementLoadsVariables}},
        {NODAL_LOADS, {NODAL_LOADS, "NodalLoadsSystem", Globals::VariableCategory::NodalLoadsVariables}},
    };

    const std::unordered_map<BoundaryConditionType, BoundaryConditionInfo> BOUNDARY_CONDITION_INFO_MAP = {
        {BoundaryConditionType::INFINITE_BOUNDARY, {BoundaryConditionType::INFINITE_BOUNDARY, "infinite_boundary", 0}},
        {BoundaryConditionType::FLUX_NORMAL, {BoundaryConditionType::FLUX_NORMAL, "flux_normal", FLUX_NORMAL_ID_BC}},
        {BoundaryConditionType::FLUX_PARALLEL, {BoundaryConditionType::FLUX_PARALLEL, "flux_parallel", FLUX_PARALLEL_ID_BC}},
        {BoundaryConditionType::UNKNOWN, {BoundaryConditionType::UNKNOWN, "UNKNOWN", 0}}};

    template <typename MapType>
    const typename MapType::mapped_type &get_info_by_key(const MapType &map, const typename MapType::key_type &key)
    {
        auto it = map.find(key);
        if (it == map.end())
            libmesh_error_msg("Key not found in info map");
        return it->second;
    }

    // Usage example
    // auto info = get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, Globals::GlobalVariableName::A_x);
    // auto sysinfo = get_info_by_key(Globals::SYSTEM_INFO_MAP, Globals::SystemType::MAIN);
    // auto bcinfo = get_info_by_key(Globals::BOUNDARY_CONDITION_INFO_MAP, Globals::BoundaryConditionType::FLUX_NORMAL);

    template <typename MapType, typename Predicate>
    std::vector<typename MapType::key_type> get_keys_by_predicate(const MapType &map, Predicate pred)
    {
        std::vector<typename MapType::key_type> result;
        for (const auto &pair : map)
            if (pred(pair.second))
                result.push_back(pair.first);
        return result;
    }
    // Usage example
    // auto keys = get_keys_by_predicate(
    //     Globals::GLOBAL_VARIABLE_INFO_MAP,
    //     [](const Globals::VariableInfo &info)
    //     { return info.api_name == Globals::APIVariableName::A; });

    // // Get all SystemType for a given VariableCategory
    // auto sys_types = get_keys_by_predicate(
    //     Globals::SYSTEM_INFO_MAP,
    //     [](const Globals::SystemInfo &info)
    //     { return info.category == Globals::VariableCategory::ElementLoadsVariables; });

    template <typename MapType, typename Predicate>
    typename MapType::key_type get_key_by_predicate_one(const MapType &map, Predicate pred)
    {
        auto keys = get_keys_by_predicate(map, pred);
        if (keys.size() == 1)
            return keys[0];
        libmesh_error_msg("get_key_by_predicate_one: mapping is not one-to-one or key not found");
    }
    // Usage example
    // // Get GlobalVariableName for API label and component
    // auto key = get_key_by_predicate_one(
    //     Globals::GLOBAL_VARIABLE_INFO_MAP,
    //     [](const Globals::VariableInfo &info)
    //     {
    //         return info.api_label == "A" && info.component_index == 2;
    //     });
}
#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>
#include <libmesh/libmesh.h>
#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>
#include <libmesh/fe.h>
// Some useful global variables
using namespace libMesh;
namespace Globals
{
    //
    // Variables globals
    //
    // Enum for global variable names
    enum class GlobalVariableName
    {
        A_x,
        A_y,
        A_z,
        Az,
        VOLT,
        EMF,
        CURR,
        Js_x,
        Js_y,
        Js_z,
        CSG_x,
        CSG_y,
        CSG_z,
        AMPS,
        VLTG,
        CURT,
        B_x,
        B_y,
        B_z,
        H_x,
        H_y,
        H_z,
        Jt_x,
        Jt_y,
        Jt_z,
        Be_x,
        Be_y,
        Be_z,
        He_x,
        He_y,
        He_z,
        UNKNOWN // Fallback for invalid strings
    };

    // Enum for API variable names
    enum class APIVariableName
    {
        A,
        Az,
        VOLT,
        EMF,
        CURR,
        CSG,
        AMPS,
        VLTG,
        CURT,
        Js,
        B,
        H,
        Jt,
        Be,
        He,
        UNKNOWN // Fallback for invalid strings
    };

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

    struct VariableInfo
    {
        GlobalVariableName global_name;
        APIVariableName api_name;
        std::string global_label;
        std::string api_label;
        unsigned int component_index;
        SystemType system_type;
        FEFamily fe_family;
        Order order;
        // Add other fields as needed
    };

    // Global variable information map
    extern const std::map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP;

    struct SystemInfo
    {
        SystemType type;
        std::string name;
        // Add other fields as needed
    };
    // Global system info map
    extern const std::map<SystemType, SystemInfo> SYSTEM_INFO_MAP;

    //
    // BCs globals infomrations
    //
    enum class BoundaryConditionType
    {
        INFINITE_BOUNDARY,
        FLUX_NORMAL,
        FLUX_PARALLEL,
        PERIODIC_BOUNDARY,
        UNKNOWN // Fallback for invalid strings
    };

    struct BoundaryConditionInfo
    {
        BoundaryConditionType array_enum;
        std::string label;
        libMesh::boundary_id_type id;
        // Add other fields as needed
    };
    extern const libMesh::boundary_id_type FLUX_PARALLEL_ID_BC;
    extern const libMesh::boundary_id_type FLUX_NORMAL_ID_BC;
    // Centralized boundary condition info map
    extern const std::map<BoundaryConditionType, BoundaryConditionInfo> BOUNDARY_CONDITION_INFO_MAP;

    // Declare physical constants
    extern const double MU_0;

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
    // auto sysinfo = get_info_by_key(Globals::SYSTEM_INFO_MAP, Globals::SystemType::DOF_SOLUTION_SYSTEM);
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

#endif // GLOBALS_H
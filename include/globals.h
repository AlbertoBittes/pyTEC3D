#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>
#include <unordered_map>
#include <vector>
#include <utility>
#include <libmesh/libmesh.h>
#include <libmesh/fe.h>
// Some useful global variables
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
        VOLT,
        Az,
        EMF,
        CURR,
        Js_x,
        Js_y,
        Js_z,
        AMPS,
        B_x,
        B_y,
        B_z,
        J_x,
        J_y,
        J_z,
        UNKNOWN // Fallback for invalid strings
    };

    // Enum for API variable names
    enum class APIVariableName
    {
        A,
        VOLT,
        EMF,
        CURR,
        Js,
        AMPS,
        B,
        J,
        UNKNOWN // Fallback for invalid strings
    };

    enum class VariableCategory
    {
        DOFVariables,
        ElementLoadsVariables,
        NodalLoadsVariables,
        NodalSolutionVariables,
        ElementSolutionVariables
    };

    struct VariableInfo
    {
        GlobalVariableName global_name;
        APIVariableName api_name;
        VariableCategory category;
        std::string global_label;
        std::string api_label;
        FEFamily fe_family;
        Order order;
        unsigned int component_index;
        bool is_active;
        unsigned int system_var_num;
        // Add other fields as needed
    };

    // Global variable information map
    extern const std::unordered_map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP;

    //
    // System type globals
    //
    enum SystemType
    {
        MAIN = 0,
        NODAL_EXPLICIT,
        ELEMENT_EXPLICIT,
        ELEMENT_LOADS,
        NODAL_LOADS,
        NUM_SYSTEMS
    };

    struct SystemInfo
    {
        SystemType type;
        std::string name;
        Globals::VariableCategory category;
        // Add other fields as needed
    };
    // Global system info map
    extern const std::unordered_map<SystemType, SystemInfo> SYSTEM_INFO_MAP;

    //
    // BCs globals infomrations
    //
    enum class BoundaryConditionType
    {
        INFINITE_BOUNDARY,
        FLUX_NORMAL,
        FLUX_PARALLEL,
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
    extern const std::unordered_map<BoundaryConditionType, BoundaryConditionInfo> BOUNDARY_CONDITION_INFO_MAP;

    // Declare physical constants
    extern const double MU_0;

    template <typename MapType>
    const typename MapType::mapped_type &get_info_by_key(const MapType &map, const typename MapType::key_type &key);
    template <typename MapType, typename Predicate>
    std::vector<typename MapType::key_type> get_keys_by_predicate(const MapType &map, Predicate pred);
    template <typename MapType, typename Predicate>
    typename MapType::key_type get_key_by_predicate_one(const MapType &map, Predicate pred);

}

#endif // GLOBALS_H
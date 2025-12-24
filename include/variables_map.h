#ifndef VARIABLES_MAP_H
#define VARIABLES_MAP_H

#include <map>
#include <string>
#include <vector>

#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>

#include "globals.h"

class ElementType;
enum PhysicalRegionEnum : unsigned int;

enum class GlobalVariableName
{
    A_x,
    A_y,
    A_z,
    Az,
    VOLT,
    EMF,
    V,
    Js_x,
    Js_y,
    Js_z,
    CSG_x,
    CSG_y,
    CSG_z,
    AMPS,
    CURT,
    B_x,
    B_y,
    B_z,
    H_x,
    H_y,
    H_z,
    PSI,
    Jt_x,
    Jt_y,
    Jt_z,
    Be_x,
    Be_y,
    Be_z,
    He_x,
    He_y,
    He_z,
    PSIe,
    UNKNOWN // Fallback for invalid strings
};

enum class APIVariableName
{
    A,
    Az,
    VOLT,
    EMF,
    V,
    CSG,
    AMPS,
    CURT,
    Js,
    B,
    H,
    PSI,
    Jt,
    Be,
    He,
    PSIe,
    UNKNOWN // Fallback for invalid strings
};

struct VariableInfo
{
    GlobalVariableName global_name;
    APIVariableName api_name;
    std::string global_label;
    std::string api_label;
    unsigned int component_index;
    Globals::SystemType system_type;
};

class VariablesMap
{
public:
    static const std::map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP;

    template <typename MapType>
    static const typename MapType::mapped_type &get_info_by_key(const MapType &map,
                                                                const typename MapType::key_type &key)
    {
        auto it = map.find(key);
        if (it == map.end())
            libmesh_error_msg("Key not found in info map");
        return it->second;
    }

    template <typename MapType, typename Predicate>
    static std::vector<typename MapType::key_type> get_keys_by_predicate(const MapType &map, Predicate pred)
    {
        std::vector<typename MapType::key_type> result;
        for (const auto &pair : map)
            if (pred(pair.second))
                result.push_back(pair.first);
        return result;
    }

    template <typename MapType, typename Predicate>
    static typename MapType::key_type get_key_by_predicate_one(const MapType &map, Predicate pred)
    {
        auto keys = get_keys_by_predicate(map, pred);
        if (keys.size() == 1)
            return keys[0];
        libmesh_error_msg("get_key_by_predicate_one: mapping is not one-to-one or key not found");
    }

    static const VariableInfo &get_info_by_key(GlobalVariableName key);
    static GlobalVariableName get_key_by_global_label(const std::string &global_label);
    static GlobalVariableName get_key_by_api_component(APIVariableName api_name, unsigned int component_index);
    static std::vector<GlobalVariableName> get_keys_by_api_label(const std::string &api_label);

    static std::vector<GlobalVariableName> get_variables_for_element(const ElementType &element);
    static std::vector<GlobalVariableName> get_active_variables_for_element(const ElementType &element,
                                                                            Globals::SystemType system_type);

    static FEFamily get_variable_fe_family(const ElementType &element,
                                           GlobalVariableName variable);
    static Order get_variable_order(const ElementType &element,
                                    GlobalVariableName variable);
    static unsigned int get_variable_physical_components(const ElementType &element,
                                                         GlobalVariableName variable);
    static bool is_variable_scalar(const ElementType &element,
                                   GlobalVariableName variable);
    static bool is_variable_vector(const ElementType &element,
                                   GlobalVariableName variable);

private:
    static std::vector<GlobalVariableName> add_av_variables();
    static std::vector<GlobalVariableName> add_av3d_variables();
    static std::vector<GlobalVariableName> add_solid97_variables(PhysicalRegionEnum physical_region);
    static std::vector<GlobalVariableName> add_solid236_variables(PhysicalRegionEnum physical_region);
    static std::vector<GlobalVariableName> add_plane233_variables(const std::vector<unsigned int> &keyopts);
};

#endif // VARIABLES_MAP_H

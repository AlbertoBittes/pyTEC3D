#ifndef VARIABLES_MAP_H
#define VARIABLES_MAP_H

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>

#include "apec_m_system.h"

namespace libMesh
{
    class System;
}

class ElementType;
enum PhysicalRegionEnum : unsigned int;

enum class GlobalVariableName
{
    A_x,
    A_y,
    A_z,
    Az,
    VOLT,
    V,
    EMF,
    Js_x,
    Js_y,
    Js_z,
    AMPS,
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
    V,
    EMF,
    AMPS,
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
    ApecMSystemBase::SystemType system_type;
};

class VariablesMap
{
public:
    static const std::map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP;

    struct VectorComponent
    {
        unsigned int var_number = libMesh::invalid_uint;
        unsigned int component_index = libMesh::invalid_uint;
        std::string api_label;
        std::string global_label;
    };

    struct PhysicalComponentsInfo
    {
        unsigned int physical_components;
        bool is_scalar;
    };

    static const VariableInfo &get_info_by_key(GlobalVariableName key);
    static std::vector<GlobalVariableName> get_keys_by_predicate(
        const std::function<bool(const VariableInfo &)> &pred);
    static GlobalVariableName get_key_by_predicate_one(
        const std::function<bool(const VariableInfo &)> &pred);
    static void sort_vector_components(
        std::map<std::string, std::vector<VectorComponent>> &vector_variables);
    static std::vector<std::string> get_global_labels_for_label(const std::string &label);

    static std::vector<GlobalVariableName> get_variables_for_element(const ElementType &element);
    static std::vector<GlobalVariableName> get_active_variables_for_element(const ElementType &element,
                                                                            ApecMSystemBase::SystemType system_type);

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

    static PhysicalComponentsInfo get_variable_components_for_system(const libMesh::System &system,
                                                                     GlobalVariableName variable);
    static std::map<std::string, std::vector<VectorComponent>> get_vector_variables_for_system(
        const libMesh::System &system);

    static std::vector<GlobalVariableName> get_system_variables(const libMesh::System &system);

private:
    static std::vector<GlobalVariableName> add_av_variables();
    static std::vector<GlobalVariableName> add_av3d_variables();
    static std::vector<GlobalVariableName> add_solid97_variables(PhysicalRegionEnum physical_region);
    static std::vector<GlobalVariableName> add_solid236_variables(PhysicalRegionEnum physical_region);
    static std::vector<GlobalVariableName> add_plane233_variables(const std::vector<unsigned int> &keyopts);
};

#endif // VARIABLES_MAP_H

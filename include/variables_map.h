/**
 * @file variables_map.h
 * @brief Physical variable registry and FE/system-type resolution utilities.
 *
 * Defines the compile-time catalogue of all physical quantities used in pyTEC3D
 * (DOF unknowns A, Az, VOLT, EMF; nodal fields B, H; element fields Be, He, Jt)
 * and provides helpers that map between enum keys, string labels, FE families,
 * integration orders, and solution system types.
 */
#ifndef VARIABLES_MAP_H
#define VARIABLES_MAP_H

#include <algorithm>
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <libmesh/enum_fe_family.h>
#include <libmesh/enum_order.h>

#include <set>

#include <libmesh/id_types.h>

#include "pytec3d_system.h"

namespace libMesh
{
    class System;
    class EquationSystems;
}

namespace pytec3d
{
    class ElementType;
}

class ElementTypesFactory;

enum class GlobalVariableName
{
    A_x,
    A_y,
    A_z,
    Az,
    V,
    VOLT,
    EMF,
    Js_x,
    Js_y,
    Js_z,
    AMPS,
    CURR,
    B_x,
    B_y,
    B_z,
    H_x,
    H_y,
    H_z,
    Be_x,
    Be_y,
    Be_z,
    He_x,
    He_y,
    He_z,
    Jt_x,
    Jt_y,
    Jt_z,
    UNKNOWN // Fallback for invalid strings
};

enum class APIVariableName
{
    A,
    Az,
    V,
    VOLT,
    EMF,
    AMPS,
    CURR,
    Js,
    B,
    H,
    Be,
    He,
    Jt,
    UNKNOWN // Fallback for invalid strings
};

using VariableToElementIdsMapType = std::map<GlobalVariableName, std::set<libMesh::subdomain_id_type>>;

struct VariableInfo
{
    GlobalVariableName global_name;
    APIVariableName api_name;
    std::string global_label;
    std::string api_label;
    unsigned int component_index;
};

class VariablesMap
{
public:
    static const std::map<GlobalVariableName, VariableInfo> GLOBAL_VARIABLE_INFO_MAP;
    static const std::map<std::string, std::string> DOF_VARIABLE_NODAL_LOAD_MAP;

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
    static const VariableInfo &resolve_variable_info(const std::string &api_label,
                                                     unsigned int component);
    static unsigned int resolve_variable_number(const libMesh::System &system,
                                                const std::string &api_label,
                                                unsigned int component);
    static void sort_vector_components(
        std::map<std::string, std::vector<VectorComponent>> &vector_variables);
    static std::vector<std::string> get_global_labels_for_label(const std::string &label);

    /** @brief One libMesh scalar (or H(curl)) variable of a system, as resolved from an element type. */
    struct SystemVariable
    {
        GlobalVariableName key;
        libMesh::FEFamily family;
        libMesh::Order order;
        std::string variable; ///< Name of the element type variable it belongs to (A for A_x, A_y, A_z).
    };

    /** @brief libMesh variables the fields of the element type active in @p domain contribute to @p system_type. */
    static std::vector<SystemVariable> get_system_variables_for_element(const pytec3d::ElementType &element,
                                                                        unsigned int domain,
                                                                        pyTEC3DSystemBase::SystemType system_type);

    /** @brief For every variable of @p system_type, the subdomains whose element type contains it. */
    static VariableToElementIdsMapType associate_variables_with_element_IDs(const ElementTypesFactory &factory,
                                                                            pyTEC3DSystemBase::SystemType system_type);

    static PhysicalComponentsInfo get_variable_components_for_system(const libMesh::System &system,
                                                                     GlobalVariableName variable);
    static std::map<std::string, std::vector<VectorComponent>> get_vector_variables_for_system(
        const libMesh::System &system);

    static std::vector<GlobalVariableName> get_system_variables(const libMesh::System &system);
    static std::map<std::string, std::string> get_dof_variable_nodal_load_map_for_system(const libMesh::System &system);
};

#endif // VARIABLES_MAP_H

/**
 * @file variables_map.cpp
 * @brief VariablesMap static data, lookup helpers, and per-element variable set builders.
 *
 * The static GLOBAL_VARIABLE_INFO_MAP table is the authoritative registry of all
 * physical quantities.  Per-formulation helpers (add_av3dnodal_variables, etc.)
 * return the variable sets required by each element type and physical-region option.
 */
#include "variables_map.h"
#include "element_type.h"
#include "element_types_factory.h"
#include "logging.h"
#include <libmesh/equation_systems.h>
#include <libmesh/system.h>

using libMesh::CONSTANT;
using libMesh::FEFamily;
using libMesh::FIRST;
using libMesh::L2_LAGRANGE;
using libMesh::L2_LAGRANGE_VEC;
using libMesh::LAGRANGE;
using libMesh::LAGRANGE_VEC;
using libMesh::MONOMIAL;
using libMesh::MONOMIAL_VEC;
using libMesh::NEDELEC_ONE;
using libMesh::Order;
using libMesh::subdomain_id_type;
using libMesh::System;

namespace
{
    bool has_label_component(const std::string &base_label,
                             const char *suffix,
                             const std::vector<GlobalVariableName> &vars)
    {
        std::string label = base_label + suffix;
        for (const auto &var : vars)
        {
            if (VariablesMap::get_info_by_key(var).global_label == label)
                return true;
        }
        return false;
    }

    unsigned int count_label_components(const std::string &base_label,
                                        const std::vector<GlobalVariableName> &vars)
    {
        unsigned int count = 0;
        if (has_label_component(base_label, "_x", vars))
            ++count;
        if (has_label_component(base_label, "_y", vars))
            ++count;
        if (has_label_component(base_label, "_z", vars))
            ++count;
        return count > 0 ? count : 1;
    }

    bool has_vector_suffix(const std::string &label)
    {
        return label.size() > 2 && label[label.size() - 2] == '_' &&
               (label[label.size() - 1] == 'x' || label[label.size() - 1] == 'y' || label[label.size() - 1] == 'z');
    }
}

const std::map<GlobalVariableName, VariableInfo> VariablesMap::GLOBAL_VARIABLE_INFO_MAP = {
    // DOF Variables
    {GlobalVariableName::A_x, {GlobalVariableName::A_x, APIVariableName::A, "A_x", "A", 0}},
    {GlobalVariableName::A_y, {GlobalVariableName::A_y, APIVariableName::A, "A_y", "A", 1}},
    {GlobalVariableName::A_z, {GlobalVariableName::A_z, APIVariableName::A, "A_z", "A", 2}},
    {GlobalVariableName::Az, {GlobalVariableName::Az, APIVariableName::Az, "Az", "Az", 0}},
    {GlobalVariableName::V, {GlobalVariableName::V, APIVariableName::V, "V", "V", 0}},
    {GlobalVariableName::VOLT, {GlobalVariableName::VOLT, APIVariableName::VOLT, "VOLT", "VOLT", 0}},
    {GlobalVariableName::EMF, {GlobalVariableName::EMF, APIVariableName::EMF, "EMF", "EMF", 0}},

    // Nodal Solution Variables
    {GlobalVariableName::B_x, {GlobalVariableName::B_x, APIVariableName::B, "B_x", "B", 0}},
    {GlobalVariableName::B_y, {GlobalVariableName::B_y, APIVariableName::B, "B_y", "B", 1}},
    {GlobalVariableName::B_z, {GlobalVariableName::B_z, APIVariableName::B, "B_z", "B", 2}},
    {GlobalVariableName::H_x, {GlobalVariableName::H_x, APIVariableName::H, "H_x", "H", 0}},
    {GlobalVariableName::H_y, {GlobalVariableName::H_y, APIVariableName::H, "H_y", "H", 1}},
    {GlobalVariableName::H_z, {GlobalVariableName::H_z, APIVariableName::H, "H_z", "H", 2}},

    // Element Solution Variables
    {GlobalVariableName::Be_x, {GlobalVariableName::Be_x, APIVariableName::Be, "Be_x", "Be", 0}},
    {GlobalVariableName::Be_y, {GlobalVariableName::Be_y, APIVariableName::Be, "Be_y", "Be", 1}},
    {GlobalVariableName::Be_z, {GlobalVariableName::Be_z, APIVariableName::Be, "Be_z", "Be", 2}},
    {GlobalVariableName::He_x, {GlobalVariableName::He_x, APIVariableName::He, "He_x", "He", 0}},
    {GlobalVariableName::He_y, {GlobalVariableName::He_y, APIVariableName::He, "He_y", "He", 1}},
    {GlobalVariableName::He_z, {GlobalVariableName::He_z, APIVariableName::He, "He_z", "He", 2}},
    {GlobalVariableName::Jt_x, {GlobalVariableName::Jt_x, APIVariableName::Jt, "Jt_x", "Jt", 0}},
    {GlobalVariableName::Jt_y, {GlobalVariableName::Jt_y, APIVariableName::Jt, "Jt_y", "Jt", 1}},
    {GlobalVariableName::Jt_z, {GlobalVariableName::Jt_z, APIVariableName::Jt, "Jt_z", "Jt", 2}},
};

const std::map<std::string, std::string> VariablesMap::DOF_VARIABLE_NODAL_LOAD_MAP = {
    {"V", "AMPS"},
    {"VOLT", "CURR"},
};

const VariableInfo &VariablesMap::get_info_by_key(GlobalVariableName key)
{
    auto it = GLOBAL_VARIABLE_INFO_MAP.find(key);
    if (it == GLOBAL_VARIABLE_INFO_MAP.end())
        pytec3d::logging::error_and_abort("Key not found in info map");
    return it->second;
}

std::vector<GlobalVariableName> VariablesMap::get_keys_by_predicate(
    const std::function<bool(const VariableInfo &)> &pred)
{
    std::vector<GlobalVariableName> result;
    for (const auto &pair : GLOBAL_VARIABLE_INFO_MAP)
        if (pred(pair.second))
            result.push_back(pair.first);
    return result;
}

GlobalVariableName VariablesMap::get_key_by_predicate_one(
    const std::function<bool(const VariableInfo &)> &pred)
{
    auto keys = get_keys_by_predicate(pred);
    if (keys.size() == 1)
        return keys[0];
    pytec3d::logging::error_and_abort("get_key_by_predicate_one: mapping is not one-to-one or key not found");
}

const VariableInfo &VariablesMap::resolve_variable_info(const std::string &api_label,
                                                        unsigned int component)
{
    const auto key = VariablesMap::get_key_by_predicate_one(
        [&api_label, component](const VariableInfo &info)
        {
            return info.api_label == api_label && info.component_index == component;
        });
    return VariablesMap::get_info_by_key(key);
}

unsigned int VariablesMap::resolve_variable_number(const libMesh::System &system,
                                                   const std::string &api_label,
                                                   unsigned int component)
{
    const auto &info = VariablesMap::resolve_variable_info(api_label, component);
    if (!system.has_variable(info.global_label))
        pytec3d::logging::error_and_abort("System does not contain variable: " + info.global_label);
    return system.variable_number(info.global_label);
}

void VariablesMap::sort_vector_components(
    std::map<std::string, std::vector<VectorComponent>> &vector_variables)
{
    for (auto &entry : vector_variables)
    {
        auto &components = entry.second;
        std::sort(components.begin(), components.end(),
                  [](const auto &a, const auto &b)
                  { return a.component_index < b.component_index; });
    }
}

std::vector<std::string> VariablesMap::get_global_labels_for_label(const std::string &label)
{
    const auto matches_global = VariablesMap::get_keys_by_predicate(
        [&label](const VariableInfo &info)
        {
            return info.global_label == label;
        });

    if (!matches_global.empty())
        return {label};

    auto matches_api = VariablesMap::get_keys_by_predicate(
        [&label](const VariableInfo &info)
        {
            return info.api_label == label;
        });

    if (matches_api.empty())
        pytec3d::logging::error_and_abort("Dof couplings: unknown variable label '" + label + "'.");

    std::sort(matches_api.begin(), matches_api.end(),
              [](GlobalVariableName lhs, GlobalVariableName rhs)
              {
                  return VariablesMap::get_info_by_key(lhs).component_index <
                         VariablesMap::get_info_by_key(rhs).component_index;
              });

    std::vector<std::string> labels;
    labels.reserve(matches_api.size());
    for (const auto &key : matches_api)
        labels.push_back(VariablesMap::get_info_by_key(key).global_label);
    return labels;
}

namespace
{
    pytec3d::FieldKind field_kind_for(pyTEC3DSystemBase::SystemType system_type)
    {
        switch (system_type)
        {
        case pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM:
            return pytec3d::FieldKind::Dof;
        case pyTEC3DSystemBase::NODAL_SOLUTION:
            return pytec3d::FieldKind::Nodal;
        case pyTEC3DSystemBase::ELEMENT_SOLUTION:
            return pytec3d::FieldKind::Element;
        default:
            pytec3d::logging::error_and_abort("VariablesMap: unknown system type.");
        }
    }
} // namespace

std::vector<VariablesMap::SystemVariable> VariablesMap::get_system_variables_for_element(
    const pytec3d::ElementType &element,
    pytec3d::registry::Domain domain,
    pyTEC3DSystemBase::SystemType system_type)
{
    std::vector<SystemVariable> out;
    for (const pytec3d::Field &field : element.fields(field_kind_for(system_type)))
    {
        if (!field.active_in(domain))
            continue;
        for (const pytec3d::LibMeshVariable &variable : field.libmesh_variables)
        {
            const auto key = get_key_by_predicate_one([&variable](const VariableInfo &info)
                                                      { return info.global_label == variable.label; });
            out.push_back({key, variable.family, variable.order, std::string(field.name())});
        }
    }
    return out;
}

VariableToElementIdsMapType VariablesMap::associate_variables_with_element_IDs(
    const ElementTypesFactory &factory,
    pyTEC3DSystemBase::SystemType system_type)
{
    VariableToElementIdsMapType variable_to_element_ids;
    for (const auto &pair : factory.entries())
        for (const SystemVariable &variable : get_system_variables_for_element(*pair.second.type, pair.second.domain, system_type))
            variable_to_element_ids[variable.key].insert(pair.first);
    return variable_to_element_ids;
}

std::vector<GlobalVariableName> VariablesMap::get_system_variables(const System &system)
{
    std::vector<GlobalVariableName> vars;
    const unsigned int n_vars = system.n_vars();
    vars.reserve(n_vars);
    for (unsigned int i = 0; i < n_vars; ++i)
    {
        const std::string &name = system.variable_name(i);
        auto keys = VariablesMap::get_keys_by_predicate(
            [&name](const VariableInfo &info)
            {
                return info.global_label == name;
            });
        if (keys.empty())
            continue;
        if (keys.size() > 1)
            pytec3d::logging::error_and_abort("Variable name maps to multiple global variables: " + name);
        vars.push_back(keys.front());
    }
    return vars;
}

std::map<std::string, std::string> VariablesMap::get_dof_variable_nodal_load_map_for_system(const System &system)
{
    std::map<std::string, std::string> system_map;
    const unsigned int n_vars = system.n_vars();
    for (unsigned int i = 0; i < n_vars; ++i)
    {
        const std::string &var_name = system.variable_name(i);
        const auto it = DOF_VARIABLE_NODAL_LOAD_MAP.find(var_name);
        if (it == DOF_VARIABLE_NODAL_LOAD_MAP.end())
            continue;
        system_map.emplace(it->first, it->second);
    }
    return system_map;
}

VariablesMap::PhysicalComponentsInfo VariablesMap::get_variable_components_for_system(const System &system,
                                                                                      GlobalVariableName variable)
{
    const std::string &label = get_info_by_key(variable).global_label;
    if (!has_vector_suffix(label))
        return {1, true};

    std::string base_label = label.substr(0, label.size() - 2);
    const auto vars = VariablesMap::get_system_variables(system);

    unsigned int count = count_label_components(base_label, vars);
    return {count, count == 1};
}

std::map<std::string, std::vector<VariablesMap::VectorComponent>>
VariablesMap::get_vector_variables_for_system(const System &system)
{
    std::map<std::string, std::vector<VectorComponent>> vector_variables;
    const unsigned int n_vars = system.n_vars();
    for (unsigned int var_num = 0; var_num < n_vars; ++var_num)
    {
        const std::string &var_name = system.variable_name(var_num);
        const GlobalVariableName global_var = get_key_by_predicate_one(
            [&var_name](const VariableInfo &info)
            {
                return info.global_label == var_name;
            });
        const auto &var_info = get_info_by_key(global_var);
        auto components_info = get_variable_components_for_system(system, global_var);
        if (components_info.is_scalar)
            continue;

        vector_variables[var_info.api_label].push_back(
            {var_num, var_info.component_index, var_info.api_label, var_info.global_label});
    }

    sort_vector_components(vector_variables);
    return vector_variables;
}

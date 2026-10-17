/**
 * @file variables_map.cpp
 * @brief VariablesMap: lookups on the generated variable registry and the variables of the element types.
 */
#include "variables_map.h"

#include <algorithm>

#include <libmesh/dof_map.h>
#include <libmesh/system.h>

#include "element_type.h"
#include "model_element_types.h"
#include "logging.h"

namespace
{
    using pytec3d::registry::VariableKind;

    const std::vector<VariableInfo> &variable_infos()
    {
        static const std::vector<VariableInfo> infos = []
        {
            std::vector<VariableInfo> out;
            for (std::size_t i = 0; i < pytec3d::registry::kNumVariableId; ++i)
            {
                const auto id = static_cast<GlobalVariableName>(i);
                const auto &info = pytec3d::registry::info(id);
                out.push_back({id, std::string(info.label), std::string(info.api_label),
                               info.kind == VariableKind::Component ? info.component : 0u, info.kind == VariableKind::Vector});
            }
            return out;
        }();
        return infos;
    }

    const VariableInfo *find_label(const std::string &label)
    {
        for (const VariableInfo &info : variable_infos())
            if (info.global_label == label)
                return &info;
        return nullptr;
    }

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

const VariableInfo &VariablesMap::get_info_by_key(GlobalVariableName key)
{
    return variable_infos().at(static_cast<std::size_t>(key));
}

std::vector<GlobalVariableName> VariablesMap::get_keys_by_predicate(
    const std::function<bool(const VariableInfo &)> &pred)
{
    std::vector<GlobalVariableName> result;
    for (const VariableInfo &info : variable_infos())
        if (pred(info))
            result.push_back(info.global_name);
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

bool VariablesMap::is_edge_variable(const libMesh::System &system, unsigned int var)
{
    return system.variable_type(var).family == libMesh::NEDELEC_ONE;
}

bool VariablesMap::has_edge_variable(const libMesh::System &system)
{
    for (unsigned int var = 0; var < system.n_vars(); ++var)
        if (is_edge_variable(system, var))
            return true;
    return false;
}

std::string VariablesMap::constraint_array_name(const libMesh::System &system, unsigned int var)
{
    const VariableInfo *info = find_label(system.variable_name(var));
    const std::string api_label = info ? info->api_label : system.variable_name(var);
    return "_" + api_label;
}

std::vector<std::string> VariablesMap::get_global_labels_for_label(const libMesh::System &system, const std::string &label)
{
    std::vector<std::string> labels;
    if (system.has_variable(label))
        return {label};

    // The variables of a field: by api label, x, y, z.
    std::vector<const VariableInfo *> matches;
    for (unsigned int var = 0; var < system.n_vars(); ++var)
    {
        const VariableInfo *info = find_label(system.variable_name(var));
        if (info && info->api_label == label)
            matches.push_back(info);
    }
    if (matches.empty())
        pytec3d::logging::error_and_abort("Dof couplings: '" + label + "' is not a variable of system '" + system.name() + "'.");
    std::sort(matches.begin(), matches.end(), [](const VariableInfo *a, const VariableInfo *b)
              { return a->component_index < b->component_index; });
    for (const VariableInfo *info : matches)
        labels.push_back(info->global_label);
    return labels;
}

std::vector<VariablesMap::SystemVariable> VariablesMap::get_system_variables_for_element(
    const pytec3d::ElementType &element,
    pytec3d::registry::Domain domain,
    pyTEC3DSystemBase::SystemType system_type)
{
    const pytec3d::FieldKind kind = field_kind_for(system_type);
    std::vector<SystemVariable> out;
    for (const pytec3d::FieldSpec &field : element.spec().fields)
    {
        if (field.kind != kind || !pytec3d::ElementType::active(field, domain))
            continue;
        const pytec3d::LibMeshVariable variable = element.variable(field);
        out.push_back({variable.id, variable.family, variable.order});
    }
    return out;
}

VariableToElementIdsMapType VariablesMap::associate_variables_with_element_IDs(
    const ModelElementTypes &factory,
    pyTEC3DSystemBase::SystemType system_type)
{
    VariableToElementIdsMapType variable_to_element_ids;
    for (const auto &pair : factory.entries())
        for (const SystemVariable &variable : get_system_variables_for_element(*pair.second.type, pair.second.domain, system_type))
            variable_to_element_ids[variable.key].insert(pair.first);
    return variable_to_element_ids;
}

std::map<std::string, std::string> VariablesMap::get_dof_variable_nodal_load_map_for_system(const libMesh::System &system)
{
    std::map<std::string, std::string> system_map;
    for (unsigned int var = 0; var < system.n_vars(); ++var)
    {
        const VariableInfo *info = find_label(system.variable_name(var));
        if (!info)
            continue;
        const std::string_view load = pytec3d::registry::info(info->global_name).nodal_load;
        if (!load.empty())
            system_map.emplace(info->global_label, std::string(load));
    }
    return system_map;
}

std::map<std::string, std::vector<VariablesMap::VectorComponent>>
VariablesMap::get_vector_variables_for_system(const libMesh::System &system)
{
    std::map<std::string, std::vector<VectorComponent>> vector_variables;
    for (unsigned int var_num = 0; var_num < system.n_vars(); ++var_num)
    {
        const VariableInfo *info = find_label(system.variable_name(var_num));
        if (!info || is_edge_variable(system, var_num))
            continue;
        if (info->vector)
        {
            for (unsigned int c = 0; c < 3; ++c)
                vector_variables[info->api_label].push_back({var_num, c, c, info->api_label, info->global_label});
        }
        else if (pytec3d::registry::info(info->global_name).kind == VariableKind::Component)
            vector_variables[info->api_label].push_back({var_num, info->component_index, 0, info->api_label, info->global_label});
    }

    // A single component is a scalar.
    for (auto it = vector_variables.begin(); it != vector_variables.end();)
    {
        if (it->second.size() < 2)
            it = vector_variables.erase(it);
        else
        {
            std::sort(it->second.begin(), it->second.end(), [](const VectorComponent &a, const VectorComponent &b)
                      { return a.component_index < b.component_index; });
            ++it;
        }
    }
    return vector_variables;
}

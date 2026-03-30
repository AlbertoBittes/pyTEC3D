#include "variables_map.h"
#include "element_types.h"
#include <libmesh/equation_systems.h>
#include <libmesh/system.h>

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

const VariableInfo &VariablesMap::get_info_by_key(GlobalVariableName key)
{
    auto it = GLOBAL_VARIABLE_INFO_MAP.find(key);
    if (it == GLOBAL_VARIABLE_INFO_MAP.end())
        libmesh_error_msg("Key not found in info map");
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
    libmesh_error_msg("get_key_by_predicate_one: mapping is not one-to-one or key not found");
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
        libmesh_error_msg("System does not contain variable: " + info.global_label);
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
        libmesh_error_msg("Dof couplings: unknown variable label '" + label + "'.");

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

VariableToElementIdsMapType VariablesMap::associate_variables_with_element_IDs(
    const ElementTypesMap &element_types_map,
    ApecMSystemBase::SystemType system_type)
{
    VariableToElementIdsMapType variable_to_element_ids;
    for (const auto &pair : element_types_map)
    {
        subdomain_id_type element_id = pair.first;
        const auto &element_type = pair.second;
        if (!element_type)
            continue;

        const auto vars = VariablesMap::get_variables_for_element(*element_type);
        for (const auto &var : vars)
        {
            if (VariablesMap::get_variable_system_type(*element_type, var) == system_type)
                variable_to_element_ids[var].insert(element_id);
        }
    }
    return variable_to_element_ids;
}

std::vector<GlobalVariableName> VariablesMap::get_variables_for_element(const ElementType &element)
{
    std::vector<GlobalVariableName> vars;
    const auto &OPTIONs = element.get_OPTIONs();
    auto type = static_cast<ElementTypeEnum>(OPTIONs.at(0));

    switch (type)
    {
    case ElementTypeEnum::AV3DNODAL:
        vars = add_av3d_variables();
        {
            auto extra = add_av3dnodal_variables(OPTIONs);
            vars.insert(vars.end(), extra.begin(), extra.end());
        }
        break;
    case ElementTypeEnum::AV3DEDGE:
        vars = add_av3d_variables();
        {
            auto extra = add_av3dedge_variables(OPTIONs);
            vars.insert(vars.end(), extra.begin(), extra.end());
        }
        break;
    case ElementTypeEnum::AV2DNODAL:
        vars = add_av_variables();
        {
            auto extra = add_av2dnodal_variables(OPTIONs);
            vars.insert(vars.end(), extra.begin(), extra.end());
        }
        break;
    default:
        libmesh_error_msg("VariablesMap: unsupported element type.");
    }
    return vars;
}

std::vector<GlobalVariableName> VariablesMap::get_active_variables_for_element(const ElementType &element,
                                                                               ApecMSystemBase::SystemType system_type)
{
    std::vector<GlobalVariableName> active_vars;
    auto vars = get_variables_for_element(element);
    for (const auto &var : vars)
    {
        if (get_variable_system_type(element, var) == system_type)
            active_vars.push_back(var);
    }
    std::sort(active_vars.begin(), active_vars.end(),
              [](GlobalVariableName lhs, GlobalVariableName rhs)
              {
                  return static_cast<unsigned int>(lhs) < static_cast<unsigned int>(rhs);
              });
    return active_vars;
}

std::vector<APIVariableName> VariablesMap::get_api_variables_for_element(const ElementType &element,
                                                                         ApecMSystemBase::SystemType system_type)
{
    auto vars = get_active_variables_for_element(element, system_type);
    std::vector<APIVariableName> api_vars;
    api_vars.reserve(vars.size());
    for (const auto &var : vars)
    {
        auto api_name = get_info_by_key(var).api_name;
        if (std::find(api_vars.begin(), api_vars.end(), api_name) == api_vars.end())
            api_vars.push_back(api_name);
    }
    return api_vars;
}

ApecMSystemBase::SystemType VariablesMap::get_variable_system_type(const ElementType &,
                                                                   GlobalVariableName variable)
{
    switch (variable)
    {
    case GlobalVariableName::A_x:
    case GlobalVariableName::A_y:
    case GlobalVariableName::A_z:
    case GlobalVariableName::Az:
    case GlobalVariableName::VOLT:
    case GlobalVariableName::EMF:
        return ApecMSystemBase::DOF_SOLUTION_SYSTEM;
    case GlobalVariableName::B_x:
    case GlobalVariableName::B_y:
    case GlobalVariableName::B_z:
    case GlobalVariableName::H_x:
    case GlobalVariableName::H_y:
    case GlobalVariableName::H_z:
        return ApecMSystemBase::NODAL_SOLUTION;
    case GlobalVariableName::Be_x:
    case GlobalVariableName::Be_y:
    case GlobalVariableName::Be_z:
    case GlobalVariableName::He_x:
    case GlobalVariableName::He_y:
    case GlobalVariableName::He_z:
    case GlobalVariableName::Jt_x:
    case GlobalVariableName::Jt_y:
    case GlobalVariableName::Jt_z:
        return ApecMSystemBase::ELEMENT_SOLUTION;
    default:
        libmesh_error_msg("VariablesMap: missing system type for variable.");
    }
}

std::vector<GlobalVariableName> VariablesMap::add_av_variables()
{
    return {GlobalVariableName::B_x,
            GlobalVariableName::B_y,
            GlobalVariableName::H_x,
            GlobalVariableName::H_y,
            GlobalVariableName::Jt_z,
            GlobalVariableName::Be_x,
            GlobalVariableName::Be_y,
            GlobalVariableName::He_x,
            GlobalVariableName::He_y};
}

std::vector<GlobalVariableName> VariablesMap::add_av3d_variables()
{
    auto vars = add_av_variables();
    vars.insert(vars.end(), {GlobalVariableName::B_z,
                             GlobalVariableName::H_z,
                             GlobalVariableName::Jt_x,
                             GlobalVariableName::Jt_y,
                             GlobalVariableName::Be_z,
                             GlobalVariableName::He_z});
    return vars;
}

std::vector<GlobalVariableName> VariablesMap::add_av3dnodal_variables(const std::vector<unsigned int> &OPTIONs)
{
    const PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(OPTIONs[1]);
    std::vector<GlobalVariableName> vars = {GlobalVariableName::A_x,
                                            GlobalVariableName::A_y,
                                            GlobalVariableName::A_z};

    switch (physical_region)
    {
    case SHORT_CIRCUIT_CONDUCTOR:
    case SOLID_CONDUCTOR:
    case STRANDED_COIL:
        if (OPTIONs[4] == 0 || OPTIONs[4] == 1)
            vars.insert(vars.end(), {GlobalVariableName::VOLT});
        if (physical_region == STRANDED_COIL)
            vars.insert(vars.end(), {GlobalVariableName::EMF});
        break;
    default:
        break;
    }
    return vars;
}

std::vector<GlobalVariableName> VariablesMap::add_av3dedge_variables(const std::vector<unsigned int> &OPTIONs)
{
    const PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(OPTIONs[1]);
    std::vector<GlobalVariableName> vars = {GlobalVariableName::Az};

    switch (physical_region)
    {
    case SHORT_CIRCUIT_CONDUCTOR:
    case SOLID_CONDUCTOR:
    case STRANDED_COIL:
        if (OPTIONs[4] == 0 || OPTIONs[4] == 1)
            vars.insert(vars.end(), {GlobalVariableName::VOLT});
        if (physical_region == STRANDED_COIL)
            vars.insert(vars.end(), {GlobalVariableName::EMF});
        break;
    default:
        break;
    }
    return vars;
}

std::vector<GlobalVariableName> VariablesMap::add_av2dnodal_variables(const std::vector<unsigned int> &OPTIONs)
{
    std::vector<GlobalVariableName> vars = {GlobalVariableName::A_z};
    const PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(OPTIONs[1]);

    switch (physical_region)
    {
    case SOLID_CONDUCTOR:
    case STRANDED_COIL:
        if (OPTIONs[4] == 0 || OPTIONs[4] == 1)
            vars.insert(vars.end(), {GlobalVariableName::VOLT});
        if (physical_region == STRANDED_COIL)
            vars.insert(vars.end(), {GlobalVariableName::EMF});
        break;
    default:
        break;
    }
    return vars;
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
            libmesh_error_msg("Variable name maps to multiple global variables: " + name);
        vars.push_back(keys.front());
    }
    return vars;
}

FEFamily VariablesMap::get_variable_fe_family(const ElementType &,
                                              GlobalVariableName variable)
{
    switch (variable)
    {
    case GlobalVariableName::A_x:
    case GlobalVariableName::A_y:
    case GlobalVariableName::A_z:
    case GlobalVariableName::VOLT:
    case GlobalVariableName::EMF:
        return LAGRANGE;
    case GlobalVariableName::Az:
        return NEDELEC_ONE;
    case GlobalVariableName::Jt_x:
    case GlobalVariableName::Jt_y:
    case GlobalVariableName::Jt_z:
    case GlobalVariableName::Be_x:
    case GlobalVariableName::Be_y:
    case GlobalVariableName::Be_z:
    case GlobalVariableName::He_x:
    case GlobalVariableName::He_y:
    case GlobalVariableName::He_z:
        return MONOMIAL;
    case GlobalVariableName::B_x:
    case GlobalVariableName::B_y:
    case GlobalVariableName::B_z:
    case GlobalVariableName::H_x:
    case GlobalVariableName::H_y:
    case GlobalVariableName::H_z:
        return L2_LAGRANGE;
    default:
        libmesh_error_msg("VariablesMap: missing FE family for variable.");
    }
}

Order VariablesMap::get_variable_order(const ElementType &element,
                                       GlobalVariableName variable)
{
    Order order = FIRST;
    switch (variable)
    {
    case GlobalVariableName::Jt_x:
    case GlobalVariableName::Jt_y:
    case GlobalVariableName::Jt_z:
    case GlobalVariableName::Be_x:
    case GlobalVariableName::Be_y:
    case GlobalVariableName::Be_z:
    case GlobalVariableName::He_x:
    case GlobalVariableName::He_y:
    case GlobalVariableName::He_z:
        order = CONSTANT;
        break;
    default:
        order = FIRST;
        break;
    }

    const auto &OPTIONs = element.get_OPTIONs();
    FEFamily family = get_variable_fe_family(element, variable);
    if (OPTIONs.size() > 2)
    {
        auto OPTION_order = static_cast<Order>(OPTIONs[2]);
        if (OPTION_order != FIRST &&
            (family == LAGRANGE || family == L2_LAGRANGE || family == LAGRANGE_VEC || family == L2_LAGRANGE_VEC))
        {
            order = OPTION_order;
        }
    }
    return order;
}

unsigned int VariablesMap::get_variable_physical_components(const ElementType &element,
                                                            GlobalVariableName variable)
{
    const std::string &label = get_info_by_key(variable).global_label;
    if (!has_vector_suffix(label))
        return 1;

    std::string base_label = label.substr(0, label.size() - 2);
    const auto vars = get_variables_for_element(element);

    return count_label_components(base_label, vars);
}

bool VariablesMap::is_variable_scalar(const ElementType &element,
                                      GlobalVariableName variable)
{
    return get_variable_physical_components(element, variable) == 1;
}

bool VariablesMap::is_variable_vector(const ElementType &element,
                                      GlobalVariableName variable)
{
    return !is_variable_scalar(element, variable);
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

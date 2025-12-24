#include "variables_map.h"
#include "element_types.h"
const std::map<GlobalVariableName, VariableInfo> VariablesMap::GLOBAL_VARIABLE_INFO_MAP = {
        // DOF Variables
        {GlobalVariableName::A_x, {GlobalVariableName::A_x, APIVariableName::A, "A_x", "A", 0, Globals::DOF_SOLUTION_SYSTEM}},
        {GlobalVariableName::A_y, {GlobalVariableName::A_y, APIVariableName::A, "A_y", "A", 1, Globals::DOF_SOLUTION_SYSTEM}},
        {GlobalVariableName::A_z, {GlobalVariableName::A_z, APIVariableName::A, "A_z", "A", 2, Globals::DOF_SOLUTION_SYSTEM}},
        {GlobalVariableName::Az, {GlobalVariableName::Az, APIVariableName::Az, "Az", "Az", 0, Globals::DOF_SOLUTION_SYSTEM}},
        {GlobalVariableName::VOLT, {GlobalVariableName::VOLT, APIVariableName::VOLT, "VOLT", "VOLT", 0, Globals::DOF_SOLUTION_SYSTEM}},
        {GlobalVariableName::EMF, {GlobalVariableName::EMF, APIVariableName::EMF, "EMF", "EMF", 0, Globals::DOF_SOLUTION_SYSTEM}},
        {GlobalVariableName::V, {GlobalVariableName::V, APIVariableName::V, "V", "V", 0, Globals::DOF_SOLUTION_SYSTEM}},

        // Nodal Loads Variables
        {GlobalVariableName::CSG_x, {GlobalVariableName::CSG_x, APIVariableName::CSG, "CSG_x", "CSG", 0, Globals::NODAL_LOADS}},
        {GlobalVariableName::CSG_y, {GlobalVariableName::CSG_y, APIVariableName::CSG, "CSG_y", "CSG", 1, Globals::NODAL_LOADS}},
        {GlobalVariableName::CSG_z, {GlobalVariableName::CSG_z, APIVariableName::CSG, "CSG_z", "CSG", 2, Globals::NODAL_LOADS}},
        {GlobalVariableName::AMPS, {GlobalVariableName::AMPS, APIVariableName::AMPS, "AMPS", "AMPS", 0, Globals::NODAL_LOADS}},
        {GlobalVariableName::CURT, {GlobalVariableName::CURT, APIVariableName::CURT, "CURT", "CURT", 0, Globals::NODAL_LOADS}},

        // Element Loads Variables
        {GlobalVariableName::Js_x, {GlobalVariableName::Js_x, APIVariableName::Js, "Js_x", "Js", 0, Globals::ELEMENT_LOADS}},
        {GlobalVariableName::Js_y, {GlobalVariableName::Js_y, APIVariableName::Js, "Js_y", "Js", 1, Globals::ELEMENT_LOADS}},
        {GlobalVariableName::Js_z, {GlobalVariableName::Js_z, APIVariableName::Js, "Js_z", "Js", 2, Globals::ELEMENT_LOADS}},

        // Nodal Solution Variables
        {GlobalVariableName::B_x, {GlobalVariableName::B_x, APIVariableName::B, "B_x", "B", 0, Globals::NODAL_SOLUTION}},
        {GlobalVariableName::B_y, {GlobalVariableName::B_y, APIVariableName::B, "B_y", "B", 1, Globals::NODAL_SOLUTION}},
        {GlobalVariableName::B_z, {GlobalVariableName::B_z, APIVariableName::B, "B_z", "B", 2, Globals::NODAL_SOLUTION}},
        {GlobalVariableName::H_x, {GlobalVariableName::H_x, APIVariableName::H, "H_x", "H", 0, Globals::NODAL_SOLUTION}},
        {GlobalVariableName::H_y, {GlobalVariableName::H_y, APIVariableName::H, "H_y", "H", 1, Globals::NODAL_SOLUTION}},
        {GlobalVariableName::H_z, {GlobalVariableName::H_z, APIVariableName::H, "H_z", "H", 2, Globals::NODAL_SOLUTION}},
        {GlobalVariableName::PSI, {GlobalVariableName::PSI, APIVariableName::PSI, "PSI", "PSI", 0, Globals::NODAL_SOLUTION}},

        // Element Solution Variables
        {GlobalVariableName::Jt_x, {GlobalVariableName::Jt_x, APIVariableName::Jt, "Jt_x", "Jt", 0, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::Jt_y, {GlobalVariableName::Jt_y, APIVariableName::Jt, "Jt_y", "Jt", 1, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::Jt_z, {GlobalVariableName::Jt_z, APIVariableName::Jt, "Jt_z", "Jt", 2, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::Be_x, {GlobalVariableName::Be_x, APIVariableName::Be, "Be_x", "Be", 0, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::Be_y, {GlobalVariableName::Be_y, APIVariableName::Be, "Be_y", "Be", 1, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::Be_z, {GlobalVariableName::Be_z, APIVariableName::Be, "Be_z", "Be", 2, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::He_x, {GlobalVariableName::He_x, APIVariableName::He, "He_x", "He", 0, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::He_y, {GlobalVariableName::He_y, APIVariableName::He, "He_y", "He", 1, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::He_z, {GlobalVariableName::He_z, APIVariableName::He, "He_z", "He", 2, Globals::ELEMENT_SOLUTION}},
        {GlobalVariableName::PSIe, {GlobalVariableName::PSIe, APIVariableName::PSIe, "PSIe", "PSIe", 0, Globals::ELEMENT_SOLUTION}},
    };

    const VariableInfo &VariablesMap::get_info_by_key(GlobalVariableName key)
    {
        return VariablesMap::get_info_by_key(GLOBAL_VARIABLE_INFO_MAP, key);
    }

    GlobalVariableName VariablesMap::get_key_by_global_label(const std::string &global_label)
    {
        return VariablesMap::get_key_by_predicate_one(
            GLOBAL_VARIABLE_INFO_MAP,
            [&global_label](const VariableInfo &info)
            {
                return info.global_label == global_label;
            });
    }

    GlobalVariableName VariablesMap::get_key_by_api_component(APIVariableName api_name, unsigned int component_index)
    {
        return VariablesMap::get_key_by_predicate_one(
            GLOBAL_VARIABLE_INFO_MAP,
            [api_name, component_index](const VariableInfo &info)
            {
                return info.api_name == api_name && info.component_index == component_index;
            });
    }

    std::vector<GlobalVariableName> VariablesMap::get_keys_by_api_label(const std::string &api_label)
    {
        return VariablesMap::get_keys_by_predicate(
            GLOBAL_VARIABLE_INFO_MAP,
            [&api_label](const VariableInfo &info)
            {
                return info.api_label == api_label;
            });
    }

    std::vector<GlobalVariableName> VariablesMap::get_variables_for_element(const ElementType &element)
    {
        std::vector<GlobalVariableName> vars;
        const auto &keyopts = element.get_keyopts();
        PhysicalRegionEnum region = UNKNOWN;
        if (keyopts.size() > 1)
            region = static_cast<PhysicalRegionEnum>(keyopts[1]);
        auto type = static_cast<ElementTypeEnum>(keyopts.at(0));

        switch (type)
        {
        case ElementTypeEnum::SOLID97:
        case ElementTypeEnum::INFIN111:
            vars = add_av3d_variables();
            {
                auto extra = add_solid97_variables(region);
                vars.insert(vars.end(), extra.begin(), extra.end());
            }
            break;
        case ElementTypeEnum::SOLID236:
        case ElementTypeEnum::INFIN112:
            vars = add_av3d_variables();
            {
                auto extra = add_solid236_variables(region);
                vars.insert(vars.end(), extra.begin(), extra.end());
            }
            break;
        case ElementTypeEnum::PLANE233:
        case ElementTypeEnum::INFIN110:
            vars = add_av_variables();
            {
                auto extra = add_plane233_variables(keyopts);
                vars.insert(vars.end(), extra.begin(), extra.end());
            }
            break;
        default:
            libmesh_error_msg("VariablesMap: unsupported element type.");
        }
        return vars;
    }

    std::vector<GlobalVariableName> VariablesMap::get_active_variables_for_element(const ElementType &element,
                                                                                   Globals::SystemType system_type)
    {
        std::vector<GlobalVariableName> active_vars;
        auto vars = get_variables_for_element(element);
        for (const auto &var : vars)
        {
            if (get_info_by_key(var).system_type == system_type)
                active_vars.push_back(var);
        }
        return active_vars;
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

    std::vector<GlobalVariableName> VariablesMap::add_solid97_variables(const PhysicalRegionEnum physical_region)
    {
        std::vector<GlobalVariableName> vars = {GlobalVariableName::A_x,
                                                GlobalVariableName::A_y,
                                                GlobalVariableName::A_z,
                                                GlobalVariableName::CSG_x,
                                                GlobalVariableName::CSG_y,
                                                GlobalVariableName::CSG_z};

        switch (physical_region)
        {
        case STRANDED_CONDUCTOR:
            vars.insert(vars.end(), {GlobalVariableName::Js_x,
                                     GlobalVariableName::Js_y,
                                     GlobalVariableName::Js_z});
            break;
        case SOLID_CONDUCTOR:
        case STRANDED_COIL:
            vars.insert(vars.end(), {GlobalVariableName::VOLT,
                                     GlobalVariableName::AMPS});
            if (physical_region == STRANDED_COIL)
                vars.insert(vars.end(), {GlobalVariableName::EMF,
                                         GlobalVariableName::CURT});
            break;
        case SOLID_CONDUCTOR_CIRCUIT_COUPLED:
        case STRANDED_COIL_CIRCUIT_COUPLED:
            vars.insert(vars.end(), {GlobalVariableName::V,
                                     GlobalVariableName::AMPS});
            if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
                vars.insert(vars.end(), {GlobalVariableName::EMF,
                                         GlobalVariableName::CURT});
            break;
        default:
            break;
        }
        return vars;
    }

    std::vector<GlobalVariableName> VariablesMap::add_solid236_variables(const PhysicalRegionEnum physical_region)
    {
        std::vector<GlobalVariableName> vars = {GlobalVariableName::Az};

        switch (physical_region)
        {
        case STRANDED_CONDUCTOR:
            vars.insert(vars.end(), {GlobalVariableName::Js_x,
                                     GlobalVariableName::Js_y,
                                     GlobalVariableName::Js_z});
            break;
        case SOLID_CONDUCTOR:
        case STRANDED_COIL:
            vars.insert(vars.end(), {GlobalVariableName::VOLT,
                                     GlobalVariableName::AMPS});
            if (physical_region == STRANDED_COIL)
                vars.insert(vars.end(), {GlobalVariableName::EMF,
                                         GlobalVariableName::CURT});
            break;
        case SOLID_CONDUCTOR_CIRCUIT_COUPLED:
        case STRANDED_COIL_CIRCUIT_COUPLED:
            vars.insert(vars.end(), {GlobalVariableName::V,
                                     GlobalVariableName::AMPS});
            if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
                vars.insert(vars.end(), {GlobalVariableName::EMF,
                                         GlobalVariableName::CURT});
            break;
        default:
            break;
        }
        return vars;
    }

    std::vector<GlobalVariableName> VariablesMap::add_plane233_variables(const std::vector<unsigned int> &keyopts)
    {
        std::vector<GlobalVariableName> vars = {GlobalVariableName::A_z,
                                                GlobalVariableName::CSG_z};

        PhysicalRegionEnum physical_region = UNKNOWN;
        if (keyopts.size() > 1)
            physical_region = static_cast<PhysicalRegionEnum>(keyopts[1]);

        switch (physical_region)
        {
        case STRANDED_CONDUCTOR:
            vars.insert(vars.end(), {GlobalVariableName::Js_z});
            break;
        case STRANDED_COIL:
            vars.insert(vars.end(), {GlobalVariableName::VOLT,
                                     GlobalVariableName::AMPS,
                                     GlobalVariableName::EMF,
                                     GlobalVariableName::CURT});
            break;
        case STRANDED_COIL_CIRCUIT_COUPLED:
        case SOLID_CONDUCTOR_CIRCUIT_COUPLED:
            vars.insert(vars.end(), {GlobalVariableName::V,
                                     GlobalVariableName::AMPS});
            if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
                vars.insert(vars.end(), {GlobalVariableName::EMF,
                                         GlobalVariableName::CURT});
            break;
        default:
            break;
        }

        if (keyopts.size() > 3 && keyopts[3] == 1)
        {
            vars.insert(vars.end(), {GlobalVariableName::PSI,
                                     GlobalVariableName::PSIe});
        }
        return vars;
    }

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

        bool has_vector_suffix(const std::string &label)
        {
            return label.size() > 2 && label[label.size() - 2] == '_' &&
                   (label[label.size() - 1] == 'x' || label[label.size() - 1] == 'y' || label[label.size() - 1] == 'z');
        }
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
        case GlobalVariableName::V:
        case GlobalVariableName::CSG_x:
        case GlobalVariableName::CSG_y:
        case GlobalVariableName::CSG_z:
        case GlobalVariableName::AMPS:
        case GlobalVariableName::CURT:
        case GlobalVariableName::PSI:
            return LAGRANGE;
        case GlobalVariableName::Az:
            return NEDELEC_ONE;
        case GlobalVariableName::Js_x:
        case GlobalVariableName::Js_y:
        case GlobalVariableName::Js_z:
        case GlobalVariableName::Jt_x:
        case GlobalVariableName::Jt_y:
        case GlobalVariableName::Jt_z:
        case GlobalVariableName::Be_x:
        case GlobalVariableName::Be_y:
        case GlobalVariableName::Be_z:
        case GlobalVariableName::He_x:
        case GlobalVariableName::He_y:
        case GlobalVariableName::He_z:
        case GlobalVariableName::PSIe:
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
        case GlobalVariableName::Js_x:
        case GlobalVariableName::Js_y:
        case GlobalVariableName::Js_z:
        case GlobalVariableName::Jt_x:
        case GlobalVariableName::Jt_y:
        case GlobalVariableName::Jt_z:
        case GlobalVariableName::Be_x:
        case GlobalVariableName::Be_y:
        case GlobalVariableName::Be_z:
        case GlobalVariableName::He_x:
        case GlobalVariableName::He_y:
        case GlobalVariableName::He_z:
        case GlobalVariableName::PSIe:
            order = CONSTANT;
            break;
        default:
            order = FIRST;
            break;
        }

        const auto &keyopts = element.get_keyopts();
        FEFamily family = get_variable_fe_family(element, variable);
        if (keyopts.size() > 2)
        {
            auto keyopt_order = static_cast<Order>(keyopts[2]);
            if (keyopt_order != FIRST &&
                (family == LAGRANGE || family == L2_LAGRANGE || family == LAGRANGE_VEC))
            {
                order = keyopt_order;
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

        unsigned int count = 0;
        if (has_label_component(base_label, "_x", vars))
            ++count;
        if (has_label_component(base_label, "_y", vars))
            ++count;
        if (has_label_component(base_label, "_z", vars))
            ++count;

        return count > 0 ? count : 1;
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

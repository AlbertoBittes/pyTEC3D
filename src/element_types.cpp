#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <vector>
#include <limits> // for std::numeric_limits
#include "element_types.h"
#include "element_builder.h"
#include <libmesh/inf_fe.h>

using namespace Globals;

ElementType::ElementType(const std::vector<unsigned int> &keyopts,
                         std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                         std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : keyopts(keyopts)
{
    build_variable_information(variable_collector);
    set_variables_order();
    build_fe_array(fe_collector);
}

void ElementType::build_variable_information(std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector)
{
    std::vector<Globals::GlobalVariableName> collected_vars;
    variable_collector(collected_vars);
    for (const auto &var : collected_vars)
        variable_info[var] = Globals::get_info_by_key(Globals::GLOBAL_VARIABLE_INFO_MAP, var);
}

// Virtual methods to be implemented by derived classes
std::vector<Globals::GlobalVariableName> ElementType::get_active_variables(Globals::SystemType variable_system_type) const
{
    std::vector<Globals::GlobalVariableName> active_vars;
    for (const auto &pair : variable_info)
    {
        const Globals::GlobalVariableName &var_name = pair.first;
        const Globals::VariableInfo &var_data = pair.second;

        if (var_data.system_type == variable_system_type)
        {
            active_vars.push_back(var_name);
        }
    }
    return active_vars;
}

void ElementType::build_fe_array(std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
{
    std::vector<FEDataShapeFunctions> fe_structs;
    fe_collector(fe_structs);

    for (const auto &fe_data : fe_structs)
    {
        FEType fe_type(fe_data.order, fe_data.family);

        std::unique_ptr<FEBase> fe_obj;
        std::unique_ptr<FEVectorBase> fe_vec_obj;

        Order quad_order = fe_data.quadrature_order != INVALID_ORDER ? fe_data.quadrature_order : fe_type.default_quadrature_order();

        QGauss *qrule = new QGauss(fe_data.quadrature_dim, quad_order);

        if (fe_data.is_vector)
        {
            if (fe_data.is_infinite)
                fe_vec_obj = FEVectorBase::build_InfFE(fe_data.dimension, fe_type);
            else
                fe_vec_obj = FEVectorBase::build(fe_data.dimension, fe_type);

            fe_vec_obj->attach_quadrature_rule(qrule);
            fe_array.push_back(std::move(fe_vec_obj));
        }
        else
        {
            if (fe_data.is_infinite)
                fe_obj = FEBase::build_InfFE(fe_data.dimension, fe_type);
            else
                fe_obj = FEBase::build(fe_data.dimension, fe_type);

            fe_obj->attach_quadrature_rule(qrule);
            fe_array.push_back(std::move(fe_obj));
        }
    }
}

void ElementType::set_variables_order()
{
    if (static_cast<Order>(keyopts[2]) != FIRST)
        for (auto &pair : variable_info)
        {
            if (pair.second.fe_family == LAGRANGE || pair.second.fe_family == L2_LAGRANGE)
            {
                pair.second.order = static_cast<Order>(keyopts[2]);
            }
        }
}
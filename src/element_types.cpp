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

ElementType::ElementType(const std::vector<unsigned int> &OPTIONs,
                         std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : OPTIONs(OPTIONs)
{
    build_fe_array(fe_collector);
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
            fe_vec_obj = FEVectorBase::build(fe_data.dimension, fe_type);

            fe_vec_obj->attach_quadrature_rule(qrule);
            fe_array.push_back(std::move(fe_vec_obj));
        }
        else
        {
            fe_obj = FEBase::build(fe_data.dimension, fe_type);

            fe_obj->attach_quadrature_rule(qrule);
            fe_array.push_back(std::move(fe_obj));
        }
    }
}

/**
 * @file element_types.cpp
 * @brief ElementType construction: FE array building and safe scalar/vector accessor implementations.
 */
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
#include "logging.h"

using libMesh::FEBase;
using libMesh::FEType;
using libMesh::FEVectorBase;
using libMesh::INVALID_ORDER;
using libMesh::Order;
using libMesh::QBase;
using libMesh::QGauss;

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

    fe_kinds.clear();
    scalar_fe_slots.clear();
    vector_fe_slots.clear();

    for (const auto &fe_data : fe_structs)
    {
        FEType fe_type(fe_data.order, fe_data.family);

        std::unique_ptr<FEBase> fe_obj;
        std::unique_ptr<FEVectorBase> fe_vec_obj;

        Order quad_order = fe_data.quadrature_order != INVALID_ORDER ? fe_data.quadrature_order : fe_type.default_quadrature_order();
        auto qrule = std::make_unique<QGauss>(fe_data.quadrature_dim, quad_order);
        QBase *qrule_ptr = qrule.get();
        quadrature_rules.push_back(std::move(qrule));

        if (fe_data.is_vector)
        {
            fe_vec_obj = FEVectorBase::build(fe_data.dimension, fe_type);

            fe_vec_obj->attach_quadrature_rule(qrule_ptr);
            FEVectorBase *raw_fe_vec = fe_vec_obj.get();
            fe_array.push_back(std::move(fe_vec_obj));
            fe_kinds.push_back(ElementType::FEKind::Vector);
            scalar_fe_slots.push_back(nullptr);
            vector_fe_slots.push_back(raw_fe_vec);
        }
        else
        {
            fe_obj = FEBase::build(fe_data.dimension, fe_type);

            fe_obj->attach_quadrature_rule(qrule_ptr);
            FEBase *raw_fe = fe_obj.get();
            fe_array.push_back(std::move(fe_obj));
            fe_kinds.push_back(ElementType::FEKind::Scalar);
            scalar_fe_slots.push_back(raw_fe);
            vector_fe_slots.push_back(nullptr);
        }
    }
}

libMesh::FEBase &ElementType::scalar_fe(std::size_t idx)
{
    if (idx >= fe_array.size())
        pytec3d::logging::error_and_abort("ElementType::scalar_fe: FE index out of range.");

    if (fe_kinds[idx] != ElementType::FEKind::Scalar || !scalar_fe_slots[idx])
        pytec3d::logging::error_and_abort("ElementType::scalar_fe: requested FE is not scalar.");
    return *scalar_fe_slots[idx];
}

const libMesh::FEBase &ElementType::scalar_fe(std::size_t idx) const
{
    if (idx >= fe_array.size())
        pytec3d::logging::error_and_abort("ElementType::scalar_fe: FE index out of range.");

    if (fe_kinds[idx] != ElementType::FEKind::Scalar || !scalar_fe_slots[idx])
        pytec3d::logging::error_and_abort("ElementType::scalar_fe: requested FE is not scalar.");
    return *scalar_fe_slots[idx];
}

libMesh::FEVectorBase &ElementType::vector_fe(std::size_t idx)
{
    if (idx >= fe_array.size())
        pytec3d::logging::error_and_abort("ElementType::vector_fe: FE index out of range.");

    if (fe_kinds[idx] != ElementType::FEKind::Vector || !vector_fe_slots[idx])
        pytec3d::logging::error_and_abort("ElementType::vector_fe: requested FE is not vector.");
    return *vector_fe_slots[idx];
}

const libMesh::FEVectorBase &ElementType::vector_fe(std::size_t idx) const
{
    if (idx >= fe_array.size())
        pytec3d::logging::error_and_abort("ElementType::vector_fe: FE index out of range.");

    if (fe_kinds[idx] != ElementType::FEKind::Vector || !vector_fe_slots[idx])
        pytec3d::logging::error_and_abort("ElementType::vector_fe: requested FE is not vector.");
    return *vector_fe_slots[idx];
}

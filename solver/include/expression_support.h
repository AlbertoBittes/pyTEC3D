/**
 * @file expression_support.h
 * @brief Stores the result of an "expressions" entry of element_types.json.
 *
 * generate_registries.py compiles each expression as a C++ expression (registries.h, evaluate_expression).
 * A scalar is a libMesh Real, a vector a RealVectorValue, a tensor a RealTensorValue; + - * / and the
 * members of libMesh (cross(), det(), inverse(), transpose(), ...) apply as usual.
 */
#pragma once

#include <type_traits>

#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

#include "parameter_value.h"

namespace pytec3d::expression_support
{
    /** @brief Wraps the result of an expression; it must be of the declared type (a compile error otherwise). */
    template <ParameterType type, typename Value>
    ParameterValue expression_value(const Value &value)
    {
        using T = std::decay_t<Value>;
        if constexpr (type == ParameterType::Scalar)
        {
            static_assert(std::is_arithmetic_v<T>, "the expression must give a scalar (Real)");
            return ParameterValue(static_cast<libMesh::Real>(value));
        }
        else if constexpr (type == ParameterType::Vector)
        {
            static_assert(std::is_base_of_v<libMesh::TypeVector<libMesh::Real>, T>, "the expression must give a vector (RealVectorValue)");
            return ParameterValue(libMesh::RealVectorValue(value));
        }
        else
        {
            static_assert(std::is_base_of_v<libMesh::TypeTensor<libMesh::Real>, T>, "the expression must give a tensor (RealTensorValue)");
            return ParameterValue(libMesh::RealTensorValue(value));
        }
    }
} // namespace pytec3d::expression_support

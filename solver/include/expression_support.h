/**
 * @file expression_support.h
 * @brief Stores the result of an "expressions" entry of element_types.json.
 *
 * generate_registries.py compiles each expression as a C++ expression (registries.h, evaluate_expression).
 * A scalar is a Real, a vector a RealVectorValue, a tensor a RealTensorValue (la.h); + - * / and the
 * members of la.h (cross(), det(), inverse(), transpose(), ...) apply as usual.
 */
#pragma once

#include "la.h"
#include <type_traits>


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
            return ParameterValue(static_cast<pytec3d::la::Real>(value));
        }
        else if constexpr (type == ParameterType::Vector)
        {
            static_assert(std::is_base_of_v<pytec3d::la::RealVectorValue, T>, "the expression must give a vector (RealVectorValue)");
            return ParameterValue(pytec3d::la::RealVectorValue(value));
        }
        else
        {
            static_assert(std::is_base_of_v<pytec3d::la::RealTensorValue, T>, "the expression must give a tensor (RealTensorValue)");
            return ParameterValue(pytec3d::la::RealTensorValue(value));
        }
    }
} // namespace pytec3d::expression_support

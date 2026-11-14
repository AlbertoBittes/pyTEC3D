/**
 * @file parameter_value.h
 * @brief Value of a kernel parameter: a scalar, a 3-vector or a 3x3 tensor, holding only its own type.
 */
#pragma once

#include "la.h"
#include <variant>


namespace pytec3d
{
    enum class ParameterType : unsigned char
    {
        Scalar, ///< pytec3d::la::Real
        Vector, ///< pytec3d::la::RealVectorValue
        Tensor  ///< pytec3d::la::RealTensorValue
    };

    class ParameterValue
    {
    public:
        /** @brief The scalar 0. */
        ParameterValue() = default;
        ParameterValue(pytec3d::la::Real value) : value_(value) {}
        ParameterValue(const pytec3d::la::RealVectorValue &value) : value_(value) {}
        ParameterValue(const pytec3d::la::RealTensorValue &value) : value_(value) {}

        /** @brief The zero of @p type. */
        static ParameterValue zero(ParameterType type);

        /** @brief The order of the alternatives is that of ParameterType. */
        ParameterType type() const { return static_cast<ParameterType>(value_.index()); }

        const pytec3d::la::Real &scalar() const { return std::get<pytec3d::la::Real>(value_); }
        const pytec3d::la::RealVectorValue &vector() const { return std::get<pytec3d::la::RealVectorValue>(value_); }
        const pytec3d::la::RealTensorValue &tensor() const { return std::get<pytec3d::la::RealTensorValue>(value_); }
        pytec3d::la::Real &scalar() { return std::get<pytec3d::la::Real>(value_); }
        pytec3d::la::RealVectorValue &vector() { return std::get<pytec3d::la::RealVectorValue>(value_); }
        pytec3d::la::RealTensorValue &tensor() { return std::get<pytec3d::la::RealTensorValue>(value_); }

        template <typename T>
        const T &get() const { return std::get<T>(value_); }
        template <typename T>
        T &get() { return std::get<T>(value_); }

        /** @brief Sum of two values of the same type. */
        ParameterValue &operator+=(const ParameterValue &other);

        /** @brief The value times @p factor. */
        ParameterValue scaled(pytec3d::la::Real factor) const;

    private:
        std::variant<pytec3d::la::Real, pytec3d::la::RealVectorValue, pytec3d::la::RealTensorValue> value_ = pytec3d::la::Real(0);
    };

    /** @brief "scalar", "vector" or "tensor". */
    const char *type_name(ParameterType type);
} // namespace pytec3d

/**
 * @file parameter_value.h
 * @brief Value of a kernel parameter: a scalar, a 3-vector or a 3x3 tensor, holding only its own type.
 */
#pragma once

#include <variant>

#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

namespace pytec3d
{
    enum class ParameterType : unsigned char
    {
        Scalar, ///< libMesh::Real
        Vector, ///< libMesh::RealVectorValue
        Tensor  ///< libMesh::RealTensorValue
    };

    class ParameterValue
    {
    public:
        /** @brief The scalar 0. */
        ParameterValue() = default;
        ParameterValue(libMesh::Real value) : value_(value) {}
        ParameterValue(const libMesh::RealVectorValue &value) : value_(value) {}
        ParameterValue(const libMesh::RealTensorValue &value) : value_(value) {}

        /** @brief The zero of @p type. */
        static ParameterValue zero(ParameterType type)
        {
            switch (type)
            {
            case ParameterType::Vector:
                return ParameterValue(libMesh::RealVectorValue());
            case ParameterType::Tensor:
                return ParameterValue(libMesh::RealTensorValue());
            case ParameterType::Scalar:
                break;
            }
            return ParameterValue(0.0);
        }

        /** @brief The order of the alternatives is that of ParameterType. */
        ParameterType type() const { return static_cast<ParameterType>(value_.index()); }

        const libMesh::Real &scalar() const { return std::get<libMesh::Real>(value_); }
        const libMesh::RealVectorValue &vector() const { return std::get<libMesh::RealVectorValue>(value_); }
        const libMesh::RealTensorValue &tensor() const { return std::get<libMesh::RealTensorValue>(value_); }
        libMesh::Real &scalar() { return std::get<libMesh::Real>(value_); }
        libMesh::RealVectorValue &vector() { return std::get<libMesh::RealVectorValue>(value_); }
        libMesh::RealTensorValue &tensor() { return std::get<libMesh::RealTensorValue>(value_); }

        template <typename T>
        const T &get() const { return std::get<T>(value_); }
        template <typename T>
        T &get() { return std::get<T>(value_); }

        /** @brief Sum of two values of the same type. */
        ParameterValue &operator+=(const ParameterValue &other)
        {
            std::visit([&other](auto &value)
                       { value += std::get<std::decay_t<decltype(value)>>(other.value_); },
                       value_);
            return *this;
        }

        /** @brief The value times @p factor. */
        ParameterValue scaled(libMesh::Real factor) const
        {
            return std::visit([factor](const auto &value)
                              { return ParameterValue(std::decay_t<decltype(value)>(factor * value)); },
                              value_);
        }

    private:
        std::variant<libMesh::Real, libMesh::RealVectorValue, libMesh::RealTensorValue> value_ = libMesh::Real(0);
    };

    inline const char *type_name(ParameterType type)
    {
        switch (type)
        {
        case ParameterType::Scalar:
            return "scalar";
        case ParameterType::Vector:
            return "vector";
        case ParameterType::Tensor:
            return "tensor";
        }
        return "scalar";
    }
} // namespace pytec3d

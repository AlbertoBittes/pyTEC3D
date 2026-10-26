/**
 * @file parameter_value.cpp
 * @brief ParameterValue: zero values, sums and scaling of scalar, vector and tensor parameters.
 */
#include "parameter_value.h"

#include <type_traits>

namespace pytec3d
{
    ParameterValue ParameterValue::zero(ParameterType type)
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

    ParameterValue &ParameterValue::operator+=(const ParameterValue &other)
    {
        std::visit([&other](auto &value)
                   { value += std::get<std::decay_t<decltype(value)>>(other.value_); },
                   value_);
        return *this;
    }

    ParameterValue ParameterValue::scaled(libMesh::Real factor) const
    {
        return std::visit([factor](const auto &value)
                          { return ParameterValue(std::decay_t<decltype(value)>(factor * value)); },
                          value_);
    }

    const char *type_name(ParameterType type)
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

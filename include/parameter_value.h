/**
 * @file parameter_value.h
 * @brief Value of a kernel parameter: a scalar, a 3-vector or a 3x3 tensor.
 */
#pragma once

#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

namespace pytec3d
{
    enum class ParameterType : unsigned char
    {
        Scalar, ///< "name"
        Vector, ///< "name[]"
        Tensor, ///< "name[][]"
        Unknown ///< "name()": a grid array whose type is its number of components.
    };

    struct ParameterValue
    {
        ParameterType type = ParameterType::Scalar;
        libMesh::Real scalar = 0.0;
        libMesh::VectorValue<libMesh::Real> vector;
        libMesh::TensorValue<libMesh::Real> tensor;
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
        case ParameterType::Unknown:
            break;
        }
        return "array";
    }

    inline ParameterValue scalar_value(libMesh::Real value)
    {
        ParameterValue out;
        out.scalar = value;
        return out;
    }

    inline ParameterValue scaled(const ParameterValue &value, libMesh::Real factor)
    {
        ParameterValue out = value;
        out.scalar *= factor;
        out.vector *= factor;
        out.tensor *= factor;
        return out;
    }
} // namespace pytec3d

/**
 * @file expression.h
 * @brief Mathematical expressions that define derived parameters from other parameters.
 *
 * An expression is parsed once and evaluated for every element. Values are scalars, vectors and tensors
 * (parameter_value.h); the result type is checked when the expression is evaluated.
 *
 *  - Operators: + - * / and parentheses. A scalar scales anything; tensor * vector and tensor * tensor
 *    are matrix products; vector * vector is an error (use dot() or outer()).
 *  - Parameters are written with their type suffix: name, name[], name[][], name().
 *  - Constants: pi, mu0, dimension (2 or 3), axisymmetric (0 or 1), ex, ey, ez (unit vectors).
 *  - Functions: abs(s) sqrt(s) cbrt(s) det(T) diag(v) outer(a, b) dot(a, b)
 *               inv(s|v|T)     reciprocal, element by element for a vector; inverse of a tensor
 *               invpos(s|v|T)  reciprocal where positive and finite, else 0 (a tensor: per entry if diagonal)
 *               global(v|T)  from the element's local axes to global axes
 *               select(c, a, b)  a if c is not zero, else b (only the one taken is evaluated)
 */
#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include <libmesh/tensor_value.h>

#include "parameter_value.h"

namespace pytec3d
{
    /** @brief What an expression reads from the element it is evaluated on. */
    class ExpressionContext
    {
    public:
        virtual ~ExpressionContext() = default;
        virtual const ParameterValue &symbol(unsigned int index) = 0;
        virtual const libMesh::TensorValue<libMesh::Real> &rotation() const = 0; ///< Local to global axes.
    };

    class Expression
    {
    public:
        /** @brief Maps a parameter as written ("nu[][]") to its index; nullopt when unknown. */
        using Resolver = std::function<std::optional<unsigned int>(const std::string &)>;

        /** @brief Parses @p text; aborts on any error. @p context names the expression in messages. */
        static Expression parse(const std::string &text,
                                const std::string &context,
                                const Resolver &resolve,
                                unsigned int dimension,
                                bool axisymmetric);

        ParameterValue evaluate(ExpressionContext &context) const;

    private:
        class Builder;

        enum class Op : unsigned char
        {
            Constant,
            Symbol,
            Negate,
            Add,
            Subtract,
            Multiply,
            Divide,
            Call
        };

        enum class Function : unsigned char
        {
            Abs,
            Sqrt,
            Cbrt,
            Det,
            Diag,
            Outer,
            Dot,
            Inv,
            Invpos,
            Global,
            Select
        };

        struct Node
        {
            Op op = Op::Constant;
            ParameterValue constant;
            unsigned int symbol = 0;
            Function function = Function::Abs;
            int a = -1;
            int b = -1;
            int c = -1;
        };

        ParameterValue evaluate_node(int index, ExpressionContext &context) const;

        std::vector<Node> nodes_;
        int root_ = -1;
        std::string description_;
    };
} // namespace pytec3d

/**
 * @file nabla.h
 * @brief Differential operators (identity, grad, curl, div) evaluated on a field's finite element.
 *
 * Nabla turns "operator applied to a variable" into one dense matrix per quadrature point. Rows are
 * the components of the operator output, columns are the element shape functions of the variable in
 * libMesh's order (for a vector variable phi_x_1 phi_y_1 phi_z_1 ... phi_x_n phi_y_n phi_z_n). Nabla knows
 * nothing about where a variable sits in the element matrix; combining operands with coefficients and
 * placing the resulting block at the right DOF indices is the job of ElementKernels.
 *
 * Coordinate conventions:
 *  - 2D elements lie in the x-y plane; a scalar that belongs to a vector field with the single
 *    active component z is the out-of-plane component (A_z, or A_theta when axisymmetric).
 *  - A scalar that is not a vector component (V, VOLT, EMF) on a 2D element is a lumped potential
 *    over the out-of-plane length L = thickness (times 2*pi*r if axisymmetric), so its gradient is
 *    the single z row -phi / L.
 *  - Axisymmetric elements use (x, y) = (r, z) and add the cylindrical terms to curl and div.
 */
#pragma once

#include <string>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/libmesh_common.h>

namespace libMesh
{
    class Elem;
}

namespace pytec3d
{
    class ElementTypeInstance;

    enum class NablaOperator : unsigned char
    {
        Identity, ///< The field itself.
        Grad,
        Curl,
        Div,
    };

    /** @brief A parsed operator expression such as "curl(dt(A))" or "grad(V)". */
    struct OperatorExpression
    {
        NablaOperator op = NablaOperator::Identity;
        bool time_derivative = false; ///< dt(...) present; selects the time-derivative matrix, not the shape data.
        std::string variable;

        /** @brief Parses @p text; dt() may wrap or be wrapped by one spatial operator. Aborts on malformed input. */
        static OperatorExpression parse(const std::string &text);
    };

    /** @brief One signed term of a field expression: sign * [coefficient *] operator(variable). */
    struct ExpressionTerm
    {
        libMesh::Real sign = 1.0;
        std::string coefficient; ///< Empty when the term has no coefficient factor.
        OperatorExpression operand;

        /**
         * @brief Parses a sum such as "-sigma[][]*dt(A) + nu[][]*curl(A)".
         *
         * Factors named in @p coefficient_names are coefficients (at most one per term); the other factor
         * is the operator expression. Aborts on malformed input.
         */
        static std::vector<ExpressionTerm> parse_sum(const std::string &text,
                                                     const std::vector<std::string> &coefficient_names);
    };

    class Nabla
    {
    public:
        /** @brief One (rows x shape functions) matrix per quadrature point. */
        using QuadratureMatrices = std::vector<libMesh::DenseMatrix<libMesh::Real>>;

        /** @brief Operator matrices plus the physical component (0=x, 1=y, 2=z) each row stands for. */
        struct Operand
        {
            QuadratureMatrices values;
            std::vector<unsigned int> components;
        };

        Nabla() = delete;

        /**
         * @brief Evaluates @p op on @p variable over @p elem.
         *
         * Output rows by case:
         *  - Identity: 1 for scalar variables, one per active component for vector variables.
         *  - Grad of a scalar: element dimension (a single z row for 2D lumped potentials).
         *  - Grad of a vector (LAGRANGE_VEC families, not axisymmetric): the tensor flattened row-major,
         *    row = component * dimension + direction, entry d phi_component / d x_direction.
         *  - Curl: 2 for an out-of-plane scalar in 2D, 3 for a vector variable in 3D.
         *  - Div: 1 (zero for an out-of-plane scalar in 2D).
         *
         * @p thickness is the out-of-plane length of planar 2D elements.
         */
        static Operand operand(const ElementTypeInstance &element_type,
                               const std::string &variable,
                               NablaOperator op,
                               const libMesh::Elem &elem,
                               libMesh::Real thickness = 1.0);

        /** @brief Same as operand(), keeping only the matrices. */
        static QuadratureMatrices evaluate(const ElementTypeInstance &element_type,
                                           const std::string &variable,
                                           NablaOperator op,
                                           const libMesh::Elem &elem,
                                           libMesh::Real thickness = 1.0);

        /** @brief Convenience overload; the time-derivative flag of @p expression is ignored here. */
        static QuadratureMatrices evaluate(const ElementTypeInstance &element_type,
                                           const OperatorExpression &expression,
                                           const libMesh::Elem &elem,
                                           libMesh::Real thickness = 1.0);

        /** @brief Quadrature weights: JxW, times thickness (2D) and 2*pi*r (axisymmetric). */
        static std::vector<libMesh::Real> integration_weights(const ElementTypeInstance &element_type,
                                                              const libMesh::Elem &elem,
                                                              libMesh::Real thickness = 1.0);
    };
} // namespace pytec3d

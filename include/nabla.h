/**
 * @file nabla.h
 * @brief Differential operators (identity, grad, curl, div) evaluated on a field's finite element.
 *
 * A Nabla object is one operator bound to one field: built once, it checks that the operator exists for the
 * field, asks the field's FE for the data it reads, and evaluates on whatever element the FE was last
 * reinit()ed on. Evaluating gives one dense matrix per quadrature point. Rows are the components of the operator output, columns are the shape
 * functions of the field in libMesh's order (for a vector field phi_x_1 phi_y_1 phi_z_1 ... phi_x_n
 * phi_y_n phi_z_n). Where a field sits in an element matrix, and every other piece of physics, is not
 * Nabla's concern.
 *
 * Conventions:
 *  - 2D elements lie in the x-y plane; on them the curl of a scalar is that of the out-of-plane vector
 *    f e_z, and the divergence of a scalar is zero.
 *  - Axisymmetric elements use (x, y) = (r, z) and add the cylindrical terms to curl and div.
 *
 * The operators themselves are the function templates of nabla_operators.h.
 */
#pragma once

#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/libmesh_common.h>

namespace libMesh
{
    class FEAbstract;
}

namespace pytec3d
{
    struct Field;

    enum class NablaOperator : unsigned char
    {
        Identity, ///< The field itself.
        Grad,
        Curl,
        Div,
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
            bool scalar = false; ///< A scalar quantity (one row, component x): a vector parameter couples it to a vector operand.
        };

        /** @brief An operator compiled for one output type, dimension and axisymmetry. */
        using Function = Operand (*)(libMesh::FEAbstract &fe, const Field &field);

        /**
         * @brief Binds @p op to @p field, whose FE must be attached; aborts if the operator does not exist for it.
         *
         * Output rows by operator:
         *  - Identity: 1 for a scalar field, one per dimension for a vector field.
         *  - Grad of a scalar: one row per dimension. Of a vector (Lagrange families, not axisymmetric): the
         *    tensor flattened row-major, row = component * dimension + direction.
         *  - Curl: 2 for a scalar in 2D, 3 for a vector in 3D.
         *  - Div: 1 (zero for a scalar in 2D).
         */
        Nabla(NablaOperator op, const Field &field, unsigned int dimension, bool axisymmetric);

        NablaOperator op() const { return op_; }
        const Field &field() const { return *field_; }

        /** @brief True when the output is a scalar quantity: the identity of a scalar field, or a divergence. */
        bool scalar_valued() const;

        /** @brief The operator on the element the field's FE was last reinit()ed on. */
        Operand evaluate() const;

    private:
        const Field *field_;
        NablaOperator op_;
        Function function_;
    };
} // namespace pytec3d

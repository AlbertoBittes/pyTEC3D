/**
 * @file element_evaluator.h
 * @brief ElementEvaluator: the evaluations of an element type on one element.
 *
 * Built for one element and one operation (it reinit()s only the FEs that operation needs), it computes on
 * demand, once each:
 *  - the nabla tables of the FEs;
 *  - the coefficients: parameters (ElementParameterValues) and law outputs, one per quadrature point, from
 *    the law's input field at the quadrature points (current solution) and the curve of the element material;
 *  - the kernels (kernels.h).
 * and places them as the evaluations of the spec say (ElementType::uses() decides which terms an operation
 * computes). The laws selected on the element (Inputs::laws) replace the terms whose coefficient they replace
 * by their variants (formulation_spec.h), with the curves of the rows that select them.
 *
 *   Block terms   Linear (K) -> K, C with a dt() trial, M with a dt(dt()) trial. Law variants: Tangent -> Jacobian
 *                 (times c_d, d the dt() on the trial), Residual -> residual (+= block * d^d u / dt^d).
 *   Rhs terms     Linear (F) -> right-hand side. Law variants: Residual -> residual (+= vector).
 *   Quadrature    a field at each quadrature point: sum of the terms' C * trial(u) (or C alone); a law input.
 *   L2Projection  M^-1 * sum of the terms' kernels (tested with the field) applied to u.
 *   Average       mean over the quadrature points of the terms, as Quadrature.
 *
 * A nodal or element field is static where none of its terms is a law variant used on the element: it is then a
 * linear operator on the solution, F u + F_dt du/dt + ..., computed once with the element matrices
 * (field_operators) and only applied at output. Dynamic fields are evaluated at output (field_values).
 *
 * Integration weights are JxW, times the out-of-plane length L of 2D elements: the parameter the element type
 * declares as its out_of_plane_length (1 if none), times 2 pi r when axisymmetric. In field terms, a dt()
 * (dt(dt())) trial acts on the first (second) time derivative of the solution.
 *
 * The work is done by an implementation compiled for each element type (element_evaluator.cpp): its terms,
 * kernels and their value types are compile-time constants of its spec, so every kernel is instantiated with
 * its own test, coefficient and trial types and nothing is dispatched at run time.
 */
#pragma once

#include <array>
#include <memory>

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/elem.h>

#include "element_parameters.h"
#include "element_type.h"

namespace pytec3d
{
    class ModelData;

    class ElementEvaluator
    {
    public:
        using Matrix = libMesh::DenseMatrix<libMesh::Number>;
        using Vector = libMesh::DenseVector<libMesh::Number>;

        struct Inputs
        {
            const ElementType *type = nullptr;
            registry::Domain domain = registry::Domain::COUNT;
            Operation operation = Operation::Parameters;          ///< What will be computed: only its FEs are reinit()ed.
            registry::LawMask laws = 0;                          ///< The laws selected on the element (model data "law").
            const libMesh::Elem *elem = nullptr;
            ElementParameterValues *parameters = nullptr;
            const ModelData *model_data = nullptr;
            unsigned int material_id = 0;                        ///< MATERIAL, CONSTANT and TYPE ids: the rows of the
            unsigned int constants_id = 0;                       ///< laws' curves, by the source of each law.
            unsigned int type_id = 0;
            const ElementType::FieldIndices *dofs = nullptr;     ///< Of the DOF fields, in the DOF system element vectors.
            const Vector *solution = nullptr;                    ///< Local DOF solution; needed by laws, residual and fields.
            const Vector *solution_dt = nullptr;                 ///< Its time derivative; needed by dt() field terms.
            const Vector *solution_ddt = nullptr;                ///< Its second time derivative; needed by dt(dt()) field terms.
            /// d(d^k u / dt^k)/du of the time scheme: {1, c_v, c_a}; the law Jacobian blocks with k dt() are scaled by it.
            std::array<libMesh::Real, 3> jacobian_scale = {1.0, 0.0, 0.0};
        };

        explicit ElementEvaluator(const Inputs &inputs);
        ~ElementEvaluator();

        /** @brief K, C and M of the element: the linear block terms by number of dt() (and the replaced ones on linear
         * elements). @p M is empty for element types of time order below 2. */
        void add_matrices(Matrix &K, Matrix &C, Matrix &M) { impl_->add_matrices(K, C, M); }

        /** @brief Right-hand side of the element: the linear rhs terms (and the replaced ones on linear elements). */
        void add_rhs(Vector &F) { impl_->add_rhs(F); }

        /** @brief Jacobian (J) and residual (R) terms of the laws selected on the element, at the current solution. */
        void add_nonlinear(Matrix &J, Vector &R) { impl_->add_nonlinear(J, R); }

        /**
         * @brief Operators of the static nodal or element fields (Matrices operation): the field values are
         * F[0] u + F[1] du/dt + F[2] d2u/dt2 + f, F of size (field system element DOFs, DOF system element DOFs).
         * All the fields of a kind collapse into these three matrices and one vector; those the spec cannot fill
         * on the element (ElementType::field_operators) are left empty.
         */
        void field_operators(FieldKind kind, const ElementType::FieldIndices &indices, std::array<Matrix, 3> &F, Vector &f)
        {
            impl_->field_operators(kind, indices, F, f);
        }

        /** @brief Values of the dynamic nodal or element fields, written at their positions @p indices in @p values. */
        void field_values(FieldKind kind, const ElementType::FieldIndices &indices, Vector &values)
        {
            impl_->field_values(kind, indices, values);
        }

        /** @brief Reads every parameter the terms active in the domain need, so that a missing one aborts now. */
        void check_parameters() { impl_->check_parameters(); }

        /** @brief The evaluator of one element type (element_evaluator.cpp). */
        class Implementation
        {
        public:
            virtual ~Implementation() = default;
            virtual void add_matrices(Matrix &K, Matrix &C, Matrix &M) = 0;
            virtual void add_rhs(Vector &F) = 0;
            virtual void add_nonlinear(Matrix &J, Vector &R) = 0;
            virtual void field_operators(FieldKind kind, const ElementType::FieldIndices &indices, std::array<Matrix, 3> &F,
                                         Vector &f) = 0;
            virtual void field_values(FieldKind kind, const ElementType::FieldIndices &indices, Vector &values) = 0;
            virtual void check_parameters() = 0;
        };

    private:
        std::unique_ptr<Implementation> impl_;
    };
} // namespace pytec3d

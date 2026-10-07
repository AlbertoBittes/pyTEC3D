/**
 * @file kernels.h
 * @brief Weak-form kernels: pure functions that combine operator matrices into element matrices.
 *
 * A kernel knows nothing but its inputs: a parameter (Real, RealVectorValue or TensorValue), the operator
 * matrices of its trial and/or test, the quadrature weights and the element matrix indices where the
 * result goes. There are four, plus a nonlinear version of the DOF block and of the two field kernels
 * in which the parameter is a law H = law(trial) with a Jacobian:
 *
 *  - Block:        K(test_ids, trial_ids) += integral of test^T * C * trial   (DOF matrices K and C)
 *  - Source:       rhs(test_ids)          += integral of test^T * C           (DOF right-hand side)
 *  - NodalField:   F(field_ids, trial_ids) += M^-1 * integral of basis^T * C * trial
 *                  the L2 projection onto the nodal basis, as a matrix acting on the DOF solution
 *  - ElementField: F(field_ids, trial_ids) += mean over the quadrature points of C * trial
 *
 * Operator rows are physical components (Nabla::Operand::components). A tensor parameter couples
 * component a of one side to component b of the other with C(a, b); a scalar parameter multiplies row
 * by row (Block) or component by component (fields). A vector parameter couples a scalar operand
 * (Nabla::Operand::scalar) to a vector one: component a of the vector side gets C(a) times the scalar.
 */
#pragma once

#include <string>
#include <type_traits>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

#include "logging.h"
#include "nabla.h"

namespace pytec3d::kernel
{
    using Real = libMesh::Real;
    using Matrix = libMesh::DenseMatrix<Real>;
    using Vector = libMesh::DenseVector<Real>;
    using Ids = std::vector<unsigned int>;
    using Components = std::vector<unsigned int>;
    using Tensor = libMesh::TensorValue<Real>;
    using Vec3 = libMesh::VectorValue<Real>;

    namespace detail
    {
        [[noreturn]] inline void fail(const std::string &message)
        {
            logging::error_and_abort("Kernel: " + message);
        }

        template <typename Parameter>
        constexpr bool is_scalar = std::is_same_v<Parameter, Real>;
        template <typename Parameter>
        constexpr bool is_vector = std::is_same_v<Parameter, Vec3>;
        template <typename Parameter>
        constexpr bool is_tensor = std::is_same_v<Parameter, Tensor>;

        inline void check_physical(const Components &components, const char *what)
        {
            for (unsigned int c : components)
                if (c > 2)
                    fail(std::string(what) + " has no physical components (x, y, z) to couple with a parameter.");
        }

        // V = parameter * operand, with one row per entry of row_components.
        template <typename Parameter>
        void apply_parameter(const Parameter &parameter,
                             const Matrix &operand,
                             const Components &operand_components,
                             const Components &row_components,
                             Matrix &V)
        {
            V.resize(static_cast<unsigned int>(row_components.size()), operand.n());
            for (unsigned int r = 0; r < V.m(); ++r)
                for (unsigned int j = 0; j < operand.n(); ++j)
                {
                    Real sum = 0.0;
                    for (unsigned int b = 0; b < operand.m(); ++b)
                    {
                        if constexpr (is_tensor<Parameter>)
                            sum += parameter(row_components[r], operand_components[b]) * operand(b, j);
                        else if constexpr (is_vector<Parameter>)
                            sum += parameter(row_components[r]) * operand(b, j);
                        else if (operand_components[b] == row_components[r])
                            sum += parameter * operand(b, j);
                    }
                    V(r, j) = sum;
                }
        }
    } // namespace detail

    /** @brief K(test_ids, trial_ids) += integral of test^T * parameter * trial. */
    template <typename Parameter>
    void Block(const Parameter &parameter,
               const Nabla::Operand &trial,
               const Nabla::Operand &test,
               const std::vector<Real> &weights,
               const Ids &trial_ids,
               const Ids &test_ids,
               Matrix &K)
    {
        if (test.values.size() != trial.values.size() || test.values.size() != weights.size())
            detail::fail("Block: test, trial and weights must have one entry per quadrature point.");
        if (test.values.empty())
            return;

        const unsigned int rows_test = test.values[0].m();
        const unsigned int rows_trial = trial.values[0].m();
        if (test.values[0].n() != test_ids.size() || trial.values[0].n() != trial_ids.size())
            detail::fail("Block: the operators do not match the DOF indices of their fields.");
        if constexpr (detail::is_tensor<Parameter>)
        {
            detail::check_physical(test.components, "Block: the test operator");
            detail::check_physical(trial.components, "Block: the trial operator");
        }
        else if constexpr (detail::is_vector<Parameter>)
        {
            if (test.scalar == trial.scalar)
                detail::fail("Block: a vector parameter couples a scalar operand with a vector operand.");
            detail::check_physical(test.scalar ? trial.components : test.components, "Block: the vector operand");
        }
        else if (rows_test != rows_trial)
            detail::fail("Block: a scalar parameter needs test and trial with the same number of rows.");

        Matrix block(test.values[0].n(), trial.values[0].n());
        std::vector<Real> cu(rows_test);
        for (std::size_t qp = 0; qp < test.values.size(); ++qp)
        {
            const Matrix &T = test.values[qp];
            const Matrix &U = trial.values[qp];
            const Real w = weights[qp];

            for (unsigned int j = 0; j < U.n(); ++j)
            {
                for (unsigned int a = 0; a < rows_test; ++a)
                {
                    cu[a] = 0.0;
                    if constexpr (detail::is_tensor<Parameter>)
                        for (unsigned int b = 0; b < rows_trial; ++b)
                            cu[a] += parameter(test.components[a], trial.components[b]) * U(b, j);
                    else if constexpr (detail::is_vector<Parameter>)
                    {
                        if (trial.scalar)
                            cu[a] = parameter(test.components[a]) * U(0, j);
                        else
                            for (unsigned int b = 0; b < rows_trial; ++b)
                                cu[a] += parameter(trial.components[b]) * U(b, j);
                    }
                    else
                        cu[a] = parameter * U(a, j);
                }
                for (unsigned int i = 0; i < T.n(); ++i)
                {
                    Real sum = 0.0;
                    for (unsigned int a = 0; a < rows_test; ++a)
                        sum += T(a, i) * cu[a];
                    block(i, j) += w * sum;
                }
            }
        }

        for (unsigned int i = 0; i < block.m(); ++i)
            for (unsigned int j = 0; j < block.n(); ++j)
                K(test_ids[i], trial_ids[j]) += block(i, j);
    }

    /** @brief rhs(test_ids) += integral of test^T * parameter. */
    template <typename Parameter>
    void Source(const Parameter &parameter,
                const Nabla::Operand &test,
                const std::vector<Real> &weights,
                const Ids &test_ids,
                Vector &rhs)
    {
        static_assert(detail::is_vector<Parameter> || detail::is_scalar<Parameter>,
                      "a source parameter is a vector or a scalar.");
        if (test.values.size() != weights.size())
            detail::fail("Source: test and weights must have one entry per quadrature point.");
        if (test.values.empty())
            return;
        if (test.values[0].n() != test_ids.size())
            detail::fail("Source: the operator does not match the DOF indices of its field.");
        if constexpr (detail::is_vector<Parameter>)
            detail::check_physical(test.components, "Source: the test operator");

        for (std::size_t qp = 0; qp < test.values.size(); ++qp)
        {
            const Matrix &T = test.values[qp];
            for (unsigned int i = 0; i < T.n(); ++i)
            {
                Real sum = 0.0;
                for (unsigned int a = 0; a < T.m(); ++a)
                {
                    if constexpr (detail::is_vector<Parameter>)
                        sum += T(a, i) * parameter(test.components[a]);
                    else
                        sum += T(a, i) * parameter;
                }
                rhs(test_ids[i]) += weights[qp] * sum;
            }
        }
    }

    /** @brief F(field_ids, trial_ids) += M^-1 * integral of basis^T * parameter * trial, M the basis mass matrix. */
    template <typename Parameter>
    void NodalField(const Parameter &parameter,
                    const Nabla::Operand &trial,
                    const Nabla::Operand &basis,
                    const std::vector<Real> &weights,
                    const Ids &trial_ids,
                    const Ids &field_ids,
                    Matrix &F)
    {
        if (basis.values.size() != trial.values.size() || basis.values.size() != weights.size())
            detail::fail("NodalField: basis, trial and weights must have one entry per quadrature point.");
        if (basis.values.empty())
            return;

        const unsigned int n = basis.values[0].n();
        const unsigned int n_trial = trial.values[0].n();
        if (n != field_ids.size() || n_trial != trial_ids.size())
            detail::fail("NodalField: the operators do not match the DOF indices of their fields.");
        if constexpr (detail::is_tensor<Parameter>)
        {
            detail::check_physical(trial.components, "NodalField: the trial operator");
            detail::check_physical(basis.components, "NodalField: the field");
        }
        else if constexpr (detail::is_vector<Parameter>)
        {
            if (!trial.scalar)
                detail::fail("NodalField: a vector parameter needs a scalar trial.");
            detail::check_physical(basis.components, "NodalField: the field");
        }

        Matrix mass(n, n);
        Matrix G(n, n_trial);
        Matrix V;
        for (std::size_t qp = 0; qp < basis.values.size(); ++qp)
        {
            const Matrix &N = basis.values[qp];
            const Real w = weights[qp];
            detail::apply_parameter(parameter, trial.values[qp], trial.components, basis.components, V);

            for (unsigned int i = 0; i < n; ++i)
            {
                for (unsigned int k = 0; k < n; ++k)
                {
                    Real sum = 0.0;
                    for (unsigned int r = 0; r < N.m(); ++r)
                        sum += N(r, i) * N(r, k);
                    mass(i, k) += w * sum;
                }
                for (unsigned int j = 0; j < n_trial; ++j)
                {
                    Real sum = 0.0;
                    for (unsigned int r = 0; r < N.m(); ++r)
                        sum += N(r, i) * V(r, j);
                    G(i, j) += w * sum;
                }
            }
        }

        Vector rhs(n);
        Vector x(n);
        for (unsigned int j = 0; j < n_trial; ++j)
        {
            for (unsigned int i = 0; i < n; ++i)
                rhs(i) = G(i, j);
            mass.cholesky_solve(rhs, x);
            for (unsigned int i = 0; i < n; ++i)
                F(field_ids[i], trial_ids[j]) += x(i);
        }
    }

    /** @brief F(field_ids, trial_ids) += mean over the quadrature points of parameter * trial, one row per field component. */
    template <typename Parameter>
    void ElementField(const Parameter &parameter,
                      const Nabla::Operand &trial,
                      const Components &field_components,
                      const Ids &trial_ids,
                      const Ids &field_ids,
                      Matrix &F)
    {
        if (trial.values.empty())
            return;
        if (field_ids.size() != field_components.size() || trial.values[0].n() != trial_ids.size())
            detail::fail("ElementField: the operator does not match the DOF indices of its field.");
        if constexpr (detail::is_tensor<Parameter>)
        {
            detail::check_physical(trial.components, "ElementField: the trial operator");
            detail::check_physical(field_components, "ElementField: the field");
        }
        else if constexpr (detail::is_vector<Parameter>)
        {
            if (!trial.scalar)
                detail::fail("ElementField: a vector parameter needs a scalar trial.");
            detail::check_physical(field_components, "ElementField: the field");
        }

        const Real inverse_count = 1.0 / static_cast<Real>(trial.values.size());
        Matrix V;
        for (std::size_t qp = 0; qp < trial.values.size(); ++qp)
        {
            detail::apply_parameter(parameter, trial.values[qp], trial.components, field_components, V);
            for (unsigned int r = 0; r < V.m(); ++r)
                for (unsigned int j = 0; j < V.n(); ++j)
                    F(field_ids[r], trial_ids[j]) += inverse_count * V(r, j);
        }
    }

    // A law maps the trial value B (3 physical components) to H and its Jacobian dH/dB: law(B, H, J).

    namespace detail
    {
        // Trial value at one quadrature point, from the local solution u.
        inline void trial_value(const Nabla::Operand &trial, std::size_t qp, const Ids &trial_ids, const Vector &u, Real B[3])
        {
            const Matrix &U = trial.values[qp];
            B[0] = B[1] = B[2] = 0.0;
            for (unsigned int r = 0; r < U.m(); ++r)
                for (unsigned int j = 0; j < U.n(); ++j)
                    B[trial.components[r]] += U(r, j) * u(trial_ids[j]);
        }
    } // namespace detail

    /** @brief tangent(test_ids, trial_ids) += integral of test^T * J * trial; residual(test_ids) += integral of test^T * H. */
    template <typename Law>
    void NonlinearBlock(const Law &law,
                        const Nabla::Operand &trial,
                        const Nabla::Operand &test,
                        const std::vector<Real> &weights,
                        const Ids &trial_ids,
                        const Ids &test_ids,
                        const Vector &u,
                        Matrix &tangent,
                        Vector &residual)
    {
        if (test.values.size() != trial.values.size() || test.values.size() != weights.size())
            detail::fail("NonlinearBlock: test, trial and weights must have one entry per quadrature point.");
        if (test.values.empty())
            return;
        if (test.values[0].n() != test_ids.size() || trial.values[0].n() != trial_ids.size())
            detail::fail("NonlinearBlock: the operators do not match the DOF indices of their fields.");
        detail::check_physical(test.components, "NonlinearBlock: the test operator");
        detail::check_physical(trial.components, "NonlinearBlock: the trial operator");

        for (std::size_t qp = 0; qp < weights.size(); ++qp)
        {
            const Matrix &T = test.values[qp];
            const Matrix &U = trial.values[qp];
            const Real w = weights[qp];

            Real B[3];
            Real H[3];
            Real J[3][3];
            detail::trial_value(trial, qp, trial_ids, u, B);
            law(B, H, J);

            for (unsigned int i = 0; i < T.n(); ++i)
            {
                Real h_term = 0.0;
                for (unsigned int a = 0; a < T.m(); ++a)
                    h_term += T(a, i) * H[test.components[a]];
                residual(test_ids[i]) += w * h_term;

                for (unsigned int j = 0; j < U.n(); ++j)
                {
                    Real tangent_term = 0.0;
                    for (unsigned int a = 0; a < T.m(); ++a)
                        for (unsigned int b = 0; b < U.m(); ++b)
                            tangent_term += T(a, i) * J[test.components[a]][trial.components[b]] * U(b, j);
                    tangent(test_ids[i], trial_ids[j]) += w * tangent_term;
                }
            }
        }
    }

    /** @brief values(field_ids) += M^-1 * integral of basis^T * H(trial(u)), M the basis mass matrix. */
    template <typename Law>
    void NonlinearNodalField(const Law &law,
                             const Nabla::Operand &trial,
                             const Nabla::Operand &basis,
                             const std::vector<Real> &weights,
                             const Ids &trial_ids,
                             const Ids &field_ids,
                             const Vector &u,
                             Vector &values)
    {
        if (basis.values.size() != trial.values.size() || basis.values.size() != weights.size())
            detail::fail("NonlinearNodalField: basis, trial and weights must have one entry per quadrature point.");
        if (basis.values.empty())
            return;
        const unsigned int n = basis.values[0].n();
        if (n != field_ids.size() || trial.values[0].n() != trial_ids.size())
            detail::fail("NonlinearNodalField: the operators do not match the DOF indices of their fields.");
        detail::check_physical(trial.components, "NonlinearNodalField: the trial operator");
        detail::check_physical(basis.components, "NonlinearNodalField: the field");

        Matrix mass(n, n);
        Vector rhs(n);
        for (std::size_t qp = 0; qp < basis.values.size(); ++qp)
        {
            const Matrix &N = basis.values[qp];
            Real B[3];
            Real H[3];
            Real J[3][3];
            detail::trial_value(trial, qp, trial_ids, u, B);
            law(B, H, J);

            for (unsigned int i = 0; i < n; ++i)
            {
                for (unsigned int k = 0; k < n; ++k)
                {
                    Real sum = 0.0;
                    for (unsigned int r = 0; r < N.m(); ++r)
                        sum += N(r, i) * N(r, k);
                    mass(i, k) += weights[qp] * sum;
                }
                Real sum = 0.0;
                for (unsigned int r = 0; r < N.m(); ++r)
                    sum += N(r, i) * H[basis.components[r]];
                rhs(i) += weights[qp] * sum;
            }
        }

        Vector x(n);
        mass.cholesky_solve(rhs, x);
        for (unsigned int i = 0; i < n; ++i)
            values(field_ids[i]) += x(i);
    }

    /** @brief values(field_ids) += mean over the quadrature points of H(trial(u)), one entry per field component. */
    template <typename Law>
    void NonlinearElementField(const Law &law,
                               const Nabla::Operand &trial,
                               const Components &field_components,
                               const Ids &trial_ids,
                               const Ids &field_ids,
                               const Vector &u,
                               Vector &values)
    {
        if (trial.values.empty())
            return;
        if (field_ids.size() != field_components.size() || trial.values[0].n() != trial_ids.size())
            detail::fail("NonlinearElementField: the operator does not match the DOF indices of its field.");
        detail::check_physical(trial.components, "NonlinearElementField: the trial operator");
        detail::check_physical(field_components, "NonlinearElementField: the field");

        const Real inverse_count = 1.0 / static_cast<Real>(trial.values.size());
        for (std::size_t qp = 0; qp < trial.values.size(); ++qp)
        {
            Real B[3];
            Real H[3];
            Real J[3][3];
            detail::trial_value(trial, qp, trial_ids, u, B);
            law(B, H, J);
            for (unsigned int r = 0; r < field_components.size(); ++r)
                values(field_ids[r]) += inverse_count * H[field_components[r]];
        }
    }
} // namespace pytec3d::kernel

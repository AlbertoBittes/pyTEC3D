/**
 * @file kernels.h
 * @brief Kernels: test^T * coefficient * trial over the quadrature points of one element.
 *
 * A kernel only sees Nabla tables, a coefficient and the integration weights; which fields it couples and
 * where the result goes is the business of the evaluation. Four shapes, by what the kernel has:
 *
 *   test and trial  Matrix(i, j) = sum_qp w_qp * dot(test_i, value(C, trial_j))   (blocks of K, C, M, mass matrices)
 *                   or, applied to a vector u, F(i) = sum_qp w_qp * dot(test_i, sum_j value(C, trial_j) u_j)
 *   test only       Vector(i)    = sum_qp w_qp * dot(test_i, C)                     (right-hand side)
 *   trial only      Table[qp][j] = value(C, trial_j) at each quadrature point, no weight (like the JVe matrices)
 *   neither         Values[qp]   = C at each quadrature point
 *
 * value(C, u) is computed once per quadrature point and trial shape function, then dotted with every test
 * shape function. The weights are JxW (times the out-of-plane length of 2D elements); the coefficient is
 * constant over the element or has one value per quadrature point (law outputs).
 *
 *   value(C, u):  Real*Real  Real*Vector  Real*Tensor  Tensor*Vector  Vector*Real  Vector.Vector (Real)
 *   dot(t, x):    Real*Real  Vector.Vector  Tensor:Tensor

 */
#pragma once

#include <type_traits>
#include <utility>
#include <vector>

#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/tensor_value.h>
#include <libmesh/vector_value.h>

#include "nabla.h"

namespace pytec3d::kernel
{
    using Real = libMesh::Real;
    using Vec = libMesh::RealVectorValue;
    using Ten = libMesh::RealTensorValue;

    /** @brief A coefficient: one value, or one per quadrature point. */
    template <typename C>
    struct Coefficient
    {
        const C *values = nullptr;
        bool per_qp = false;
        const C &operator[](std::size_t qp) const { return per_qp ? values[qp] : values[0]; }
    };

    // ---- products -------------------------------------------------------------------------------------------

    inline Real value(Real c, Real u) { return c * u; }
    inline Vec value(Real c, const Vec &u) { return c * u; }
    inline Ten value(Real c, const Ten &u) { return c * u; }
    inline Vec value(const Ten &c, const Vec &u) { return c * u; }
    inline Vec value(const Vec &c, Real u) { return c * u; }
    inline Real value(const Vec &c, const Vec &u) { return c * u; }

    inline Real dot(Real t, Real x) { return t * x; }
    inline Real dot(const Vec &t, const Vec &x) { return t * x; }
    inline Real dot(const Ten &t, const Ten &x) { return t.contract(x); }

    inline Real dense(Real x) { return x; }
    inline Vec dense(const Vec &x) { return x; }
    inline Ten dense(const Ten &x) { return x; }

    // ---- what exists, at compile time -----------------------------------------------------------------------

    namespace detail
    {
        template <typename, typename = void>
        struct valid : std::false_type
        {
        };
        template <typename... A>
        struct valid<std::tuple<A...>, std::void_t<decltype(value(std::declval<A>()...))>> : std::true_type
        {
        };
        template <typename T, typename X, typename = void>
        struct dottable : std::false_type
        {
        };
        template <typename T, typename X>
        struct dottable<T, X, std::void_t<decltype(dot(std::declval<T>(), std::declval<X>()))>> : std::true_type
        {
        };
    } // namespace detail

    /** @brief value(C, U) exists. */
    template <typename C, typename U>
    constexpr bool has_value = detail::valid<std::tuple<C, U>>::value;

    /** @brief A block test^T * C * trial exists. */
    template <typename T, typename C, typename U>
    constexpr bool has_block = [] {
        if constexpr (has_value<C, U>)
            return detail::dottable<T, decltype(value(std::declval<C>(), std::declval<U>()))>::value;
        else
            return false;
    }();

    /** @brief A source test^T * C exists. */
    template <typename T, typename C>
    constexpr bool has_source = detail::dottable<T, C>::value;

    // ---- kernels ------------------------------------------------------------------------------------------

    /**
     * @brief M(i, j) = sum_qp w * dot(test_i, value(C, trial_j)); M is resized to (test dofs, trial dofs).
     * @tparam Symmetric The kernel is symmetric (a mass matrix: same table on both sides, scalar coefficient): only
     * the upper triangle is computed.
     */
    template <bool Symmetric = false, typename Test, typename C, typename Trial>
    void matrix(const QpDofTable<Test> &test, const Coefficient<C> &coefficient, const QpDofTable<Trial> &trial,
                const std::vector<Real> &weights, libMesh::DenseMatrix<Real> &out)
    {
        using Value = decltype(value(std::declval<C>(), std::declval<Trial>()));
        const std::size_t nq = weights.size();
        const std::size_t n_test = test.empty() ? 0 : test[0].size();
        const std::size_t n_trial = trial.empty() ? 0 : trial[0].size();
        out.resize(static_cast<unsigned int>(n_test), static_cast<unsigned int>(n_trial));
        std::vector<Value> cu(n_trial);
        if constexpr (Symmetric && std::is_same_v<Test, Trial> && std::is_same_v<C, Real>)
            if (&test == &trial)
            {
                for (std::size_t qp = 0; qp < nq; ++qp)
                {
                    const Real wc = weights[qp] * coefficient[qp];
                    const auto &t = test[qp];
                    for (std::size_t i = 0; i < n_test; ++i)
                        for (std::size_t j = i; j < n_trial; ++j)
                            out(static_cast<unsigned int>(i), static_cast<unsigned int>(j)) += wc * dot(t[i], t[j]);
                }
                for (std::size_t i = 0; i < n_test; ++i)
                    for (std::size_t j = 0; j < i; ++j)
                        out(static_cast<unsigned int>(i), static_cast<unsigned int>(j)) =
                            out(static_cast<unsigned int>(j), static_cast<unsigned int>(i));
                return;
            }
        for (std::size_t qp = 0; qp < nq; ++qp)
        {
            const Real w = weights[qp];
            const C &c = coefficient[qp];
            for (std::size_t j = 0; j < n_trial; ++j)
                cu[j] = value(c, trial[qp][j]);
            const auto &t = test[qp];
            for (std::size_t i = 0; i < n_test; ++i)
                for (std::size_t j = 0; j < n_trial; ++j)
                    out(static_cast<unsigned int>(i), static_cast<unsigned int>(j)) += w * dot(t[i], cu[j]);
        }
    }

    /** @brief F(i) = sum_qp w * dot(test_i, C); F is resized to the test dofs. */
    template <typename Test, typename C>
    void vector(const QpDofTable<Test> &test, const Coefficient<C> &coefficient, const std::vector<Real> &weights,
                libMesh::DenseVector<Real> &out)
    {
        const std::size_t nq = weights.size();
        out.resize(test.empty() ? 0 : static_cast<unsigned int>(test[0].size()));
        for (std::size_t qp = 0; qp < nq; ++qp)
        {
            const Real w = weights[qp];
            const C &c = coefficient[qp];
            const auto &t = test[qp];
            for (std::size_t i = 0; i < t.size(); ++i)
                out(static_cast<unsigned int>(i)) += w * dot(t[i], c);
        }
    }

    /**
     * @brief F(i) += sum_qp w * dot(test_i, sum_j value(C, trial_j) u_j): the matrix kernel applied to @p u, computed
     * at the quadrature points (O(nq (n_test + n_trial)) instead of building the matrix).
     */
    template <typename Test, typename C, typename Trial>
    void applied(const QpDofTable<Test> &test, const Coefficient<C> &coefficient, const QpDofTable<Trial> &trial,
                 const libMesh::DenseVector<Real> &u, const std::vector<Real> &weights, libMesh::DenseVector<Real> &out)
    {
        using Value = decltype(dense(value(std::declval<C>(), std::declval<Trial>())));
        for (std::size_t qp = 0; qp < weights.size(); ++qp)
        {
            const C &c = coefficient[qp];
            Value x{};
            const auto &trial_qp = trial[qp];
            for (std::size_t j = 0; j < trial_qp.size(); ++j)
                x += dense(value(c, trial_qp[j])) * u(static_cast<unsigned int>(j));
            const Real w = weights[qp];
            const auto &t = test[qp];
            for (std::size_t i = 0; i < t.size(); ++i)
                out(static_cast<unsigned int>(i)) += w * dot(t[i], x);
        }
    }

    /** @brief Table[qp][j] = value(C, trial_j): the coefficient times the trial at each point, no weight. */
    template <typename C, typename Trial, typename Out>
    void table(const Coefficient<C> &coefficient, const QpDofTable<Trial> &trial, QpDofTable<Out> &out)
    {
        out.assign(trial.size(), {});
        for (std::size_t qp = 0; qp < trial.size(); ++qp)
        {
            const C &c = coefficient[qp];
            out[qp].reserve(trial[qp].size());
            for (const Trial &u : trial[qp])
                out[qp].push_back(dense(value(c, u)));
        }
    }
} // namespace pytec3d::kernel

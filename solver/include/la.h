/**
 * @file la.h
 * @brief The small linear algebra of the solver core: scalars, 3-vectors, 3x3 tensors and dense element matrices.
 *
 * The kernels, the expressions of element_types.json (compiled by generate_registries.py), the laws and the element
 * evaluations compute with these types only; the finite element library (MFEM) and PETSc stay below them. They keep
 * the operations the element library uses:
 *
 *   VectorValue  v(i), +, -, scalar * v, v / scalar, v * w (dot), cross(), norm(), norm_sq(), unit()
 *   TensorValue  T(i, j), +, -, scalar * T, T * v, T * S, contract() (T:S), transpose(), inverse(), det(), tr()
 *   outer_product(v, w), pi
 *   DenseMatrix  m(), n(), resize(), zero(), (i, j), +=, add(), scale(), vector_mult_add(), cholesky_solve(),
 *                get_values() (row-major)
 *   DenseVector  size(), resize(), zero(), (i), +=, add(), scale(), get_values()
 */
#pragma once

#include <cmath>
#include <cstddef>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace pytec3d::la
{
    using Real = double;
    using Number = double;
    using dof_id_type = std::int64_t;
    using subdomain_id_type = std::int32_t; ///< A TYPE id.
    constexpr unsigned int invalid_uint = std::numeric_limits<unsigned int>::max();
    constexpr Real pi = 3.14159265358979323846264338327950288;

    // ---- 3-vectors -------------------------------------------------------------------------------------------

    template <typename T>
    class VectorValue
    {
    public:
        VectorValue() = default;
        VectorValue(T x, T y = 0, T z = 0) : v_{x, y, z} {}

        T &operator()(unsigned int i) { return v_[i]; }
        const T &operator()(unsigned int i) const { return v_[i]; }
        T &operator[](unsigned int i) { return v_[i]; }
        const T &operator[](unsigned int i) const { return v_[i]; }

        VectorValue &operator+=(const VectorValue &o)
        {
            for (int i = 0; i < 3; ++i)
                v_[i] += o.v_[i];
            return *this;
        }
        VectorValue &operator-=(const VectorValue &o)
        {
            for (int i = 0; i < 3; ++i)
                v_[i] -= o.v_[i];
            return *this;
        }
        VectorValue &operator*=(T s)
        {
            for (int i = 0; i < 3; ++i)
                v_[i] *= s;
            return *this;
        }
        VectorValue &operator/=(T s)
        {
            for (int i = 0; i < 3; ++i)
                v_[i] /= s;
            return *this;
        }
        VectorValue operator-() const { return {-v_[0], -v_[1], -v_[2]}; }
        VectorValue operator+(const VectorValue &o) const { return {v_[0] + o.v_[0], v_[1] + o.v_[1], v_[2] + o.v_[2]}; }
        VectorValue operator-(const VectorValue &o) const { return {v_[0] - o.v_[0], v_[1] - o.v_[1], v_[2] - o.v_[2]}; }
        VectorValue operator*(T s) const { return {v_[0] * s, v_[1] * s, v_[2] * s}; }
        VectorValue operator/(T s) const { return {v_[0] / s, v_[1] / s, v_[2] / s}; }
        /** @brief The dot product. */
        T operator*(const VectorValue &o) const { return v_[0] * o.v_[0] + v_[1] * o.v_[1] + v_[2] * o.v_[2]; }
        bool operator==(const VectorValue &o) const { return v_[0] == o.v_[0] && v_[1] == o.v_[1] && v_[2] == o.v_[2]; }
        bool operator!=(const VectorValue &o) const { return !(*this == o); }

        VectorValue cross(const VectorValue &o) const
        {
            return {v_[1] * o.v_[2] - v_[2] * o.v_[1], v_[2] * o.v_[0] - v_[0] * o.v_[2], v_[0] * o.v_[1] - v_[1] * o.v_[0]};
        }
        T norm_sq() const { return (*this) * (*this); }
        T norm() const { return std::sqrt(norm_sq()); }
        VectorValue unit() const { return *this / norm(); }
        void zero() { v_[0] = v_[1] = v_[2] = 0; }
        void add_scaled(const VectorValue &o, T s) { *this += o * s; }

    private:
        T v_[3] = {0, 0, 0};
    };

    template <typename T>
    VectorValue<T> operator*(T s, const VectorValue<T> &v) { return v * s; }
    template <typename T>
    VectorValue<T> operator*(int s, const VectorValue<T> &v) { return v * static_cast<T>(s); }

    // ---- 3x3 tensors -----------------------------------------------------------------------------------------

    template <typename T>
    class TensorValue
    {
    public:
        TensorValue() = default;
        TensorValue(T xx, T xy = 0, T xz = 0, T yx = 0, T yy = 0, T yz = 0, T zx = 0, T zy = 0, T zz = 0)
            : t_{xx, xy, xz, yx, yy, yz, zx, zy, zz}
        {
        }

        T &operator()(unsigned int i, unsigned int j) { return t_[3 * i + j]; }
        const T &operator()(unsigned int i, unsigned int j) const { return t_[3 * i + j]; }

        TensorValue &operator+=(const TensorValue &o)
        {
            for (int i = 0; i < 9; ++i)
                t_[i] += o.t_[i];
            return *this;
        }
        TensorValue &operator-=(const TensorValue &o)
        {
            for (int i = 0; i < 9; ++i)
                t_[i] -= o.t_[i];
            return *this;
        }
        TensorValue &operator*=(T s)
        {
            for (int i = 0; i < 9; ++i)
                t_[i] *= s;
            return *this;
        }
        TensorValue &operator/=(T s)
        {
            for (int i = 0; i < 9; ++i)
                t_[i] /= s;
            return *this;
        }
        TensorValue operator-() const
        {
            TensorValue r(*this);
            r *= T(-1);
            return r;
        }
        TensorValue operator+(const TensorValue &o) const
        {
            TensorValue r(*this);
            r += o;
            return r;
        }
        TensorValue operator-(const TensorValue &o) const
        {
            TensorValue r(*this);
            r -= o;
            return r;
        }
        TensorValue operator*(T s) const
        {
            TensorValue r(*this);
            r *= s;
            return r;
        }
        TensorValue operator/(T s) const
        {
            TensorValue r(*this);
            r /= s;
            return r;
        }
        VectorValue<T> operator*(const VectorValue<T> &v) const
        {
            return {t_[0] * v(0) + t_[1] * v(1) + t_[2] * v(2), t_[3] * v(0) + t_[4] * v(1) + t_[5] * v(2),
                    t_[6] * v(0) + t_[7] * v(1) + t_[8] * v(2)};
        }
        TensorValue operator*(const TensorValue &o) const
        {
            TensorValue r;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    for (unsigned int k = 0; k < 3; ++k)
                        r(i, j) += (*this)(i, k) * o(k, j);
            return r;
        }
        bool operator==(const TensorValue &o) const
        {
            for (int i = 0; i < 9; ++i)
                if (t_[i] != o.t_[i])
                    return false;
            return true;
        }

        /** @brief T : S, the double contraction. */
        T contract(const TensorValue &o) const
        {
            T s = 0;
            for (int i = 0; i < 9; ++i)
                s += t_[i] * o.t_[i];
            return s;
        }
        TensorValue transpose() const
        {
            return {t_[0], t_[3], t_[6], t_[1], t_[4], t_[7], t_[2], t_[5], t_[8]};
        }
        T det() const
        {
            return t_[0] * (t_[4] * t_[8] - t_[5] * t_[7]) - t_[1] * (t_[3] * t_[8] - t_[5] * t_[6]) +
                   t_[2] * (t_[3] * t_[7] - t_[4] * t_[6]);
        }
        T tr() const { return t_[0] + t_[4] + t_[8]; }
        TensorValue inverse() const
        {
            const T d = det();
            if (d == T(0))
                throw std::runtime_error("TensorValue::inverse(): singular tensor");
            return TensorValue(t_[4] * t_[8] - t_[5] * t_[7], t_[2] * t_[7] - t_[1] * t_[8], t_[1] * t_[5] - t_[2] * t_[4],
                               t_[5] * t_[6] - t_[3] * t_[8], t_[0] * t_[8] - t_[2] * t_[6], t_[2] * t_[3] - t_[0] * t_[5],
                               t_[3] * t_[7] - t_[4] * t_[6], t_[1] * t_[6] - t_[0] * t_[7], t_[0] * t_[4] - t_[1] * t_[3]) /
                   d;
        }
        T norm() const { return std::sqrt(contract(*this)); }
        void zero()
        {
            for (int i = 0; i < 9; ++i)
                t_[i] = 0;
        }

    private:
        T t_[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    };

    template <typename T>
    TensorValue<T> operator*(T s, const TensorValue<T> &t) { return t * s; }
    template <typename T>
    TensorValue<T> operator*(int s, const TensorValue<T> &t) { return t * static_cast<T>(s); }

    template <typename T>
    TensorValue<T> outer_product(const VectorValue<T> &a, const VectorValue<T> &b)
    {
        return {a(0) * b(0), a(0) * b(1), a(0) * b(2), a(1) * b(0), a(1) * b(1), a(1) * b(2),
                a(2) * b(0), a(2) * b(1), a(2) * b(2)};
    }

    using RealVectorValue = VectorValue<Real>;
    using RealTensorValue = TensorValue<Real>;
    using RealGradient = RealVectorValue;
    using Point = RealVectorValue;

    // ---- dense element vectors and matrices --------------------------------------------------------------------

    template <typename T>
    class DenseVector
    {
    public:
        DenseVector() = default;
        explicit DenseVector(unsigned int n) : v_(n, T(0)) {}

        unsigned int size() const { return static_cast<unsigned int>(v_.size()); }
        bool empty() const { return v_.empty(); }
        void resize(unsigned int n) { v_.assign(n, T(0)); }
        void zero() { std::fill(v_.begin(), v_.end(), T(0)); }
        T &operator()(unsigned int i) { return v_[i]; }
        const T &operator()(unsigned int i) const { return v_[i]; }
        DenseVector &operator+=(const DenseVector &o)
        {
            for (std::size_t i = 0; i < v_.size(); ++i)
                v_[i] += o.v_[i];
            return *this;
        }
        void add(T s, const DenseVector &o)
        {
            for (std::size_t i = 0; i < v_.size(); ++i)
                v_[i] += s * o.v_[i];
        }
        void scale(T s)
        {
            for (T &x : v_)
                x *= s;
        }
        T l2_norm() const
        {
            T s = 0;
            for (const T x : v_)
                s += x * x;
            return std::sqrt(s);
        }
        std::vector<T> &get_values() { return v_; }
        const std::vector<T> &get_values() const { return v_; }
        T *data() { return v_.data(); }
        const T *data() const { return v_.data(); }

    private:
        std::vector<T> v_;
    };

    /** @brief A dense row-major matrix. */
    template <typename T>
    class DenseMatrix
    {
    public:
        DenseMatrix() = default;
        DenseMatrix(unsigned int m, unsigned int n) { resize(m, n); }

        unsigned int m() const { return m_; }
        unsigned int n() const { return n_; }
        void resize(unsigned int m, unsigned int n)
        {
            m_ = m;
            n_ = n;
            v_.assign(static_cast<std::size_t>(m) * n, T(0));
            factored_ = false;
        }
        void zero()
        {
            std::fill(v_.begin(), v_.end(), T(0));
            factored_ = false;
        }
        T &operator()(unsigned int i, unsigned int j) { return v_[static_cast<std::size_t>(i) * n_ + j]; }
        const T &operator()(unsigned int i, unsigned int j) const { return v_[static_cast<std::size_t>(i) * n_ + j]; }
        DenseMatrix &operator+=(const DenseMatrix &o)
        {
            check_same(o);
            for (std::size_t i = 0; i < v_.size(); ++i)
                v_[i] += o.v_[i];
            return *this;
        }
        void add(T s, const DenseMatrix &o)
        {
            check_same(o);
            for (std::size_t i = 0; i < v_.size(); ++i)
                v_[i] += s * o.v_[i];
        }
        void scale(T s)
        {
            for (T &x : v_)
                x *= s;
        }
        /** @brief out += factor * this * in. */
        template <typename In>
        void vector_mult_add(DenseVector<T> &out, T factor, const In &in) const
        {
            for (unsigned int i = 0; i < m_; ++i)
            {
                T s = 0;
                const T *row = &v_[static_cast<std::size_t>(i) * n_];
                for (unsigned int j = 0; j < n_; ++j)
                    s += row[j] * value(in, j);
                out(i) += factor * s;
            }
        }
        /**
         * @brief Solves this * x = b for a symmetric positive definite matrix (Cholesky). The factorization is computed
         * at the first call and kept: the matrix must not change afterwards.
         */
        void cholesky_solve(const DenseVector<T> &b, DenseVector<T> &x)
        {
            if (m_ != n_)
                throw std::runtime_error("DenseMatrix::cholesky_solve(): not square");
            if (!factored_)
                factor();
            x.resize(m_);
            // L y = b, L^T x = y; L is stored in the lower triangle.
            for (unsigned int i = 0; i < m_; ++i)
            {
                T s = b(i);
                for (unsigned int k = 0; k < i; ++k)
                    s -= (*this)(i, k) * x(k);
                x(i) = s / (*this)(i, i);
            }
            for (unsigned int ii = m_; ii-- > 0;)
            {
                T s = x(ii);
                for (unsigned int k = ii + 1; k < m_; ++k)
                    s -= (*this)(k, ii) * x(k);
                x(ii) = s / (*this)(ii, ii);
            }
        }
        std::vector<T> &get_values() { return v_; }
        const std::vector<T> &get_values() const { return v_; }

    private:
        static T value(const DenseVector<T> &v, unsigned int j) { return v(j); }
        static T value(const std::vector<T> &v, unsigned int j) { return v[j]; }

        void check_same(const DenseMatrix &o) const
        {
            if (o.m_ != m_ || o.n_ != n_)
                throw std::runtime_error("DenseMatrix: sizes differ");
        }
        void factor()
        {
            for (unsigned int j = 0; j < n_; ++j)
            {
                T d = (*this)(j, j);
                for (unsigned int k = 0; k < j; ++k)
                    d -= (*this)(j, k) * (*this)(j, k);
                if (!(d > T(0)))
                    throw std::runtime_error("DenseMatrix::cholesky_solve(): matrix not positive definite");
                d = std::sqrt(d);
                (*this)(j, j) = d;
                for (unsigned int i = j + 1; i < m_; ++i)
                {
                    T s = (*this)(i, j);
                    for (unsigned int k = 0; k < j; ++k)
                        s -= (*this)(i, k) * (*this)(j, k);
                    (*this)(i, j) = s / d;
                }
            }
            factored_ = true;
        }

        unsigned int m_ = 0, n_ = 0;
        std::vector<T> v_;
        bool factored_ = false;
    };
} // namespace pytec3d::la

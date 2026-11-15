/**
 * @file curve.h
 * @brief A curve of the model data, as the nonlinear laws read it: y(x) or y(x, p) or y(x, p, q), given as points,
 *        as a table on a grid, or as an analytic law with its derivative.
 *
 * Three forms (model data "nonlinear_laws", the keys of the law's curve, nonlinear_laws.json):
 *
 *   - points: "<x>": [...], "<y>": [...]: piecewise linear; outside the points the curve goes on with its first and
 *     last segments, so the points must already describe the behaviour wanted there (e.g. end on a segment of slope
 *     1/mu0 for a BH curve that saturates to vacuum);
 *   - table (a curve with arguments, "args" in nonlinear_laws.json): "<x>": [...], "<arg>": [...] for each argument,
 *     "<y>": nested arrays [x][arg1][arg2]: multilinear on the grid, extended linearly outside as the points are;
 *   - analytic: "<y>": "<expression>" of the variables <x> (and the arguments), with "d<y>": "<expression>" its
 *     derivative in x (expression.h: + - * / ^, functions, pi, e, mu0).
 *
 * value(x, ...) is y, slope(x, ...) is dy/dx (the derivative the Newton tangent needs; the arguments are held).
 */
#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

#include <libmesh/libmesh_common.h>

#include "expression.h"

namespace pytec3d
{
    class CurveData
    {
    public:
        CurveData() = default;

        /** @brief Points: x strictly increasing, at least two (checked by the model data reader). */
        CurveData(std::vector<libMesh::Real> x, std::vector<libMesh::Real> y);

        /** @brief A table on the grid @p axes (axes[0] the x of the curve, the others its arguments; each strictly
         * increasing, at least two values): @p values in row-major order, x the slowest index. */
        static CurveData table(std::vector<std::vector<libMesh::Real>> axes, std::vector<libMesh::Real> values);

        /** @brief An analytic law: @p value y(x, args...) and @p derivative dy/dx, of the same variables. */
        static CurveData analytic(Expression value, Expression derivative);

        /** @brief The points (or the x axis of a table); empty for an analytic curve. */
        const std::vector<libMesh::Real> &x() const { return axes_.empty() ? empty_ : axes_[0]; }
        const std::vector<libMesh::Real> &y() const { return values_; }
        /** @brief The largest x given (points, table); +inf for an analytic curve (no extrapolation). */
        libMesh::Real x_max() const
        {
            return axes_.empty() ? std::numeric_limits<libMesh::Real>::infinity() : axes_[0].back();
        }
        /** @brief Number of arguments besides x (0 for points). */
        unsigned int n_args() const { return n_args_; }
        bool is_analytic() const { return analytic_ != nullptr; }

        libMesh::Real value(libMesh::Real x) const { return evaluate(&x, 1, false); }
        libMesh::Real slope(libMesh::Real x) const { return evaluate(&x, 1, true); }
        libMesh::Real value(libMesh::Real x, libMesh::Real p) const
        {
            const libMesh::Real a[2] = {x, p};
            return evaluate(a, 2, false);
        }
        libMesh::Real slope(libMesh::Real x, libMesh::Real p) const
        {
            const libMesh::Real a[2] = {x, p};
            return evaluate(a, 2, true);
        }
        libMesh::Real value(libMesh::Real x, libMesh::Real p, libMesh::Real q) const
        {
            const libMesh::Real a[3] = {x, p, q};
            return evaluate(a, 3, false);
        }
        libMesh::Real slope(libMesh::Real x, libMesh::Real p, libMesh::Real q) const
        {
            const libMesh::Real a[3] = {x, p, q};
            return evaluate(a, 3, true);
        }

    private:
        /** @brief y (or dy/dx) at the @p n coordinates @p a (x first); extra coordinates beyond the curve's are
         * an error, missing ones too. */
        libMesh::Real evaluate(const libMesh::Real *a, unsigned int n, bool derivative) const;

        /** @brief Index of the interval [v_i, v_i+1] of @p axis used at @p x, the first or last one outside. */
        static std::size_t segment(const std::vector<libMesh::Real> &axis, libMesh::Real x);

        std::vector<std::vector<libMesh::Real>> axes_; ///< Points / table: the x axis, then the arguments.
        std::vector<libMesh::Real> values_;            ///< y on the grid, row-major (x slowest).
        unsigned int n_args_ = 0;
        std::shared_ptr<const Expression> analytic_, derivative_;
        static const std::vector<libMesh::Real> empty_;
    };
} // namespace pytec3d

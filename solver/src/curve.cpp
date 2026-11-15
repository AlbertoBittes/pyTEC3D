/**
 * @file curve.cpp
 * @brief CurveData: points (piecewise linear), tables (multilinear) and analytic laws of the model data.
 */
#include "curve.h"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

namespace pytec3d
{
    const std::vector<pytec3d::la::Real> CurveData::empty_;

    CurveData::CurveData(std::vector<pytec3d::la::Real> x, std::vector<pytec3d::la::Real> y)
    {
        axes_.push_back(std::move(x));
        values_ = std::move(y);
    }

    CurveData CurveData::table(std::vector<std::vector<pytec3d::la::Real>> axes, std::vector<pytec3d::la::Real> values)
    {
        CurveData curve;
        curve.n_args_ = static_cast<unsigned int>(axes.size()) - 1;
        curve.axes_ = std::move(axes);
        curve.values_ = std::move(values);
        return curve;
    }

    CurveData CurveData::analytic(Expression value, Expression derivative)
    {
        CurveData curve;
        curve.n_args_ = static_cast<unsigned int>(value.variables().size()) - 1;
        curve.analytic_ = std::make_shared<const Expression>(std::move(value));
        curve.derivative_ = std::make_shared<const Expression>(std::move(derivative));
        return curve;
    }

    std::size_t CurveData::segment(const std::vector<pytec3d::la::Real> &axis, pytec3d::la::Real x)
    {
        if (x <= axis.front())
            return 0;
        if (x >= axis.back())
            return axis.size() - 2;
        const auto upper = std::upper_bound(axis.begin(), axis.end(), x);
        return static_cast<std::size_t>(std::distance(axis.begin(), upper) - 1);
    }

    pytec3d::la::Real CurveData::evaluate(const pytec3d::la::Real *a, unsigned int n, bool derivative) const
    {
        if (n != n_args_ + 1)
            throw std::logic_error("curve evaluated with " + std::to_string(n) + " coordinate(s); it has " +
                                   std::to_string(n_args_ + 1));
        if (analytic_)
        {
            double args[3] = {0.0, 0.0, 0.0};
            for (unsigned int k = 0; k < n; ++k)
                args[k] = static_cast<double>(a[k]);
            return (derivative ? *derivative_ : *analytic_).evaluate(args);
        }
        // Multilinear on the cell of the grid holding a (the end cells outside): sum over its 2^n corners.
        std::size_t lower[3] = {0, 0, 0};
        pytec3d::la::Real t[3] = {0.0, 0.0, 0.0}, h[3] = {1.0, 1.0, 1.0};
        for (unsigned int k = 0; k < n; ++k)
        {
            const std::vector<pytec3d::la::Real> &axis = axes_[k];
            lower[k] = segment(axis, a[k]);
            h[k] = axis[lower[k] + 1] - axis[lower[k]];
            t[k] = (a[k] - axis[lower[k]]) / h[k];
        }
        pytec3d::la::Real out = 0.0;
        for (unsigned int corner = 0; corner < (1u << n); ++corner)
        {
            std::size_t flat = 0;
            pytec3d::la::Real w = 1.0;
            for (unsigned int k = 0; k < n; ++k)
            {
                const unsigned int bit = (corner >> k) & 1u;
                flat = flat * axes_[k].size() + lower[k] + bit;
                if (k == 0 && derivative)
                    w *= (bit ? 1.0 : -1.0) / h[0]; // d/dx of the linear weight in x
                else
                    w *= bit ? t[k] : 1.0 - t[k];
            }
            out += w * values_[flat];
        }
        return out;
    }
} // namespace pytec3d

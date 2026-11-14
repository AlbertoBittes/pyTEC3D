/**
 * @file curve.cpp
 * @brief CurveData: piecewise-linear interpolation of a curve of the model data.
 */
#include "curve.h"

#include <algorithm>
#include <iterator>
#include <utility>

namespace pytec3d
{
    CurveData::CurveData(std::vector<pytec3d::la::Real> x, std::vector<pytec3d::la::Real> y)
        : x_(std::move(x)), y_(std::move(y))
    {
    }

    pytec3d::la::Real CurveData::value(pytec3d::la::Real x) const
    {
        const std::size_t i = segment(x);
        return y_[i] + (x - x_[i]) / (x_[i + 1] - x_[i]) * (y_[i + 1] - y_[i]);
    }

    pytec3d::la::Real CurveData::slope(pytec3d::la::Real x) const
    {
        const std::size_t i = segment(x);
        return (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    }

    std::size_t CurveData::segment(pytec3d::la::Real x) const
    {
        if (x <= x_.front())
            return 0;
        if (x >= x_.back())
            return x_.size() - 2;
        const auto upper = std::upper_bound(x_.begin(), x_.end(), x);
        return static_cast<std::size_t>(std::distance(x_.begin(), upper) - 1);
    }
} // namespace pytec3d

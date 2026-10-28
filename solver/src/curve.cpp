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
    CurveData::CurveData(std::vector<libMesh::Real> x, std::vector<libMesh::Real> y)
        : x_(std::move(x)), y_(std::move(y))
    {
    }

    libMesh::Real CurveData::value(libMesh::Real x) const
    {
        const std::size_t i = segment(x);
        return y_[i] + (x - x_[i]) / (x_[i + 1] - x_[i]) * (y_[i + 1] - y_[i]);
    }

    libMesh::Real CurveData::slope(libMesh::Real x) const
    {
        const std::size_t i = segment(x);
        return (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
    }

    std::size_t CurveData::segment(libMesh::Real x) const
    {
        if (x <= x_.front())
            return 0;
        if (x >= x_.back())
            return x_.size() - 2;
        const auto upper = std::upper_bound(x_.begin(), x_.end(), x);
        return static_cast<std::size_t>(std::distance(x_.begin(), upper) - 1);
    }
} // namespace pytec3d

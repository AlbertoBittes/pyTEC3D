/**
 * @file curve.h
 * @brief A piecewise-linear curve y(x) of the model data, as the nonlinear laws read it.
 *
 * The user gives the points (x, y), e.g. the B and H of a BH curve, and nothing else: outside the points the
 * curve goes on with its first and last segments, so the points must already describe the behaviour wanted
 * there (e.g. end on a segment of slope 1/mu0 for a BH curve that saturates to vacuum).
 */
#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include <libmesh/libmesh_common.h>

namespace pytec3d
{
    class CurveData
    {
    public:
        CurveData() = default;

        /** @brief The points of the curve; x increasing, at least two (checked by the model data reader). */
        CurveData(std::vector<libMesh::Real> x, std::vector<libMesh::Real> y) : x_(std::move(x)), y_(std::move(y)) {}

        const std::vector<libMesh::Real> &x() const { return x_; }
        const std::vector<libMesh::Real> &y() const { return y_; }

        /** @brief y(x), linear on each segment and along the first and last segments outside the points. */
        libMesh::Real value(libMesh::Real x) const
        {
            const std::size_t i = segment(x);
            return y_[i] + (x - x_[i]) / (x_[i + 1] - x_[i]) * (y_[i + 1] - y_[i]);
        }

        /** @brief dy/dx at x: the slope of the segment used at x. */
        libMesh::Real slope(libMesh::Real x) const
        {
            const std::size_t i = segment(x);
            return (y_[i + 1] - y_[i]) / (x_[i + 1] - x_[i]);
        }

    private:
        /** @brief Index of the segment [x_i, x_i+1] used at @p x, the first or last one outside the points. */
        std::size_t segment(libMesh::Real x) const
        {
            if (x <= x_.front())
                return 0;
            if (x >= x_.back())
                return x_.size() - 2;
            const auto upper = std::upper_bound(x_.begin(), x_.end(), x);
            return static_cast<std::size_t>(std::distance(x_.begin(), upper) - 1);
        }

        std::vector<libMesh::Real> x_;
        std::vector<libMesh::Real> y_;
    };
} // namespace pytec3d

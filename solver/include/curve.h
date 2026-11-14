/**
 * @file curve.h
 * @brief A piecewise-linear curve y(x) of the model data, as the nonlinear laws read it.
 *
 * The user gives the points (x, y), e.g. the B and H of a BH curve, and nothing else: outside the points the
 * curve goes on with its first and last segments, so the points must already describe the behaviour wanted
 * there (e.g. end on a segment of slope 1/mu0 for a BH curve that saturates to vacuum).
 */
#pragma once

#include "la.h"
#include <cstddef>
#include <vector>


namespace pytec3d
{
    class CurveData
    {
    public:
        CurveData() = default;

        /** @brief The points of the curve; x increasing, at least two (checked by the model data reader). */
        CurveData(std::vector<pytec3d::la::Real> x, std::vector<pytec3d::la::Real> y);

        const std::vector<pytec3d::la::Real> &x() const { return x_; }
        const std::vector<pytec3d::la::Real> &y() const { return y_; }

        /** @brief y(x), linear on each segment and along the first and last segments outside the points. */
        pytec3d::la::Real value(pytec3d::la::Real x) const;

        /** @brief dy/dx at x: the slope of the segment used at x. */
        pytec3d::la::Real slope(pytec3d::la::Real x) const;

    private:
        /** @brief Index of the segment [x_i, x_i+1] used at @p x, the first or last one outside the points. */
        std::size_t segment(pytec3d::la::Real x) const;

        std::vector<pytec3d::la::Real> x_;
        std::vector<pytec3d::la::Real> y_;
    };
} // namespace pytec3d

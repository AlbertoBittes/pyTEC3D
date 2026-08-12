#include <libmesh/libmesh.h>
#include "bh_curve.h"

#include <algorithm>

namespace
{
    constexpr libMesh::Real kZeroB = 1e-14;
    // Vacuum permeability: commercial BH-curve implementations (e.g. ANSYS Maxwell)
    // extrapolate beyond the last table point with at least this slope, so deeply
    // saturated iron asymptotically behaves like vacuum rather than staying "soft".
    const libMesh::Real kMu0 = 4.0 * libMesh::pi * 1e-7;
    const libMesh::Real kVacuumReluctivity = 1.0 / kMu0;
}

BHCurve::BHCurve(std::vector<libMesh::Real> h_values, std::vector<libMesh::Real> b_values)
    : h_values_(std::move(h_values)),
      b_values_(std::move(b_values))
{
}

std::size_t BHCurve::segment_index(const std::vector<libMesh::Real> &x_values, libMesh::Real x)
{
    if (x <= x_values.front())
        return 0;

    if (x >= x_values.back())
        return x_values.size() - 2;

    const auto upper = std::upper_bound(x_values.begin(), x_values.end(), x);
    return static_cast<std::size_t>(std::distance(x_values.begin(), upper) - 1);
}

libMesh::Real BHCurve::evaluate_reluctivity_from_b(libMesh::Real b_value) const
{
    if (b_value <= kZeroB)
        return slope_h_vs_b(0.0);

    const libMesh::Real h_value = evaluate_h_from_b(b_value);
    return h_value / b_value;
}

libMesh::Real BHCurve::evaluate_h_from_b(libMesh::Real b_value) const
{
    const std::size_t idx = segment_index(b_values_, b_value);
    const libMesh::Real b0 = b_values_[idx];
    const libMesh::Real b1 = b_values_[idx + 1];
    const libMesh::Real h0 = h_values_[idx];
    const libMesh::Real h1 = h_values_[idx + 1];

    if (b_value >= b_values_.back())
    {
        const libMesh::Real slope = std::max((h1 - h0) / (b1 - b0), kVacuumReluctivity);
        return h1 + slope * (b_value - b1);
    }

    const libMesh::Real t = (b_value - b0) / (b1 - b0);
    return h0 + t * (h1 - h0);
}

libMesh::Real BHCurve::slope_h_vs_b(libMesh::Real b_value) const
{
    const std::size_t idx = segment_index(b_values_, b_value);
    const libMesh::Real b0 = b_values_[idx];
    const libMesh::Real b1 = b_values_[idx + 1];
    const libMesh::Real h0 = h_values_[idx];
    const libMesh::Real h1 = h_values_[idx + 1];

    const libMesh::Real slope = (h1 - h0) / (b1 - b0);
    if (b_value >= b_values_.back())
        return std::max(slope, kVacuumReluctivity);
    return slope;
}

/**
 * @file bh_curve.h
 * @brief Monotonic B-H curve storage and interpolation helpers.
 */
#ifndef BH_CURVE_H
#define BH_CURVE_H

#include <vector>

#include "apec_m_visibility.h"
#include <libmesh/libmesh.h>

class APEC_M_HIDDEN BHCurve
{
public:
    BHCurve() = default;
    BHCurve(std::vector<libMesh::Real> h_values, std::vector<libMesh::Real> b_values);

    const std::vector<libMesh::Real> &h_values() const
    {
        return h_values_;
    }

    const std::vector<libMesh::Real> &b_values() const
    {
        return b_values_;
    }

    libMesh::Real evaluate_reluctivity_from_b(libMesh::Real b_value) const;
    libMesh::Real evaluate_h_from_b(libMesh::Real b_value) const;
    libMesh::Real slope_h_vs_b(libMesh::Real b_value) const;

private:
    std::vector<libMesh::Real> h_values_;
    std::vector<libMesh::Real> b_values_;

    static std::size_t segment_index(const std::vector<libMesh::Real> &x_values, libMesh::Real x);
};

#endif // BH_CURVE_H

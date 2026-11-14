/**
 * @file mesh_quality.h
 * @brief Shape quality of the mesh elements, by element shape (TRI3, HEX20, ...), checked once when the model is
 * loaded.
 *
 * Per shape: the largest edge length ratio, the smallest and largest interior angle, and the scaled Jacobian, the
 * smallest over the elements of min / max of det(J) over the quadrature points (1 for an undistorted element,
 * <= 0 for an inverted or degenerate one).
 */
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <cstdint>

namespace pytec3d::fem
{
    class FEMesh;
}

namespace pytec3d
{
    struct ShapeQuality
    {
        std::string shape;                          ///< Name of the cell type, e.g. "HEX20".
        std::size_t elements = 0;
        double max_edge_ratio = 0.0;
        double min_angle = 0.0;              ///< Degrees; 0 for 1D elements.
        double max_angle = 0.0;
        double min_scaled_jacobian = 1.0;
        std::int64_t worst_element = 0;     ///< Of the smallest scaled Jacobian (cell id).
        std::size_t distorted = 0;                  ///< Scaled Jacobian below kDistortedScaledJacobian.
        std::size_t degenerate = 0;                 ///< det(J) <= 0 at a quadrature point.
    };

    /** @brief Elements with a scaled Jacobian below this are reported as distorted. */
    inline constexpr double kDistortedScaledJacobian = 0.2;

    /** @brief The quality of the active elements, by shape. Collective. */
    std::vector<ShapeQuality> compute_mesh_quality(const fem::FEMesh &mesh);

    /** @brief Warns about degenerate and distorted elements. */
    void check_mesh_quality(const std::vector<ShapeQuality> &quality);
} // namespace pytec3d

/**
 * @file coordinate_rotations.h
 * @brief Local-to-global rotations of the coordinate systems, as callbacks provided by the caller.
 *
 * The coordinate systems are defined by the Python model (pytec3d.coordinate_system); python_bindings.cpp
 * implements these callbacks with it, so that the solver core does not depend on Python.
 */
#pragma once

#include <array>
#include <functional>
#include <vector>

#include <libmesh/point.h>
#include <libmesh/tensor_value.h>

struct CoordinateRotations
{
    using Rotations = std::vector<libMesh::RealTensorValue>;

    /** @brief Rotation of the coordinate system @p cs_id (COORDINATE_SYSTEM id) at each of @p points. */
    std::function<Rotations(unsigned int cs_id, const std::vector<libMesh::Point> &points)> at_points;

    /** @brief Rotation given by each triplet of angles in degrees (the "node_angles" point array). */
    std::function<Rotations(const std::vector<std::array<libMesh::Real, 3>> &angles)> from_angles;
};

#include "coordinate_system.h"

#include <cmath>

#include <libmesh/libmesh.h>

namespace
{
    ApecMTensorValue build_rotation_from_basis(const VectorValue<Real> &e0,
                                               const VectorValue<Real> &e1,
                                               const VectorValue<Real> &e2)
    {
        ApecMTensorValue rotation;
        rotation(0, 0) = e0(0);
        rotation(0, 1) = e0(1);
        rotation(0, 2) = e0(2);
        rotation(1, 0) = e1(0);
        rotation(1, 1) = e1(1);
        rotation(1, 2) = e1(2);
        rotation(2, 0) = e2(0);
        rotation(2, 1) = e2(1);
        rotation(2, 2) = e2(2);
        return rotation;
    }

    ApecMTensorValue build_local_to_global_from_axes(const VectorValue<Real> &x_axis,
                                                     const VectorValue<Real> &y_axis,
                                                     const VectorValue<Real> &z_axis)
    {
        ApecMTensorValue rotation;
        rotation(0, 0) = x_axis(0);
        rotation(1, 0) = x_axis(1);
        rotation(2, 0) = x_axis(2);
        rotation(0, 1) = y_axis(0);
        rotation(1, 1) = y_axis(1);
        rotation(2, 1) = y_axis(2);
        rotation(0, 2) = z_axis(0);
        rotation(1, 2) = z_axis(1);
        rotation(2, 2) = z_axis(2);
        return rotation;
    }

    VectorValue<Real> to_vector(const Point &point)
    {
        return VectorValue<Real>(point(0), point(1), point(2));
    }

    VectorValue<Real> normalize_or_error(const VectorValue<Real> &vec,
                                         const char *error_message)
    {
        const Real norm = vec.norm();
        if (norm <= 0.0)
            libmesh_error_msg(error_message);
        return vec / norm;
    }

    ApecMTensorValue build_local_to_global_from_points(const Point &origin,
                                                       const Point &nxax,
                                                       const Point &nxypl)
    {
        const VectorValue<Real> origin_vec = to_vector(origin);
        const VectorValue<Real> x_axis =
            normalize_or_error(to_vector(nxax) - origin_vec, "Coordinate system x-axis is invalid.");
        const VectorValue<Real> plane_vec =
            normalize_or_error(to_vector(nxypl) - origin_vec, "Coordinate system xy-plane vector is invalid.");
        const VectorValue<Real> z_axis =
            normalize_or_error(x_axis.cross(plane_vec), "Coordinate system axes are collinear.");
        const VectorValue<Real> y_axis =
            normalize_or_error(z_axis.cross(x_axis), "Coordinate system y-axis is invalid.");
        return build_local_to_global_from_axes(x_axis, y_axis, z_axis);
    }
} // namespace

CoordinateSystem::CoordinateSystem(Type type,
                                   const Point &origin,
                                   const Point &nxax,
                                   const Point &nxypl)
    : type_(type),
      origin_(origin),
      local_to_global_cart_(build_local_to_global_from_points(origin, nxax, nxypl))
{
}

CoordinateSystem::CoordinateSystem(Type type,
                                   const Point &origin,
                                   Real thxy_deg,
                                   Real thyz_deg,
                                   Real thzx_deg)
    : type_(type),
      origin_(origin),
      local_to_global_cart_(ApecMTensorValue::build_ansys_rotation_matrix(thxy_deg, thyz_deg, thzx_deg))
{
}

Point CoordinateSystem::axis_direction() const
{
    return Point(local_to_global_cart_(0, 2),
                 local_to_global_cart_(1, 2),
                 local_to_global_cart_(2, 2));
}

VectorValue<Real> CoordinateSystem::to_local_vector(const Point &global_point) const
{
    const VectorValue<Real> rel = to_vector(global_point - origin_);
    return local_to_global_cart_.transpose() * rel;
}

CartesianCoordinateSystem::CartesianCoordinateSystem(const Point &origin,
                                                     const Point &nxax,
                                                     const Point &nxypl)
    : CoordinateSystem(Type::Cartesian, origin, nxax, nxypl)
{
}

CartesianCoordinateSystem::CartesianCoordinateSystem(const Point &origin,
                                                     Real thxy_deg,
                                                     Real thyz_deg,
                                                     Real thzx_deg)
    : CoordinateSystem(Type::Cartesian, origin, thxy_deg, thyz_deg, thzx_deg)
{
}

ApecMTensorValue CartesianCoordinateSystem::global_to_local_rotation(const Point &global_point) const
{
    static_cast<void>(global_point);
    return local_to_global_cartesian().transpose();
}

Point CartesianCoordinateSystem::to_local_coordinates(const Point &global_point) const
{
    const VectorValue<Real> local = to_local_vector(global_point);
    return Point(local(0), local(1), local(2));
}

CylindricalCoordinateSystem::CylindricalCoordinateSystem(const Point &origin,
                                                         const Point &nxax,
                                                         const Point &nxypl)
    : CoordinateSystem(Type::Cylindrical, origin, nxax, nxypl)
{
}

CylindricalCoordinateSystem::CylindricalCoordinateSystem(const Point &origin,
                                                         Real thxy_deg,
                                                         Real thyz_deg,
                                                         Real thzx_deg)
    : CoordinateSystem(Type::Cylindrical, origin, thxy_deg, thyz_deg, thzx_deg)
{
}

ApecMTensorValue CylindricalCoordinateSystem::global_to_local_rotation(const Point &global_point) const
{
    const VectorValue<Real> local = to_local_vector(global_point);
    const Real x = local(0);
    const Real z = local(2);
    const Real r = std::sqrt(x * x + z * z);

    VectorValue<Real> e_r_local(1.0, 0.0, 0.0);
    VectorValue<Real> e_theta_local(0.0, 0.0, -1.0);
    if (r > 0.0)
    {
        const Real inv_r = 1.0 / r;
        e_r_local = VectorValue<Real>(x * inv_r, 0.0, z * inv_r);
        e_theta_local = VectorValue<Real>(z * inv_r, 0.0, -x * inv_r);
    }
    const VectorValue<Real> e_y_local(0.0, 1.0, 0.0);
    const VectorValue<Real> e_z_local = -e_theta_local;

    const auto &local_to_global = local_to_global_cartesian();
    const VectorValue<Real> e_r_global = local_to_global * e_r_local;
    const VectorValue<Real> e_y_global = local_to_global * e_y_local;
    const VectorValue<Real> e_z_global = local_to_global * e_z_local;

    return build_rotation_from_basis(e_r_global, e_y_global, e_z_global);
}

Point CylindricalCoordinateSystem::to_local_coordinates(const Point &global_point) const
{
    const VectorValue<Real> local = to_local_vector(global_point);
    const Real x = local(0);
    const Real z = local(2);
    const Real r = std::sqrt(x * x + z * z);
    const Real theta = (r > 0.0) ? std::atan2(-z, x) : 0.0;
    return Point(r, theta, local(1));
}

SphericalCoordinateSystem::SphericalCoordinateSystem(const Point &origin,
                                                     const Point &nxax,
                                                     const Point &nxypl)
    : CoordinateSystem(Type::Spherical, origin, nxax, nxypl)
{
}

SphericalCoordinateSystem::SphericalCoordinateSystem(const Point &origin,
                                                     Real thxy_deg,
                                                     Real thyz_deg,
                                                     Real thzx_deg)
    : CoordinateSystem(Type::Spherical, origin, thxy_deg, thyz_deg, thzx_deg)
{
}

ApecMTensorValue SphericalCoordinateSystem::global_to_local_rotation(const Point &global_point) const
{
    const VectorValue<Real> local = to_local_vector(global_point);
    const Real x = local(0);
    const Real y = local(1);
    const Real z = local(2);
    const Real r = std::sqrt(x * x + y * y + z * z);
    const Real rho = std::sqrt(x * x + y * y);

    VectorValue<Real> e_r_local(0.0, 0.0, 1.0);
    VectorValue<Real> e_theta_local(1.0, 0.0, 0.0);
    VectorValue<Real> e_phi_local(0.0, 1.0, 0.0);

    if (r > 0.0)
    {
        const Real inv_r = 1.0 / r;
        e_r_local = VectorValue<Real>(x * inv_r, y * inv_r, z * inv_r);

        if (rho > 0.0)
        {
            const Real inv_rho = 1.0 / rho;
            const Real cos_theta = z * inv_r;
            const Real sin_theta = rho * inv_r;
            const Real cos_phi = x * inv_rho;
            const Real sin_phi = y * inv_rho;

            e_theta_local = VectorValue<Real>(cos_theta * cos_phi,
                                              cos_theta * sin_phi,
                                              -sin_theta);
            e_phi_local = VectorValue<Real>(-sin_phi, cos_phi, 0.0);
        }
    }

    const auto &local_to_global = local_to_global_cartesian();
    const VectorValue<Real> e_r_global = local_to_global * e_r_local;
    const VectorValue<Real> e_theta_global = local_to_global * e_theta_local;
    const VectorValue<Real> e_phi_global = local_to_global * e_phi_local;

    return build_rotation_from_basis(e_r_global, e_theta_global, e_phi_global);
}

Point SphericalCoordinateSystem::to_local_coordinates(const Point &global_point) const
{
    const VectorValue<Real> local = to_local_vector(global_point);
    const Real x = local(0);
    const Real y = local(1);
    const Real z = local(2);
    const Real r = std::sqrt(x * x + y * y + z * z);
    if (r <= 0.0)
        return Point(0.0, 0.0, 0.0);
    const Real theta = std::acos(z / r);
    const Real phi = std::atan2(y, x);
    return Point(r, theta, phi);
}

std::unique_ptr<CoordinateSystem> make_coordinate_system_local(int kcs,
                                                               const Point &origin,
                                                               Real thxy_deg,
                                                               Real thyz_deg,
                                                               Real thzx_deg)
{
    switch (kcs)
    {
    case 0:
        return std::make_unique<CartesianCoordinateSystem>(origin, thxy_deg, thyz_deg, thzx_deg);
    case 1:
        return std::make_unique<CylindricalCoordinateSystem>(origin, thxy_deg, thyz_deg, thzx_deg);
    case 2:
        return std::make_unique<SphericalCoordinateSystem>(origin, thxy_deg, thyz_deg, thzx_deg);
    default:
        libmesh_error_msg("Unsupported coordinate system type.");
    }
}

std::unique_ptr<CoordinateSystem> make_coordinate_system_cs(int kcs,
                                                            const Point &norig,
                                                            const Point &nxax,
                                                            const Point &nxypl)
{
    switch (kcs)
    {
    case 0:
        return std::make_unique<CartesianCoordinateSystem>(norig, nxax, nxypl);
    case 1:
        return std::make_unique<CylindricalCoordinateSystem>(norig, nxax, nxypl);
    case 2:
        return std::make_unique<SphericalCoordinateSystem>(norig, nxax, nxypl);
    default:
        libmesh_error_msg("Unsupported coordinate system type.");
    }
}

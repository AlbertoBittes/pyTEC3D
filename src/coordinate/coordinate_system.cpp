#include "apec_m_coordinate_system/coordinate_system.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace apec_m
{
    namespace
    {
        constexpr double kPi = 3.141592653589793238462643383279502884;
        constexpr double kOrthonormalTol = 1e-8;
        constexpr double kFiniteTolerance = 1e-12;
        constexpr double kDefaultOrigin[3] = {0.0, 0.0, 0.0};
        // THXY, THYZ, THZX (deg) rotating the reference (Z-aligned) axis onto each global axis.
        constexpr double kXAxisAngles[3] = {0.0, 0.0, 90.0};
        constexpr double kYAxisAngles[3] = {0.0, -90.0, 0.0};
        constexpr double kZAxisAngles[3] = {0.0, 0.0, 0.0};

        inline void mat_t_vec3(const double *m, const double *v, double *out) noexcept
        {
            // Use temporaries so calls with out == v remain correct.
            const double x = m[0] * v[0] + m[3] * v[1] + m[6] * v[2];
            const double y = m[1] * v[0] + m[4] * v[1] + m[7] * v[2];
            const double z = m[2] * v[0] + m[5] * v[1] + m[8] * v[2];
            out[0] = x;
            out[1] = y;
            out[2] = z;
        }

        inline void mat_vec3(const double *m, const double *v, double *out) noexcept
        {
            // Use temporaries so calls with out == v remain correct.
            const double x = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
            const double y = m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
            const double z = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
            out[0] = x;
            out[1] = y;
            out[2] = z;
        }

        inline double snap_trig(double v) noexcept
        {
            if (std::abs(v) < kFiniteTolerance)
                return 0.0;
            if (std::abs(v - 1.0) < kFiniteTolerance)
                return 1.0;
            if (std::abs(v + 1.0) < kFiniteTolerance)
                return -1.0;
            return v;
        }

    } // namespace

    CoordinateSystem::CoordinateSystem(const double *origin,
                                       const double *local_to_global_cart,
                                       std::string cs_type)
        : _type(std::move(cs_type))
    {
        const double default_origin[3] = {0.0, 0.0, 0.0};
        const double default_rot[9] = {
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0};

        const double *origin_ptr = origin ? origin : default_origin;
        const double *rot_ptr = local_to_global_cart ? local_to_global_cart : default_rot;

        if (_type.empty())
            throw std::invalid_argument("cs_type must be a non-empty string.");

        _origin[0] = origin_ptr[0];
        _origin[1] = origin_ptr[1];
        _origin[2] = origin_ptr[2];
        for (int i = 0; i < 9; ++i)
            _local_to_global_cart[i] = rot_ptr[i];
    }

    const std::string &CoordinateSystem::type() const noexcept
    {
        return _type;
    }

    void CoordinateSystem::origin_raw(double *origin_out) const noexcept
    {
        origin_out[0] = _origin[0];
        origin_out[1] = _origin[1];
        origin_out[2] = _origin[2];
    }

    void CoordinateSystem::local_to_global_cart_raw(double *matrix_out) const noexcept
    {
        for (int i = 0; i < 9; ++i)
            matrix_out[i] = _local_to_global_cart[i];
    }

    void CoordinateSystem::axis_direction_raw(double *axis_out) const noexcept
    {
        axis_out[0] = _local_to_global_cart[2];
        axis_out[1] = _local_to_global_cart[5];
        axis_out[2] = _local_to_global_cart[8];
    }

    void CoordinateSystem::to_local_raw(
        const double *global_points,
        const std::size_t n_points,
        double *local_points_out) const
    {
        for (std::size_t i = 0; i < n_points; ++i)
        {
            const std::size_t base = i * 3;
            const double delta[3] = {
                global_points[base + 0] - _origin[0],
                global_points[base + 1] - _origin[1],
                global_points[base + 2] - _origin[2]};
            mat_t_vec3(_local_to_global_cart, delta, &local_points_out[base]);
        }
    }

    void CoordinateSystem::to_global_raw(
        const double *local_points,
        const std::size_t n_points,
        double *global_points_out) const
    {
        for (std::size_t i = 0; i < n_points; ++i)
        {
            const std::size_t base = i * 3;
            double global[3] = {0.0, 0.0, 0.0};
            mat_vec3(_local_to_global_cart, &local_points[base], global);
            global_points_out[base + 0] = _origin[0] + global[0];
            global_points_out[base + 1] = _origin[1] + global[1];
            global_points_out[base + 2] = _origin[2] + global[2];
        }
    }

    void CoordinateSystem::vector_components_to_local_raw(
        const double *global_vectors,
        const double *global_positions,
        const std::size_t n_vectors,
        double *local_components_out) const
    {
        for (std::size_t i = 0; i < n_vectors; ++i)
        {
            const std::size_t base = i * 3;
            mat_t_vec3(_local_to_global_cart, &global_vectors[base], &local_components_out[base]);
        }
    }

    void CoordinateSystem::vector_components_to_global_raw(
        const double *local_components,
        const double *global_positions,
        const std::size_t n_vectors,
        double *global_vectors_out) const
    {
        for (std::size_t i = 0; i < n_vectors; ++i)
        {
            const std::size_t base = i * 3;
            mat_vec3(_local_to_global_cart, &local_components[base], &global_vectors_out[base]);
        }
    }

    void CoordinateSystem::local_to_global_matrix_raw(
        const double *global_position,
        double *matrix_out) const
    {
        for (std::size_t a = 0; a < 3; ++a)
        {
            const double local_basis[3] = {
                (a == 0) ? 1.0 : 0.0,
                (a == 1) ? 1.0 : 0.0,
                (a == 2) ? 1.0 : 0.0};
            double global_vec[3] = {0.0, 0.0, 0.0};
            this->vector_components_to_global_raw(local_basis, global_position, 1, global_vec);
            for (std::size_t i = 0; i < 3; ++i)
                matrix_out[i * 3 + a] = global_vec[i];
        }
    }

    void CoordinateSystem::rotate_vector_to_global_raw(
        const double *local_components,
        const double *global_position,
        double *global_components_out) const
    {
        this->vector_components_to_global_raw(local_components, global_position, 1, global_components_out);
    }

    void CoordinateSystem::rotate_tensor_to_global_raw(
        const double *local_tensor,
        const double *global_position,
        double *global_tensor_out) const
    {
        double rotation[9] = {0.0};
        this->local_to_global_matrix_raw(global_position, rotation);

        double tmp[9] = {0.0};
        for (std::size_t i = 0; i < 3; ++i)
        {
            for (std::size_t j = 0; j < 3; ++j)
            {
                double sum = 0.0;
                for (std::size_t k = 0; k < 3; ++k)
                    sum += rotation[i * 3 + k] * local_tensor[k * 3 + j];
                tmp[i * 3 + j] = sum;
            }
        }

        for (std::size_t i = 0; i < 3; ++i)
        {
            for (std::size_t j = 0; j < 3; ++j)
            {
                double sum = 0.0;
                for (std::size_t k = 0; k < 3; ++k)
                    sum += tmp[i * 3 + k] * rotation[j * 3 + k];
                global_tensor_out[i * 3 + j] = sum;
            }
        }
    }

    void CoordinateSystem::build_local_to_global_matrix(
        const double angle_xy_deg,
        const double angle_yz_deg,
        const double angle_xz_deg,
        double *matrix_out)
    {
        const double angles[3] = {angle_xy_deg, angle_yz_deg, angle_xz_deg};
        build_local_to_global_coefficients_raw(angles, 1, matrix_out);
    }

    void CoordinateSystem::build_local_to_global_coefficients_raw(
        const double *angles_deg,
        const std::size_t n_rows,
        double *coeffs_out)
    {
        for (std::size_t i = 0; i < n_rows; ++i)
        {
            const std::size_t in = i * 3;
            const std::size_t out = i * 9;

            const double a1 = deg_to_rad(angles_deg[in + 0]);
            const double a2 = deg_to_rad(angles_deg[in + 1]);
            const double a3 = deg_to_rad(angles_deg[in + 2]);

            const double c1 = snap_trig(std::cos(a1));
            const double s1 = snap_trig(std::sin(a1));
            const double c2 = snap_trig(std::cos(a2));
            const double s2 = snap_trig(std::sin(a2));
            const double c3 = snap_trig(std::cos(a3));
            const double s3 = snap_trig(std::sin(a3));

            coeffs_out[out + 0] = c1 * c3 - s1 * s2 * s3;
            coeffs_out[out + 1] = -s1 * c2;
            coeffs_out[out + 2] = c1 * s3 + s1 * s2 * c3;
            coeffs_out[out + 3] = s1 * c3 + c1 * s2 * s3;
            coeffs_out[out + 4] = c1 * c2;
            coeffs_out[out + 5] = s1 * s3 - c1 * s2 * c3;
            coeffs_out[out + 6] = -c2 * s3;
            coeffs_out[out + 7] = s2;
            coeffs_out[out + 8] = c2 * c3;
        }
    }

    void CoordinateSystem::build_global_to_local_matrix(
        const double angle_xy_deg,
        const double angle_yz_deg,
        const double angle_xz_deg,
        double *matrix_out)
    {
        double local_to_global[9] = {0.0};
        build_local_to_global_matrix(angle_xy_deg, angle_yz_deg, angle_xz_deg, local_to_global);
        matrix_out[0] = local_to_global[0];
        matrix_out[1] = local_to_global[3];
        matrix_out[2] = local_to_global[6];
        matrix_out[3] = local_to_global[1];
        matrix_out[4] = local_to_global[4];
        matrix_out[5] = local_to_global[7];
        matrix_out[6] = local_to_global[2];
        matrix_out[7] = local_to_global[5];
        matrix_out[8] = local_to_global[8];
    }

    void CoordinateSystem::build_global_to_local(
        const double angle_xy_deg,
        const double angle_yz_deg,
        const double angle_xz_deg,
        double *matrix_out)
    {
        build_global_to_local_matrix(angle_xy_deg, angle_yz_deg, angle_xz_deg, matrix_out);
    }

    double CoordinateSystem::deg_to_rad(const double value) noexcept
    {
        return value * (kPi / 180.0);
    }

    double CoordinateSystem::rad_to_deg(const double value) noexcept
    {
        return value * (180.0 / kPi);
    }

    void CoordinateSystem::normalize3_or_error(
        const double *vec,
        double *out,
        const std::string &error_message)
    {
        const double n = std::sqrt(dot3(vec, vec));
        if (n <= 0.0)
            throw std::invalid_argument(error_message);
        out[0] = vec[0] / n;
        out[1] = vec[1] / n;
        out[2] = vec[2] / n;
    }

    void CoordinateSystem::build_local_to_global_from_points_raw(
        const double *origin,
        const double *nxax,
        const double *nxypl,
        double *matrix_out)
    {
        const double x_raw[3] = {
            nxax[0] - origin[0],
            nxax[1] - origin[1],
            nxax[2] - origin[2]};
        const double plane_raw[3] = {
            nxypl[0] - origin[0],
            nxypl[1] - origin[1],
            nxypl[2] - origin[2]};

        double x_axis[3] = {0.0, 0.0, 0.0};
        double plane_vec[3] = {0.0, 0.0, 0.0};
        normalize3_or_error(x_raw, x_axis, "Coordinate system x-axis is invalid.");
        normalize3_or_error(plane_raw, plane_vec, "Coordinate system xy-plane vector is invalid.");

        double z_cross[3] = {0.0, 0.0, 0.0};
        cross3(x_axis, plane_vec, z_cross);
        double z_axis[3] = {0.0, 0.0, 0.0};
        normalize3_or_error(z_cross, z_axis, "Coordinate system axes are collinear.");

        double y_cross[3] = {0.0, 0.0, 0.0};
        cross3(z_axis, x_axis, y_cross);
        double y_axis[3] = {0.0, 0.0, 0.0};
        normalize3_or_error(y_cross, y_axis, "Coordinate system y-axis is invalid.");

        matrix_out[0] = x_axis[0];
        matrix_out[1] = y_axis[0];
        matrix_out[2] = z_axis[0];
        matrix_out[3] = x_axis[1];
        matrix_out[4] = y_axis[1];
        matrix_out[5] = z_axis[1];
        matrix_out[6] = x_axis[2];
        matrix_out[7] = y_axis[2];
        matrix_out[8] = z_axis[2];
    }

    double CoordinateSystem::dot3(const double *a, const double *b) noexcept
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    void CoordinateSystem::cross3(const double *a, const double *b, double *out) noexcept
    {
        out[0] = a[1] * b[2] - a[2] * b[1];
        out[1] = a[2] * b[0] - a[0] * b[2];
        out[2] = a[0] * b[1] - a[1] * b[0];
    }

    CartesianCoordinateSystem::CartesianCoordinateSystem(
        const double *origin,
        const double *nxax,
        const double *nxypl,
        const double *angles_deg,
        const double *local_to_global_cart)
        : CoordinateSystem(
              origin ? origin : kDefaultOrigin,
              nullptr,
              "cartesian")
    {
        const double default_origin[3] = {0.0, 0.0, 0.0};
        const double default_nxax[3] = {1.0, 0.0, 0.0};
        const double default_nxypl[3] = {0.0, 1.0, 0.0};

        const double *o = origin ? origin : default_origin;
        _origin[0] = o[0];
        _origin[1] = o[1];
        _origin[2] = o[2];

        if (local_to_global_cart)
        {
            for (int i = 0; i < 9; ++i)
                _local_to_global_cart[i] = local_to_global_cart[i];
        }
        else if (angles_deg)
        {
            build_local_to_global_matrix(angles_deg[0], angles_deg[1], angles_deg[2], _local_to_global_cart);
        }
        else
        {
            const double *x = nxax ? nxax : default_nxax;
            const double *p = nxypl ? nxypl : default_nxypl;
            build_local_to_global_from_points_raw(_origin, x, p, _local_to_global_cart);
        }
    }

    CartesianCoordinateSystem CartesianCoordinateSystem::from_angles(
        const double *origin,
        const double *angles_deg)
    {
        return CartesianCoordinateSystem(origin, nullptr, nullptr, angles_deg);
    }

    CartesianCoordinateSystem CartesianCoordinateSystem::from_points(
        const double *origin,
        const double *nxax,
        const double *nxypl)
    {
        return CartesianCoordinateSystem(origin, nxax, nxypl, nullptr);
    }

    CylindricalCoordinateSystem::CylindricalCoordinateSystem(
        const double *origin,
        const double *nxax,
        const double *nxypl,
        const double *angles_deg,
        const double *local_to_global_cart)
        : CoordinateSystem(
              origin ? origin : kDefaultOrigin,
              nullptr,
              "cylindrical")
    {
        const double default_origin[3] = {0.0, 0.0, 0.0};
        const double default_nxax[3] = {1.0, 0.0, 0.0};
        const double default_nxypl[3] = {0.0, 0.0, -1.0};

        const double *o = origin ? origin : default_origin;
        _origin[0] = o[0];
        _origin[1] = o[1];
        _origin[2] = o[2];

        if (local_to_global_cart)
        {
            for (int i = 0; i < 9; ++i)
                _local_to_global_cart[i] = local_to_global_cart[i];
        }
        else if (angles_deg)
        {
            build_local_to_global_matrix(angles_deg[0], angles_deg[1], angles_deg[2], _local_to_global_cart);
        }
        else
        {
            const double *x = nxax ? nxax : default_nxax;
            const double *p = nxypl ? nxypl : default_nxypl;
            build_local_to_global_from_points_raw(_origin, x, p, _local_to_global_cart);
        }
    }

    CylindricalCoordinateSystem CylindricalCoordinateSystem::from_angles(
        const double *origin,
        const double *angles_deg)
    {
        return CylindricalCoordinateSystem(origin, nullptr, nullptr, angles_deg);
    }

    CylindricalCoordinateSystem CylindricalCoordinateSystem::from_points(
        const double *origin,
        const double *nxax,
        const double *nxypl)
    {
        return CylindricalCoordinateSystem(origin, nxax, nxypl, nullptr);
    }

    XAxisCylindricalCoordinateSystem::XAxisCylindricalCoordinateSystem(
        const double *origin,
        const double *local_to_global_cart)
        : CylindricalCoordinateSystem(origin, nullptr, nullptr, local_to_global_cart ? nullptr : kXAxisAngles, local_to_global_cart)
    {
    }

    YAxisCylindricalCoordinateSystem::YAxisCylindricalCoordinateSystem(
        const double *origin,
        const double *local_to_global_cart)
        : CylindricalCoordinateSystem(origin, nullptr, nullptr, local_to_global_cart ? nullptr : kYAxisAngles, local_to_global_cart)
    {
    }

    ZAxisCylindricalCoordinateSystem::ZAxisCylindricalCoordinateSystem(
        const double *origin,
        const double *local_to_global_cart)
        : CylindricalCoordinateSystem(origin, nullptr, nullptr, local_to_global_cart ? nullptr : kZAxisAngles, local_to_global_cart)
    {
    }

    void CylindricalCoordinateSystem::vector_components_to_local_raw(
        const double *global_vectors,
        const double *global_positions,
        const std::size_t n_vectors,
        double *local_components_out) const
    {
        for (std::size_t i = 0; i < n_vectors; ++i)
        {
            const std::size_t b = i * 3;
            double lp[3] = {
                global_positions[b + 0] - _origin[0],
                global_positions[b + 1] - _origin[1],
                global_positions[b + 2] - _origin[2]};
            mat_t_vec3(_local_to_global_cart, lp, lp);

            double lv[3] = {global_vectors[b + 0], global_vectors[b + 1], global_vectors[b + 2]};
            mat_t_vec3(_local_to_global_cart, lv, lv);

            const double x = lp[0];
            const double y = lp[1];
            const double r = std::sqrt(x * x + y * y);
            double comp_r = lv[0];
            double comp_theta = lv[1];
            if (r > 0.0)
            {
                const double inv_r = 1.0 / r;
                comp_r = (x * lv[0] + y * lv[1]) * inv_r;
                comp_theta = (-y * lv[0] + x * lv[1]) * inv_r;
            }

            local_components_out[b + 0] = comp_r;
            local_components_out[b + 1] = comp_theta;
            local_components_out[b + 2] = lv[2];
        }
    }

    void CylindricalCoordinateSystem::vector_components_to_global_raw(
        const double *local_components,
        const double *global_positions,
        const std::size_t n_vectors,
        double *global_vectors_out) const
    {
        for (std::size_t i = 0; i < n_vectors; ++i)
        {
            const std::size_t b = i * 3;
            double lp[3] = {
                global_positions[b + 0] - _origin[0],
                global_positions[b + 1] - _origin[1],
                global_positions[b + 2] - _origin[2]};
            mat_t_vec3(_local_to_global_cart, lp, lp);

            const double x = lp[0];
            const double y = lp[1];
            const double r = std::sqrt(x * x + y * y);

            double lv[3] = {local_components[b + 0], local_components[b + 1], local_components[b + 2]};
            if (r > 0.0)
            {
                const double inv_r = 1.0 / r;
                const double xr = x * inv_r;
                const double yr = y * inv_r;
                const double comp_r = lv[0];
                const double comp_theta = lv[1];
                lv[0] = xr * comp_r - yr * comp_theta;
                lv[1] = yr * comp_r + xr * comp_theta;
            }

            mat_vec3(_local_to_global_cart, lv, &global_vectors_out[b]);
        }
    }

    void CylindricalCoordinateSystem::to_local_raw(
        const double *global_points,
        const std::size_t n_points,
        double *local_points_out) const
    {
        CoordinateSystem::to_local_raw(global_points, n_points, local_points_out);
        for (std::size_t i = 0; i < n_points; ++i)
        {
            const std::size_t b = i * 3;
            const double x = local_points_out[b + 0];
            const double y = local_points_out[b + 1];
            const double z = local_points_out[b + 2];
            const double r = std::sqrt(x * x + y * y);
            local_points_out[b + 0] = r;
            local_points_out[b + 1] = (r > 0.0) ? rad_to_deg(std::atan2(y, x)) : 0.0;
            local_points_out[b + 2] = z;
        }
    }

    void CylindricalCoordinateSystem::to_global_raw(
        const double *local_points,
        const std::size_t n_points,
        double *global_points_out) const
    {
        for (std::size_t i = 0; i < n_points; ++i)
        {
            const std::size_t b = i * 3;
            const double r = local_points[b + 0];
            const double theta = deg_to_rad(local_points[b + 1]);
            const double local_cart[3] = {r * std::cos(theta), r * std::sin(theta), local_points[b + 2]};
            double gv[3] = {0.0, 0.0, 0.0};
            mat_vec3(_local_to_global_cart, local_cart, gv);
            global_points_out[b + 0] = _origin[0] + gv[0];
            global_points_out[b + 1] = _origin[1] + gv[1];
            global_points_out[b + 2] = _origin[2] + gv[2];
        }
    }

    SphericalCoordinateSystem::SphericalCoordinateSystem(
        const double *origin,
        const double *nxax,
        const double *nxypl,
        const double *angles_deg,
        const double *local_to_global_cart)
        : CoordinateSystem(
              origin ? origin : kDefaultOrigin,
              nullptr,
              "spherical")
    {
        const double default_origin[3] = {0.0, 0.0, 0.0};
        const double default_nxax[3] = {1.0, 0.0, 0.0};
        const double default_nxypl[3] = {0.0, 0.0, -1.0};

        const double *o = origin ? origin : default_origin;
        _origin[0] = o[0];
        _origin[1] = o[1];
        _origin[2] = o[2];

        if (local_to_global_cart)
        {
            for (int i = 0; i < 9; ++i)
                _local_to_global_cart[i] = local_to_global_cart[i];
        }
        else if (angles_deg)
        {
            build_local_to_global_matrix(angles_deg[0], angles_deg[1], angles_deg[2], _local_to_global_cart);
        }
        else
        {
            const double *x = nxax ? nxax : default_nxax;
            const double *p = nxypl ? nxypl : default_nxypl;
            build_local_to_global_from_points_raw(_origin, x, p, _local_to_global_cart);
        }
    }

    SphericalCoordinateSystem SphericalCoordinateSystem::from_angles(
        const double *origin,
        const double *angles_deg)
    {
        return SphericalCoordinateSystem(origin, nullptr, nullptr, angles_deg);
    }

    SphericalCoordinateSystem SphericalCoordinateSystem::from_points(
        const double *origin,
        const double *nxax,
        const double *nxypl)
    {
        return SphericalCoordinateSystem(origin, nxax, nxypl, nullptr);
    }

    void SphericalCoordinateSystem::vector_components_to_local_raw(
        const double *global_vectors,
        const double *global_positions,
        const std::size_t n_vectors,
        double *local_components_out) const
    {
        for (std::size_t i = 0; i < n_vectors; ++i)
        {
            const std::size_t b = i * 3;
            double lp[3] = {
                global_positions[b + 0] - _origin[0],
                global_positions[b + 1] - _origin[1],
                global_positions[b + 2] - _origin[2]};
            mat_t_vec3(_local_to_global_cart, lp, lp);

            double lv[3] = {global_vectors[b + 0], global_vectors[b + 1], global_vectors[b + 2]};
            mat_t_vec3(_local_to_global_cart, lv, lv);

            const double x = lp[0];
            const double y = lp[1];
            const double z = lp[2];
            const double r = std::sqrt(x * x + y * y + z * z);
            const double rho = std::sqrt(x * x + y * y);

            double e_r[3] = {0.0, 0.0, 1.0};
            double e_theta[3] = {1.0, 0.0, 0.0};
            double e_phi[3] = {0.0, 1.0, 0.0};

            if (r > 0.0)
            {
                const double inv_r = 1.0 / r;
                e_r[0] = x * inv_r;
                e_r[1] = y * inv_r;
                e_r[2] = z * inv_r;
                if (rho > 0.0)
                {
                    const double inv_rho = 1.0 / rho;
                    const double cos_theta = z * inv_r;
                    const double sin_theta = rho * inv_r;
                    const double cos_phi = x * inv_rho;
                    const double sin_phi = y * inv_rho;
                    e_theta[0] = cos_theta * cos_phi;
                    e_theta[1] = cos_theta * sin_phi;
                    e_theta[2] = -sin_theta;
                    e_phi[0] = -sin_phi;
                    e_phi[1] = cos_phi;
                    e_phi[2] = 0.0;
                }
            }

            local_components_out[b + 0] = dot3(e_r, lv);
            local_components_out[b + 1] = dot3(e_theta, lv);
            local_components_out[b + 2] = dot3(e_phi, lv);
        }
    }

    void SphericalCoordinateSystem::vector_components_to_global_raw(
        const double *local_components,
        const double *global_positions,
        const std::size_t n_vectors,
        double *global_vectors_out) const
    {
        for (std::size_t i = 0; i < n_vectors; ++i)
        {
            const std::size_t b = i * 3;
            double lp[3] = {
                global_positions[b + 0] - _origin[0],
                global_positions[b + 1] - _origin[1],
                global_positions[b + 2] - _origin[2]};
            mat_t_vec3(_local_to_global_cart, lp, lp);

            const double x = lp[0];
            const double y = lp[1];
            const double z = lp[2];
            const double r = std::sqrt(x * x + y * y + z * z);
            const double rho = std::sqrt(x * x + y * y);

            double e_r[3] = {0.0, 0.0, 1.0};
            double e_theta[3] = {1.0, 0.0, 0.0};
            double e_phi[3] = {0.0, 1.0, 0.0};

            if (r > 0.0)
            {
                const double inv_r = 1.0 / r;
                e_r[0] = x * inv_r;
                e_r[1] = y * inv_r;
                e_r[2] = z * inv_r;
                if (rho > 0.0)
                {
                    const double inv_rho = 1.0 / rho;
                    const double cos_theta = z * inv_r;
                    const double sin_theta = rho * inv_r;
                    const double cos_phi = x * inv_rho;
                    const double sin_phi = y * inv_rho;
                    e_theta[0] = cos_theta * cos_phi;
                    e_theta[1] = cos_theta * sin_phi;
                    e_theta[2] = -sin_theta;
                    e_phi[0] = -sin_phi;
                    e_phi[1] = cos_phi;
                    e_phi[2] = 0.0;
                }
            }

            const double lv[3] = {
                e_r[0] * local_components[b + 0] + e_theta[0] * local_components[b + 1] + e_phi[0] * local_components[b + 2],
                e_r[1] * local_components[b + 0] + e_theta[1] * local_components[b + 1] + e_phi[1] * local_components[b + 2],
                e_r[2] * local_components[b + 0] + e_theta[2] * local_components[b + 1] + e_phi[2] * local_components[b + 2]};

            mat_vec3(_local_to_global_cart, lv, &global_vectors_out[b]);
        }
    }

    void SphericalCoordinateSystem::to_local_raw(
        const double *global_points,
        const std::size_t n_points,
        double *local_points_out) const
    {
        CoordinateSystem::to_local_raw(global_points, n_points, local_points_out);
        for (std::size_t i = 0; i < n_points; ++i)
        {
            const std::size_t b = i * 3;
            const double x = local_points_out[b + 0];
            const double y = local_points_out[b + 1];
            const double z = local_points_out[b + 2];
            const double r = std::sqrt(x * x + y * y + z * z);
            double theta = 0.0;
            double phi = 0.0;
            if (r > 0.0)
            {
                theta = std::acos(std::clamp(z / r, -1.0, 1.0));
                phi = std::atan2(y, x);
            }
            local_points_out[b + 0] = r;
            local_points_out[b + 1] = rad_to_deg(theta);
            local_points_out[b + 2] = rad_to_deg(phi);
        }
    }

    void SphericalCoordinateSystem::to_global_raw(
        const double *local_points,
        const std::size_t n_points,
        double *global_points_out) const
    {
        for (std::size_t i = 0; i < n_points; ++i)
        {
            const std::size_t b = i * 3;
            const double r = local_points[b + 0];
            const double theta = deg_to_rad(local_points[b + 1]);
            const double phi = deg_to_rad(local_points[b + 2]);
            const double sin_theta = std::sin(theta);
            const double local_cart[3] = {
                r * sin_theta * std::cos(phi),
                r * sin_theta * std::sin(phi),
                r * std::cos(theta)};

            double gv[3] = {0.0, 0.0, 0.0};
            mat_vec3(_local_to_global_cart, local_cart, gv);
            global_points_out[b + 0] = _origin[0] + gv[0];
            global_points_out[b + 1] = _origin[1] + gv[1];
            global_points_out[b + 2] = _origin[2] + gv[2];
        }
    }

} // namespace apec_m

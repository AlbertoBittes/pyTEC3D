#ifndef APEC_M_COORDINATE_SYSTEM_H
#define APEC_M_COORDINATE_SYSTEM_H

#include <cstddef>
#include <string>

namespace apec_m
{

    class CoordinateSystem
    {
    public:
        CoordinateSystem(const double *origin,
                         const double *local_to_global_cart,
                         std::string cs_type);
        virtual ~CoordinateSystem() = default;

        const std::string &type() const noexcept;
        void origin_raw(double *origin_out) const noexcept;
        void local_to_global_cart_raw(double *matrix_out) const noexcept;
        void axis_direction_raw(double *axis_out) const noexcept;

        virtual void to_local_raw(
            const double *global_points,
            std::size_t n_points,
            double *local_points_out) const;

        virtual void to_global_raw(
            const double *local_points,
            std::size_t n_points,
            double *global_points_out) const;

        virtual void vector_components_to_local_raw(
            const double *global_vectors,
            const double *global_positions,
            std::size_t n_vectors,
            double *local_components_out) const;

        virtual void vector_components_to_global_raw(
            const double *local_components,
            const double *global_positions,
            std::size_t n_vectors,
            double *global_vectors_out) const;

        void local_to_global_matrix_raw(
            const double *global_position,
            double *matrix_out) const;

        void rotate_vector_to_global_raw(
            const double *local_components,
            const double *global_position,
            double *global_components_out) const;

        void rotate_tensor_to_global_raw(
            const double *local_tensor,
            const double *global_position,
            double *global_tensor_out) const;

        static void build_local_to_global_matrix(
            double angle_xy_deg,
            double angle_yz_deg,
            double angle_xz_deg,
            double *matrix_out);

        static void build_local_to_global_coefficients_raw(
            const double *angles_deg,
            std::size_t n_rows,
            double *coeffs_out);

        static void build_global_to_local_matrix(
            double angle_xy_deg,
            double angle_yz_deg,
            double angle_xz_deg,
            double *matrix_out);

        static void build_global_to_local(
            double angle_xy_deg,
            double angle_yz_deg,
            double angle_xz_deg,
            double *matrix_out);

    protected:
        static constexpr double kFiniteTolerance = 1e-12;

        static double deg_to_rad(double value) noexcept;
        static double rad_to_deg(double value) noexcept;

        static void normalize3_or_error(
            const double *vec,
            double *out,
            const std::string &error_message);

        static void build_local_to_global_from_points_raw(
            const double *origin,
            const double *nxax,
            const double *nxypl,
            double *matrix_out);

        static double dot3(const double *a, const double *b) noexcept;
        static void cross3(const double *a, const double *b, double *out) noexcept;

        double _origin[3];
        double _local_to_global_cart[9];
        std::string _type;
    };

    class CartesianCoordinateSystem : public CoordinateSystem
    {
    public:
        explicit CartesianCoordinateSystem(
            const double *origin = nullptr,
            const double *nxax = nullptr,
            const double *nxypl = nullptr,
            const double *angles_deg = nullptr,
            const double *local_to_global_cart = nullptr);

        static CartesianCoordinateSystem from_angles(
            const double *origin,
            const double *angles_deg);

        static CartesianCoordinateSystem from_points(
            const double *origin,
            const double *nxax,
            const double *nxypl);
    };

    class CylindricalCoordinateSystem : public CoordinateSystem
    {
    public:
        explicit CylindricalCoordinateSystem(
            const double *origin = nullptr,
            const double *nxax = nullptr,
            const double *nxypl = nullptr,
            const double *angles_deg = nullptr,
            const double *local_to_global_cart = nullptr);

        static CylindricalCoordinateSystem from_angles(
            const double *origin,
            const double *angles_deg);

        static CylindricalCoordinateSystem from_points(
            const double *origin,
            const double *nxax,
            const double *nxypl);

        void vector_components_to_local_raw(
            const double *global_vectors,
            const double *global_positions,
            std::size_t n_vectors,
            double *local_components_out) const override;

        void vector_components_to_global_raw(
            const double *local_components,
            const double *global_positions,
            std::size_t n_vectors,
            double *global_vectors_out) const override;

        void to_local_raw(
            const double *global_points,
            std::size_t n_points,
            double *local_points_out) const override;

        void to_global_raw(
            const double *local_points,
            std::size_t n_points,
            double *global_points_out) const override;
    };

    // Rotation axis (local Z) aligned with global X.
    class XAxisCylindricalCoordinateSystem : public CylindricalCoordinateSystem
    {
    public:
        explicit XAxisCylindricalCoordinateSystem(
            const double *origin = nullptr,
            const double *local_to_global_cart = nullptr);
    };

    // Rotation axis (local Z) aligned with global Y.
    class YAxisCylindricalCoordinateSystem : public CylindricalCoordinateSystem
    {
    public:
        explicit YAxisCylindricalCoordinateSystem(
            const double *origin = nullptr,
            const double *local_to_global_cart = nullptr);
    };

    // Rotation axis (local Z) aligned with global Z (ANSYS APDL LOCAL default).
    class ZAxisCylindricalCoordinateSystem : public CylindricalCoordinateSystem
    {
    public:
        explicit ZAxisCylindricalCoordinateSystem(
            const double *origin = nullptr,
            const double *local_to_global_cart = nullptr);
    };

    class SphericalCoordinateSystem : public CoordinateSystem
    {
    public:
        explicit SphericalCoordinateSystem(
            const double *origin = nullptr,
            const double *nxax = nullptr,
            const double *nxypl = nullptr,
            const double *angles_deg = nullptr,
            const double *local_to_global_cart = nullptr);

        static SphericalCoordinateSystem from_angles(
            const double *origin,
            const double *angles_deg);

        static SphericalCoordinateSystem from_points(
            const double *origin,
            const double *nxax,
            const double *nxypl);

        void vector_components_to_local_raw(
            const double *global_vectors,
            const double *global_positions,
            std::size_t n_vectors,
            double *local_components_out) const override;

        void vector_components_to_global_raw(
            const double *local_components,
            const double *global_positions,
            std::size_t n_vectors,
            double *global_vectors_out) const override;

        void to_local_raw(
            const double *global_points,
            std::size_t n_points,
            double *local_points_out) const override;

        void to_global_raw(
            const double *local_points,
            std::size_t n_points,
            double *global_points_out) const override;
    };

} // namespace apec_m

#endif // APEC_M_COORDINATE_SYSTEM_H

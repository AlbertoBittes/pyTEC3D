#include "apec_m_coordinate/coordinate_system.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <memory>
#include <optional>
#include <sstream>
#include <string>

namespace py = pybind11;
using namespace pybind11::literals;

namespace
{
    std::size_t points_count_or_throw(
        const py::array_t<double, py::array::c_style | py::array::forcecast> &array,
        const std::string &arg_name)
    {
        if (array.ndim() == 1)
        {
            if (array.shape(0) != 3)
            {
                std::ostringstream msg;
                msg << arg_name << " must have shape (N, 3) or (3,), got (" << array.shape(0) << ").";
                throw std::invalid_argument(msg.str());
            }
            return 1;
        }

        if (array.ndim() != 2 || array.shape(1) != 3)
        {
            std::ostringstream msg;
            msg << arg_name << " must have shape (N, 3) or (3,), got ndim=" << array.ndim() << ".";
            throw std::invalid_argument(msg.str());
        }

        return static_cast<std::size_t>(array.shape(0));
    }

    py::array_t<double> allocate_points3_array(const std::size_t n_rows)
    {
        return py::array_t<double>({static_cast<py::ssize_t>(n_rows), static_cast<py::ssize_t>(3)});
    }

    py::array_t<double> point3_out_from_raw(const double *raw3)
    {
        py::array_t<double> out({static_cast<py::ssize_t>(3)});
        auto view = out.mutable_unchecked<1>();
        view(0) = raw3[0];
        view(1) = raw3[1];
        view(2) = raw3[2];
        return out;
    }

    py::array_t<double> matrix3_out_from_raw(const double *raw9)
    {
        py::array_t<double> out({static_cast<py::ssize_t>(3), static_cast<py::ssize_t>(3)});
        auto view = out.mutable_unchecked<2>();
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
                view(i, j) = raw9[i * 3 + j];
        }
        return out;
    }

    const double *optional_point3_ptr(const py::object &obj, double *storage)
    {
        if (obj.is_none())
            return nullptr;

        const auto arr = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(obj);
        const std::size_t n = points_count_or_throw(arr, "point");
        if (n != 1)
            throw std::invalid_argument("Point arguments must have shape (3,).\n");

        storage[0] = arr.data()[0];
        storage[1] = arr.data()[1];
        storage[2] = arr.data()[2];
        return storage;
    }

    const double *optional_matrix3_ptr(const py::object &obj, double *storage)
    {
        if (obj.is_none())
            return nullptr;

        const auto arr = py::cast<py::array_t<double, py::array::c_style | py::array::forcecast>>(obj);
        if (arr.ndim() == 1)
        {
            if (arr.shape(0) != 9)
                throw std::invalid_argument("local_to_global_cart must have shape (3, 3) or (9,).");
            for (int i = 0; i < 9; ++i)
                storage[i] = arr.data()[i];
            return storage;
        }

        if (arr.ndim() != 2 || arr.shape(0) != 3 || arr.shape(1) != 3)
            throw std::invalid_argument("local_to_global_cart must have shape (3, 3) or (9,).");

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                storage[i * 3 + j] = arr.data()[i * 3 + j];

        return storage;
    }

    const double *required_point3_ptr(
        const py::array_t<double, py::array::c_style | py::array::forcecast> &arr,
        const std::string &arg_name,
        double *storage)
    {
        const std::size_t n = points_count_or_throw(arr, arg_name);
        if (n != 1)
            throw std::invalid_argument(arg_name + " must have shape (3,).\n");

        storage[0] = arr.data()[0];
        storage[1] = arr.data()[1];
        storage[2] = arr.data()[2];
        return storage;
    }

} // namespace

void bind_coordinate_system(py::module_ &module)
{
    py::class_<apec_m::CoordinateSystem, std::shared_ptr<apec_m::CoordinateSystem>> coord(module, "CoordinateSystem");

    coord
        .def(py::init([](const py::object &origin, const py::object &local_to_global_cart, const std::string &cs_type)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);
                          return apec_m::CoordinateSystem(o_ptr, m_ptr, cs_type); }),
             "origin"_a = py::none(),
             "local_to_global_cart"_a = py::none(),
             "cs_type"_a = "cartesian")
        .def_property_readonly("type", &apec_m::CoordinateSystem::type)
        .def_property_readonly("origin", [](const apec_m::CoordinateSystem &self)
                               {
                                   double out[3] = {0.0, 0.0, 0.0};
                                   self.origin_raw(out);
                                   return point3_out_from_raw(out); })
        .def_property_readonly("local_to_global_cart", [](const apec_m::CoordinateSystem &self)
                               {
                                   double out[9] = {0.0};
                                   self.local_to_global_cart_raw(out);
                                   return matrix3_out_from_raw(out); })
        .def("axis_direction", [](const apec_m::CoordinateSystem &self)
             {
                 double out[3] = {0.0, 0.0, 0.0};
                 self.axis_direction_raw(out);
                 return point3_out_from_raw(out); })
        .def("to_local", [](const apec_m::CoordinateSystem &self, const py::array_t<double, py::array::c_style | py::array::forcecast> &global_points)
             {
                 const std::size_t n_points = points_count_or_throw(global_points, "global_points");
                 py::array_t<double> out = allocate_points3_array(n_points);
                 self.to_local_raw(global_points.data(), n_points, out.mutable_data());
                 return out; }, "global_points"_a)
        .def("to_global", [](const apec_m::CoordinateSystem &self, const py::array_t<double, py::array::c_style | py::array::forcecast> &local_points)
             {
                 const std::size_t n_points = points_count_or_throw(local_points, "local_points");
                 py::array_t<double> out = allocate_points3_array(n_points);
                 self.to_global_raw(local_points.data(), n_points, out.mutable_data());
                 return out; }, "local_points"_a)
        .def("vector_components_to_local", [](const apec_m::CoordinateSystem &self, const py::array_t<double, py::array::c_style | py::array::forcecast> &global_vectors, const py::array_t<double, py::array::c_style | py::array::forcecast> &global_positions)
             {
                 const std::size_t n_vectors = points_count_or_throw(global_vectors, "global_vectors");
                 const std::size_t n_positions = points_count_or_throw(global_positions, "global_positions");
                 if (n_vectors != n_positions)
                     throw std::invalid_argument("global_vectors and global_positions must have the same number of rows.");
                 py::array_t<double> out = allocate_points3_array(n_vectors);
                 self.vector_components_to_local_raw(global_vectors.data(), global_positions.data(), n_vectors, out.mutable_data());
                 return out; }, "global_vectors"_a, "global_positions"_a)
        .def("vector_components_to_global", [](const apec_m::CoordinateSystem &self, const py::array_t<double, py::array::c_style | py::array::forcecast> &local_components, const py::array_t<double, py::array::c_style | py::array::forcecast> &global_positions)
             {
                 const std::size_t n_vectors = points_count_or_throw(local_components, "local_components");
                 const std::size_t n_positions = points_count_or_throw(global_positions, "global_positions");
                 if (n_vectors != n_positions)
                     throw std::invalid_argument("local_components and global_positions must have the same number of rows.");
                 py::array_t<double> out = allocate_points3_array(n_vectors);
                 self.vector_components_to_global_raw(local_components.data(), global_positions.data(), n_vectors, out.mutable_data());
                 return out; }, "local_components"_a, "global_positions"_a)
        .def_static("build_local_to_global_matrix", [](const double angle_xy_deg, const double angle_yz_deg, const double angle_xz_deg)
                    {
                        double out[9] = {0.0};
                        apec_m::CoordinateSystem::build_local_to_global_matrix(angle_xy_deg, angle_yz_deg, angle_xz_deg, out);
                        return matrix3_out_from_raw(out); }, "angle_xy_deg"_a, "angle_yz_deg"_a, "angle_xz_deg"_a)
        .def_static("build_local_to_global_coefficients", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &angles_deg)
                    {
                        const std::size_t n_rows = points_count_or_throw(angles_deg, "angles_deg");
                        py::array_t<double> out({static_cast<py::ssize_t>(n_rows), static_cast<py::ssize_t>(9)});
                        apec_m::CoordinateSystem::build_local_to_global_coefficients_raw(angles_deg.data(), n_rows, out.mutable_data());
                        return out; }, "angles_deg"_a)
        .def_static("build_global_to_local_matrix", [](const double angle_xy_deg, const double angle_yz_deg, const double angle_xz_deg)
                    {
                        double out[9] = {0.0};
                        apec_m::CoordinateSystem::build_global_to_local_matrix(angle_xy_deg, angle_yz_deg, angle_xz_deg, out);
                        return matrix3_out_from_raw(out); }, "angle_xy_deg"_a, "angle_yz_deg"_a, "angle_xz_deg"_a)
        .def_static("build_global_to_local", [](const double angle_xy_deg, const double angle_yz_deg, const double angle_xz_deg)
                    {
                        double out[9] = {0.0};
                        apec_m::CoordinateSystem::build_global_to_local(angle_xy_deg, angle_yz_deg, angle_xz_deg, out);
                        return matrix3_out_from_raw(out); }, "angle_xy_deg"_a, "angle_yz_deg"_a, "angle_xz_deg"_a);

    py::class_<apec_m::CartesianCoordinateSystem, apec_m::CoordinateSystem, std::shared_ptr<apec_m::CartesianCoordinateSystem>>(module, "CartesianCoordinateSystem")
        .def(py::init([](const py::object &origin,
                         const py::object &nxax,
                         const py::object &nxypl,
                         const py::object &local_to_global_cart,
                         const py::object &thxy_deg,
                         const py::object &thyz_deg,
                         const py::object &thzx_deg)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double x[3] = {0.0, 0.0, 0.0};
                          double p[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *x_ptr = optional_point3_ptr(nxax, x);
                          const double *p_ptr = optional_point3_ptr(nxypl, p);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);

                          double angles[3] = {0.0, 0.0, 0.0};
                          const bool has_xy = !thxy_deg.is_none();
                          const bool has_yz = !thyz_deg.is_none();
                          const bool has_zx = !thzx_deg.is_none();
                          const bool any = has_xy || has_yz || has_zx;
                          const bool all = has_xy && has_yz && has_zx;
                          if (any && !all)
                              throw std::invalid_argument("thxy_deg, thyz_deg and thzx_deg must all be provided together.");
                          if (m_ptr && (x_ptr || p_ptr || any))
                              throw std::invalid_argument("local_to_global_cart cannot be combined with nxax/nxypl/angles.");
                          const double *angles_ptr = nullptr;
                          if (all)
                          {
                              angles[0] = py::cast<double>(thxy_deg);
                              angles[1] = py::cast<double>(thyz_deg);
                              angles[2] = py::cast<double>(thzx_deg);
                              angles_ptr = angles;
                          }

                          return apec_m::CartesianCoordinateSystem(o_ptr, x_ptr, p_ptr, angles_ptr, m_ptr); }),
             "origin"_a = py::none(),
             "nxax"_a = py::none(),
             "nxypl"_a = py::none(),
             "local_to_global_cart"_a = py::none(),
             "thxy_deg"_a = py::none(),
             "thyz_deg"_a = py::none(),
             "thzx_deg"_a = py::none())
        .def_static("from_angles", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &origin, const py::array_t<double, py::array::c_style | py::array::forcecast> &angles_deg)
                    {
                        double o[3] = {0.0, 0.0, 0.0};
                        double a[3] = {0.0, 0.0, 0.0};
                        return apec_m::CartesianCoordinateSystem::from_angles(
                            required_point3_ptr(origin, "origin", o),
                            required_point3_ptr(angles_deg, "angles_deg", a)); }, "origin"_a, "angles_deg"_a)
        .def_static("from_points", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &origin, const py::array_t<double, py::array::c_style | py::array::forcecast> &nxax, const py::array_t<double, py::array::c_style | py::array::forcecast> &nxypl)
                    {
                        double o[3] = {0.0, 0.0, 0.0};
                        double x[3] = {0.0, 0.0, 0.0};
                        double p[3] = {0.0, 0.0, 0.0};
                        return apec_m::CartesianCoordinateSystem::from_points(
                            required_point3_ptr(origin, "origin", o),
                            required_point3_ptr(nxax, "nxax", x),
                            required_point3_ptr(nxypl, "nxypl", p)); }, "origin"_a, "nxax"_a, "nxypl"_a);

    py::class_<apec_m::CylindricalCoordinateSystem, apec_m::CoordinateSystem, std::shared_ptr<apec_m::CylindricalCoordinateSystem>>(module, "CylindricalCoordinateSystem")
        .def(py::init([](const py::object &origin,
                         const py::object &nxax,
                         const py::object &nxypl,
                         const py::object &local_to_global_cart,
                         const py::object &thxy_deg,
                         const py::object &thyz_deg,
                         const py::object &thzx_deg)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double x[3] = {0.0, 0.0, 0.0};
                          double p[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *x_ptr = optional_point3_ptr(nxax, x);
                          const double *p_ptr = optional_point3_ptr(nxypl, p);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);

                          double angles[3] = {0.0, 0.0, 0.0};
                          const bool has_xy = !thxy_deg.is_none();
                          const bool has_yz = !thyz_deg.is_none();
                          const bool has_zx = !thzx_deg.is_none();
                          const bool any = has_xy || has_yz || has_zx;
                          const bool all = has_xy && has_yz && has_zx;
                          if (any && !all)
                              throw std::invalid_argument("thxy_deg, thyz_deg and thzx_deg must all be provided together.");
                          if (m_ptr && (x_ptr || p_ptr || any))
                              throw std::invalid_argument("local_to_global_cart cannot be combined with nxax/nxypl/angles.");
                          const double *angles_ptr = nullptr;
                          if (all)
                          {
                              angles[0] = py::cast<double>(thxy_deg);
                              angles[1] = py::cast<double>(thyz_deg);
                              angles[2] = py::cast<double>(thzx_deg);
                              angles_ptr = angles;
                          }

                          return apec_m::CylindricalCoordinateSystem(o_ptr, x_ptr, p_ptr, angles_ptr, m_ptr); }),
             "origin"_a = py::none(),
             "nxax"_a = py::none(),
             "nxypl"_a = py::none(),
             "local_to_global_cart"_a = py::none(),
             "thxy_deg"_a = py::none(),
             "thyz_deg"_a = py::none(),
             "thzx_deg"_a = py::none())
        .def_static("from_angles", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &origin, const py::array_t<double, py::array::c_style | py::array::forcecast> &angles_deg)
                    {
                        double o[3] = {0.0, 0.0, 0.0};
                        double a[3] = {0.0, 0.0, 0.0};
                        return apec_m::CylindricalCoordinateSystem::from_angles(
                            required_point3_ptr(origin, "origin", o),
                            required_point3_ptr(angles_deg, "angles_deg", a)); }, "origin"_a, "angles_deg"_a)
        .def_static("from_points", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &origin, const py::array_t<double, py::array::c_style | py::array::forcecast> &nxax, const py::array_t<double, py::array::c_style | py::array::forcecast> &nxypl)
                    {
                        double o[3] = {0.0, 0.0, 0.0};
                        double x[3] = {0.0, 0.0, 0.0};
                        double p[3] = {0.0, 0.0, 0.0};
                        return apec_m::CylindricalCoordinateSystem::from_points(
                            required_point3_ptr(origin, "origin", o),
                            required_point3_ptr(nxax, "nxax", x),
                            required_point3_ptr(nxypl, "nxypl", p)); }, "origin"_a, "nxax"_a, "nxypl"_a);

    py::class_<apec_m::XAxisCylindricalCoordinateSystem, apec_m::CylindricalCoordinateSystem, std::shared_ptr<apec_m::XAxisCylindricalCoordinateSystem>>(module, "XAxisCylindricalCoordinateSystem")
        .def(py::init([](const py::object &origin, const py::object &local_to_global_cart)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);
                          return apec_m::XAxisCylindricalCoordinateSystem(o_ptr, m_ptr); }),
             "origin"_a = py::none(),
             "local_to_global_cart"_a = py::none());

    py::class_<apec_m::YAxisCylindricalCoordinateSystem, apec_m::CylindricalCoordinateSystem, std::shared_ptr<apec_m::YAxisCylindricalCoordinateSystem>>(module, "YAxisCylindricalCoordinateSystem")
        .def(py::init([](const py::object &origin, const py::object &local_to_global_cart)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);
                          return apec_m::YAxisCylindricalCoordinateSystem(o_ptr, m_ptr); }),
             "origin"_a = py::none(),
             "local_to_global_cart"_a = py::none());

    py::class_<apec_m::ZAxisCylindricalCoordinateSystem, apec_m::CylindricalCoordinateSystem, std::shared_ptr<apec_m::ZAxisCylindricalCoordinateSystem>>(module, "ZAxisCylindricalCoordinateSystem")
        .def(py::init([](const py::object &origin, const py::object &local_to_global_cart)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);
                          return apec_m::ZAxisCylindricalCoordinateSystem(o_ptr, m_ptr); }),
             "origin"_a = py::none(),
             "local_to_global_cart"_a = py::none());

    py::class_<apec_m::SphericalCoordinateSystem, apec_m::CoordinateSystem, std::shared_ptr<apec_m::SphericalCoordinateSystem>>(module, "SphericalCoordinateSystem")
        .def(py::init([](const py::object &origin,
                         const py::object &nxax,
                         const py::object &nxypl,
                         const py::object &local_to_global_cart,
                         const py::object &thxy_deg,
                         const py::object &thyz_deg,
                         const py::object &thzx_deg)
                      {
                          double o[3] = {0.0, 0.0, 0.0};
                          double x[3] = {0.0, 0.0, 0.0};
                          double p[3] = {0.0, 0.0, 0.0};
                          double m[9] = {0.0};
                          const double *o_ptr = optional_point3_ptr(origin, o);
                          const double *x_ptr = optional_point3_ptr(nxax, x);
                          const double *p_ptr = optional_point3_ptr(nxypl, p);
                          const double *m_ptr = optional_matrix3_ptr(local_to_global_cart, m);

                          double angles[3] = {0.0, 0.0, 0.0};
                          const bool has_xy = !thxy_deg.is_none();
                          const bool has_yz = !thyz_deg.is_none();
                          const bool has_zx = !thzx_deg.is_none();
                          const bool any = has_xy || has_yz || has_zx;
                          const bool all = has_xy && has_yz && has_zx;
                          if (any && !all)
                              throw std::invalid_argument("thxy_deg, thyz_deg and thzx_deg must all be provided together.");
                          if (m_ptr && (x_ptr || p_ptr || any))
                              throw std::invalid_argument("local_to_global_cart cannot be combined with nxax/nxypl/angles.");
                          const double *angles_ptr = nullptr;
                          if (all)
                          {
                              angles[0] = py::cast<double>(thxy_deg);
                              angles[1] = py::cast<double>(thyz_deg);
                              angles[2] = py::cast<double>(thzx_deg);
                              angles_ptr = angles;
                          }

                          return apec_m::SphericalCoordinateSystem(o_ptr, x_ptr, p_ptr, angles_ptr, m_ptr); }),
             "origin"_a = py::none(),
             "nxax"_a = py::none(),
             "nxypl"_a = py::none(),
             "local_to_global_cart"_a = py::none(),
             "thxy_deg"_a = py::none(),
             "thyz_deg"_a = py::none(),
             "thzx_deg"_a = py::none())
        .def_static("from_angles", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &origin, const py::array_t<double, py::array::c_style | py::array::forcecast> &angles_deg)
                    {
                        double o[3] = {0.0, 0.0, 0.0};
                        double a[3] = {0.0, 0.0, 0.0};
                        return apec_m::SphericalCoordinateSystem::from_angles(
                            required_point3_ptr(origin, "origin", o),
                            required_point3_ptr(angles_deg, "angles_deg", a)); }, "origin"_a, "angles_deg"_a)
        .def_static("from_points", [](const py::array_t<double, py::array::c_style | py::array::forcecast> &origin, const py::array_t<double, py::array::c_style | py::array::forcecast> &nxax, const py::array_t<double, py::array::c_style | py::array::forcecast> &nxypl)
                    {
                        double o[3] = {0.0, 0.0, 0.0};
                        double x[3] = {0.0, 0.0, 0.0};
                        double p[3] = {0.0, 0.0, 0.0};
                        return apec_m::SphericalCoordinateSystem::from_points(
                            required_point3_ptr(origin, "origin", o),
                            required_point3_ptr(nxax, "nxax", x),
                            required_point3_ptr(nxypl, "nxypl", p)); }, "origin"_a, "nxax"_a, "nxypl"_a);
}

/**
 * @file main.cpp
 * @brief pybind11 module entry point that exposes ApecMSolver to Python.
 *
 * Defines the @c apec_m_solver Python extension module.  The module exports:
 *   - @c ApecMSolver class with @c load_model(), @c solve(), @c close(), and
 *     @c is_loaded() bindings.
 *   - @c shutdown_runtime() static function to cleanly finalise the libMesh
 *     and spdlog runtimes after all solver instances have been destroyed.
 */
#include "apec_m_solver.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_coordinate_system(py::module_ &module);

PYBIND11_MODULE(apec_m_solver, module)
{
     module.doc() = "APEC-M solver module (Python-owned model, C++ execution core).";
     module.def("shutdown_runtime", &ApecMSolver::shutdown_runtime);
     bind_coordinate_system(module);

     py::class_<ApecMSolver>(module, "ApecMSolver")
         .def(py::init<py::object>(),
              py::arg("mpi_comm") = py::none())
         .def("load_model",
              &ApecMSolver::load_model,
              py::arg("apec_object"))
         .def("close", &ApecMSolver::close)
         .def("solve", &ApecMSolver::solve)
         .def("is_loaded", &ApecMSolver::is_loaded);
}

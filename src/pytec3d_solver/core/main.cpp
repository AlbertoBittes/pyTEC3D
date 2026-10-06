/**
 * @file main.cpp
 * @brief pybind11 module entry point that exposes pyTEC3DSolver to Python.
 */
#include "pytec3d_solver.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(pytec3d_solver, module)
{
     module.doc() = "pyTEC3D solver module (Python-owned model, C++ execution core).";
     module.def("shutdown_runtime", &pyTEC3DSolver::shutdown_runtime);
     module.def("registries", &pyTEC3DSolver::registries,
                "Domains, fields, element types, material properties and real constants compiled into the module.");

     py::class_<pyTEC3DSolver>(module, "pyTEC3DSolver")
         .def(py::init<int>(),
              py::arg("num_threads_per_mpi_process") = 1)
         .def("load_model",
              &pyTEC3DSolver::load_model,
              py::arg("pytec3d_object"),
              py::arg("solver_settings"))
         .def("close", &pyTEC3DSolver::close)
         .def("solve", &pyTEC3DSolver::solve, py::arg("time_step"))
         .def("is_loaded", &pyTEC3DSolver::is_loaded);
}

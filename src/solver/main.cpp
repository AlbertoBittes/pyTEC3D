/**
 * @file main.cpp
 * @brief pybind11 module entry point that exposes ApecMSolver to Python.
 */
#include "apec_m_solver/apec_m_solver.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(apec_m_solver, module)
{
     module.doc() = "APEC-M solver module (Python-owned model, C++ execution core).";
     module.def("shutdown_runtime", &ApecMSolver::shutdown_runtime);

     py::class_<ApecMSolver>(module, "ApecMSolver")
         .def(py::init<int>(),
              py::arg("num_threads_per_mpi_process") = 1)
         .def("load_model",
              &ApecMSolver::load_model,
              py::arg("apec_object"))
         .def("close", &ApecMSolver::close)
         .def("solve", &ApecMSolver::solve)
         .def("is_loaded", &ApecMSolver::is_loaded);
}

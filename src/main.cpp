#include "apec_m_solver.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(apec_m_solver, module)
{
    module.doc() = "APEC-M solver module (Python-owned model, C++ execution core).";

    py::class_<ApecMSolver>(module, "ApecMSolver")
        .def(py::init<>())
        .def("load_model",
             &ApecMSolver::load_model,
             py::arg("apec_object"),
             py::arg("mpi_comm") = py::none())
        .def("solve", &ApecMSolver::solve)
        .def("is_loaded", &ApecMSolver::is_loaded);
}

#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_coordinate_system(py::module_ &module);

PYBIND11_MODULE(apec_m_coordinate_system, module)
{
    module.doc() = "APEC-M coordinate system module.";
    bind_coordinate_system(module);
}

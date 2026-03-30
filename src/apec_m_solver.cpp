#include "apec_m_solver.h"

#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "augment_sparsity_on_nodes.h"
#include "element_types_factory.h"
#include "vtkio_updated.h"

#include <libmesh/distributed_mesh.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_modification.h>

#include <mpi.h>

#include <vtkUnstructuredGrid.h>

#include <pybind11/stl.h>

#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;

namespace
{
    py::object active_vtk_grid_from_apec(const py::object &apec_object)
    {
        if (!py::hasattr(apec_object, "mesh"))
            throw std::runtime_error("Expected an APEC-like Python object exposing a 'mesh' property.");

        return apec_object.attr("mesh");
    }

    vtkUnstructuredGrid *vtk_grid_from_python(const py::object &vtk_grid)
    {
        if (!py::hasattr(vtk_grid, "GetAddressAsString"))
            throw std::runtime_error("Expected a vtkUnstructuredGrid Python object with GetAddressAsString().");

        const std::string address = py::cast<std::string>(vtk_grid.attr("GetAddressAsString")("vtkUnstructuredGrid"));
        const std::size_t prefix_pos = address.find("0x");
        if (prefix_pos == std::string::npos)
            throw std::runtime_error("Could not parse VTK pointer address from GetAddressAsString().");

        const char *hex_begin = address.c_str() + prefix_pos + 2;
        char *hex_end = nullptr;
        const auto ptr_value = std::strtoull(hex_begin, &hex_end, 16);
        if (ptr_value == 0)
            throw std::runtime_error("Parsed null VTK pointer from Python object.");

        return reinterpret_cast<vtkUnstructuredGrid *>(ptr_value);
    }

    MPI_Comm mpi_comm_from_python(const py::object &mpi_comm)
    {
        if (mpi_comm.is_none())
            return MPI_COMM_WORLD;

        if (!py::hasattr(mpi_comm, "py2f"))
            throw std::runtime_error("Expected an mpi4py communicator with py2f().");

        const auto f_comm = static_cast<MPI_Fint>(py::cast<int>(mpi_comm.attr("py2f")()));
        return MPI_Comm_f2c(f_comm);
    }
} // namespace

void ApecMSolver::ensure_libmesh_initialized(MPI_Comm mpi_comm)
{
    if (libmesh_init)
        return;

    libmesh_init = std::make_unique<libMesh::LibMeshInit>(argc, argv, mpi_comm);
}

ApecMSolver::ApecMSolver()
{
}

ApecMSolver::~ApecMSolver() = default;

void ApecMSolver::load_model(py::object apec_object,
                             py::object mpi_comm)
{
    mpi_comm_object = mpi_comm;
    this->apec_object = apec_object;
    element_types_factory = std::make_shared<ElementTypesFactory>(apec_object.attr("element_types"));
    const MPI_Comm active_comm = mpi_comm_from_python(mpi_comm);
    ensure_libmesh_initialized(active_comm);

    mesh = std::make_unique<libMesh::DistributedMesh>(libmesh_init->comm());
    py::object vtk_grid = active_vtk_grid_from_apec(apec_object);
    vtkUnstructuredGrid *grid_ptr = vtk_grid_from_python(vtk_grid);
    active_vtk_grid = grid_ptr;
    VTKIO_Updated vtkio(*mesh, grid_ptr);
    vtkio.read(grid_ptr);

    libMesh::MeshTools::Modification::orient_elements(*mesh);

    equation_system = std::make_shared<ApecMEquationSystems>(*mesh);
    equation_system->set_parent_object(this);
    equation_system->configure_systems();
    equation_system->initialize_systems();

    auto &main_system = equation_system->get_system<ApecMLinearSystem>(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    coupling_functor = std::make_unique<AugmentSparsityOnNodes>(*mesh,
                                                                grid_ptr);
    main_system.get_dof_map().add_coupling_functor(*coupling_functor);
    libMesh::out << "Added DOF coupling sparsity functor for active mesh." << std::endl;

    mesh->prepare_for_use(true);
    mesh->print_info();

    equation_system->initialize_physics();

    loaded = true;
}

py::dict ApecMSolver::solve()
{
    if (!loaded || !mesh)
        throw std::runtime_error("Solver has no loaded model. Call load_model() first.");

    py::dict result;
    result["status"] = "model_loaded";
    result["n_nodes"] = mesh->n_nodes();
    result["n_elements"] = mesh->n_elem();
    result["mesh_dimension"] = mesh->mesh_dimension();
    if (apec_object && !apec_object.is_none())
    {
        if (py::hasattr(apec_object, "active_mesh_name"))
            result["active_mesh_name"] = apec_object.attr("active_mesh_name");
        if (py::hasattr(apec_object, "element_types"))
            result["n_element_types"] = py::len(apec_object.attr("element_types"));
        if (py::hasattr(apec_object, "materials"))
            result["n_materials"] = py::len(apec_object.attr("materials"));
        if (py::hasattr(apec_object, "real_constants"))
            result["n_real_constants"] = py::len(apec_object.attr("real_constants"));
        if (py::hasattr(apec_object, "dof_couplings"))
            result["n_dof_couplings"] = py::len(apec_object.attr("dof_couplings"));
    }

    write_system_to_active_vtk_grid(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    write_system_to_active_vtk_grid(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::NODAL_SOLUTION));
    write_system_to_active_vtk_grid(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::ELEMENT_SOLUTION));

    // Next step: initialize ApecMEquationSystems from apec_object,
    // then run the real solve loop and write back output arrays.
    return result;
}

void ApecMSolver::write_system_to_active_vtk_grid(const std::string &system_name)
{
    if (!loaded || !mesh || !equation_system || !active_vtk_grid)
        throw std::runtime_error("Solver is not fully initialized. Call load_model() first.");

    VTKIO_Updated vtkio(*mesh, equation_system.get(), active_vtk_grid);
    vtkio.write_system_grid(system_name);
}

bool ApecMSolver::is_loaded() const
{
    return loaded;
}

py::object ApecMSolver::get_apec_object() const
{
    return apec_object;
}

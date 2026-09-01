#include "apec_m_model.h"

#include "apec_m_equation_systems.h"
#include "apec_m_logging.h"
#include "apec_m_system.h"
#include "augment_sparsity_on_nodes.h"
#include "vtkio_updated.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include <vtkObjectBase.h>
#include <vtkPythonUtil.h>
#include <vtkUnstructuredGrid.h>

namespace py = pybind11;

namespace
{
    py::object active_vtk_grid_from_apec(const py::object &apec_object)
    {
        if (!py::hasattr(apec_object, "mesh"))
            apec::logging::error_and_abort("Expected an APEC-like Python model object exposing a 'mesh' property.");

        return apec_object.attr("mesh");
    }

    vtkUnstructuredGrid *vtk_grid_from_python(const py::object &vtk_grid)
    {
        if (vtk_grid.is_none())
            apec::logging::error_and_abort("Expected a vtkUnstructuredGrid Python object, received None.");

        vtkObjectBase *vtk_base = vtkPythonUtil::GetPointerFromObject(vtk_grid.ptr(), "vtkUnstructuredGrid");
        if (!vtk_base)
            apec::logging::error_and_abort("Expected a vtkUnstructuredGrid Python object compatible with vtkPythonUtil.");

        auto *grid = vtkUnstructuredGrid::SafeDownCast(vtk_base);
        if (!grid)
            apec::logging::error_and_abort("Python object is not a valid vtkUnstructuredGrid instance.");

        return grid;
    }
}

ApecMModel::~ApecMModel()
{
    clear();
}

void ApecMModel::load_into_solver(const std::shared_ptr<libMesh::LibMeshInit> &libmesh_init,
                                  py::object apec_model_object)
{
    if (!libmesh_init)
        apec::logging::error_and_abort("Cannot load model without an initialized libMesh runtime.");

    clear();
    closed_ = false;

    apec::logging::info("Loading model into solver runtime");
    apec_object = std::move(apec_model_object);

    mesh = std::make_unique<libMesh::DistributedMesh>(libmesh_init->comm());
    py::object vtk_grid = active_vtk_grid_from_apec(apec_object);
    vtkUnstructuredGrid *grid_ptr = vtk_grid_from_python(vtk_grid);
    active_vtk_grid = grid_ptr;

    VTKIO_Updated vtkio(*mesh, grid_ptr);
    vtkio.read(grid_ptr);
    apec::logging::info("VTK mesh successfully imported into solver");

    apec::logging::info("Checking orientation of mesh elements");
    libMesh::MeshTools::Modification::orient_elements(*mesh);

    apec::logging::info("Adding physical variables to systems");
    equation_system = std::make_shared<ApecMEquationSystems>(*mesh);
    equation_system->set_parent_object(this);
    equation_system->load_model_data(apec_object);
    equation_system->configure_systems();
    equation_system->initialize_systems();

    apec::logging::info("Preparing to distribute Mesh and Degrees-Of-Freedom across processors");
    auto &main_system = equation_system->get_system<libMesh::System>(
        ApecMSystemBase::system_name_from_type(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    coupling_functor = std::make_unique<AugmentSparsityOnNodes>(*mesh, grid_ptr);
    main_system.get_dof_map().add_coupling_functor(*coupling_functor);

    apec::logging::info("Finding elements' neighbors and distributing mesh across processors");
    mesh->prepare_for_use(true);
    apec::logging::table_root(apec::logging::format_mesh_summary_table(*mesh));

    apec::logging::info("Distributing Degrees-Of-Freedom");
    equation_system->initialize_physics();

    loaded = true;
}

void ApecMModel::solve()
{
    if (!loaded || !mesh || !equation_system || apec_object.is_none())
        apec::logging::error_and_abort("Solver model is not fully initialized. Call load_model() first.");

    apec::logging::info("Starting system solve");

    equation_system->solve_equation_systems(
        py::cast<py::dict>(apec_object.attr("time_step")));

    write_system_to_vtk_grid(
        ApecMSystemBase::system_name_from_type(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    write_system_to_vtk_grid(
        ApecMSystemBase::system_name_from_type(ApecMSystemBase::NODAL_SOLUTION));
    write_system_to_vtk_grid(
        ApecMSystemBase::system_name_from_type(ApecMSystemBase::ELEMENT_SOLUTION));

    apec::logging::info("Solve complete and solution fields written back to VTK grid");
}

void ApecMModel::write_system_to_vtk_grid(const std::string &system_name)
{
    if (!loaded || !mesh || !equation_system || !active_vtk_grid)
        apec::logging::error_and_abort("Solver model is not fully initialized. Call load_model() first.");

    VTKIO_Updated vtkio(*mesh, equation_system.get(), active_vtk_grid);
    vtkio.write_system_grid(system_name);
}

void ApecMModel::clear()
{
    if (closed_)
        return;

    if (equation_system)
        equation_system->clear();

    coupling_functor.reset();
    equation_system.reset();

    if (mesh)
        mesh->clear();

    mesh.reset();

    active_vtk_grid = nullptr;
    apec_object = py::none();
    loaded = false;
    closed_ = true;
}

bool ApecMModel::is_loaded() const
{
    return loaded;
}

py::object ApecMModel::get_apec_object() const
{
    return apec_object;
}

vtkUnstructuredGrid *ApecMModel::get_active_vtk_grid() const
{
    return active_vtk_grid;
}

const AugmentSparsityOnNodes *ApecMModel::get_coupling_functor() const
{
    if (!coupling_functor)
        apec::logging::error_and_abort("Solver coupling functor is not initialized.");
    return coupling_functor.get();
}

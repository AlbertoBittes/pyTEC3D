#include "pytec3d_model.h"

#include "pytec3d_equation_systems.h"
#include "logging.h"
#include "pytec3d_system.h"
#include "augment_sparsity_on_nodes.h"
#include "vtkio.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include <vtkObjectBase.h>
#include <vtkPythonUtil.h>
#include <vtkUnstructuredGrid.h>

namespace py = pybind11;

namespace
{
    py::object active_vtk_grid_from_pytec3d(const py::object &pytec3d_object)
    {
        if (!py::hasattr(pytec3d_object, "mesh"))
            pytec3d::logging::error_and_abort("Expected an pyTEC3D-like Python model object exposing a 'mesh' property.");

        return pytec3d_object.attr("mesh");
    }

    vtkUnstructuredGrid *vtk_grid_from_python(const py::object &vtk_grid)
    {
        if (vtk_grid.is_none())
            pytec3d::logging::error_and_abort("Expected a vtkUnstructuredGrid Python object, received None.");

        vtkObjectBase *vtk_base = vtkPythonUtil::GetPointerFromObject(vtk_grid.ptr(), "vtkUnstructuredGrid");
        if (!vtk_base)
            pytec3d::logging::error_and_abort("Expected a vtkUnstructuredGrid Python object compatible with vtkPythonUtil.");

        auto *grid = vtkUnstructuredGrid::SafeDownCast(vtk_base);
        if (!grid)
            pytec3d::logging::error_and_abort("Python object is not a valid vtkUnstructuredGrid instance.");

        return grid;
    }
}

pyTEC3DModel::~pyTEC3DModel()
{
    clear();
}

void pyTEC3DModel::load_into_solver(const std::shared_ptr<libMesh::LibMeshInit> &libmesh_init,
                                    const std::shared_ptr<const pytec3d::ElementTypeRegistry> &registry,
                                    py::object pytec3dmodel_object)
{
    if (!libmesh_init)
        pytec3d::logging::error_and_abort("Cannot load model without an initialized libMesh runtime.");

    clear();
    closed_ = false;

    pytec3d::logging::info("Loading model into solver runtime");
    pytec3d_object = std::move(pytec3dmodel_object);

    mesh = std::make_unique<libMesh::DistributedMesh>(libmesh_init->comm());
    py::object vtk_grid = active_vtk_grid_from_pytec3d(pytec3d_object);
    vtkUnstructuredGrid *grid_ptr = vtk_grid_from_python(vtk_grid);
    active_vtk_grid = grid_ptr;

    VTKIO vtkio(*mesh, grid_ptr);
    vtkio.read(grid_ptr);
    pytec3d::logging::info("VTK mesh successfully imported into solver");

    pytec3d::logging::info("Checking orientation of mesh elements");
    libMesh::MeshTools::Modification::orient_elements(*mesh);

    pytec3d::logging::info("Adding physical variables to systems");
    equation_system = std::make_shared<pyTEC3DEquationSystems>(*mesh);
    equation_system->set_parent_object(this);
    equation_system->load_model_data(pytec3d_object, registry);
    equation_system->configure_systems();
    equation_system->initialize_systems();

    pytec3d::logging::info("Preparing to distribute Mesh and Degrees-Of-Freedom across processors");
    auto &main_system = equation_system->get_system<libMesh::System>(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM));
    coupling_functor = std::make_unique<AugmentSparsityOnNodes>(*mesh, grid_ptr);
    main_system.get_dof_map().add_coupling_functor(*coupling_functor);

    pytec3d::logging::info("Finding elements' neighbors and distributing mesh across processors");
    mesh->prepare_for_use(true);
    pytec3d::logging::table_root(pytec3d::logging::format_mesh_summary_table(*mesh));

    equation_system->build_rotation_caches(pytec3d_object);

    pytec3d::logging::info("Distributing Degrees-Of-Freedom");
    equation_system->initialize_physics();

    loaded = true;
}

void pyTEC3DModel::solve()
{
    if (!loaded || !mesh || !equation_system || pytec3d_object.is_none())
        pytec3d::logging::error_and_abort("Solver model is not fully initialized. Call load_model() first.");

    pytec3d::logging::info("Starting system solve");

    equation_system->solve_equation_systems(
        py::cast<py::dict>(pytec3d_object.attr("time_step")));

    write_system_to_vtk_grid(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM));
    write_system_to_vtk_grid(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::NODAL_SOLUTION));
    write_system_to_vtk_grid(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::ELEMENT_SOLUTION));

    pytec3d::logging::info("Solve complete and solution fields written back to VTK grid");
}

void pyTEC3DModel::write_system_to_vtk_grid(const std::string &system_name)
{
    if (!loaded || !mesh || !equation_system || !active_vtk_grid)
        pytec3d::logging::error_and_abort("Solver model is not fully initialized. Call load_model() first.");

    VTKIO vtkio(*mesh, equation_system.get(), active_vtk_grid);
    vtkio.write_system_grid(system_name);
}

void pyTEC3DModel::clear()
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
    pytec3d_object = py::none();
    loaded = false;
    closed_ = true;
}

bool pyTEC3DModel::is_loaded() const
{
    return loaded;
}

py::object pyTEC3DModel::get_pytec3d_object() const
{
    return pytec3d_object;
}

vtkUnstructuredGrid *pyTEC3DModel::get_active_vtk_grid() const
{
    return active_vtk_grid;
}

const AugmentSparsityOnNodes *pyTEC3DModel::get_coupling_functor() const
{
    if (!coupling_functor)
        pytec3d::logging::error_and_abort("Solver coupling functor is not initialized.");
    return coupling_functor.get();
}

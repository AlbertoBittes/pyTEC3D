/**
 * @file apec_m_solver.cpp
 * @brief ApecMSolver implementation: runtime lifecycle, model loading, and solve dispatch.
 *
 * The file-local RuntimeManager singleton ensures that libMesh and the spdlog
 * logger are initialised exactly once per process, reference-counted across
 * multiple ApecMSolver instances, and finalised cleanly when shutdown_runtime()
 * is called after all instances have been closed.
 */
#include "apec_m_solver.h"

#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "augment_sparsity_on_nodes.h"
#include "vtkio_updated.h"
#include "apec_m_logging.h"

#include <libmesh/distributed_mesh.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_modification.h>

#include <mpi.h>

#include <vtkUnstructuredGrid.h>
#include <vtkObjectBase.h>
#include <vtkPythonUtil.h>

#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;

namespace
{
    enum class RuntimeState
    {
        Uninitialized,
        Initialized,
        Closed
    };

    class RuntimeManager
    {
    public:
        static RuntimeManager &instance()
        {
            static RuntimeManager manager;
            return manager;
        }

        std::shared_ptr<libMesh::LibMeshInit> acquire(MPI_Comm mpi_comm, int argc, char **argv)
        {
            std::lock_guard<std::mutex> guard(mutex_);

            if (state_ == RuntimeState::Closed)
                apec::logging::error_and_abort("ApecMSolver runtime is permanently closed in this process.");

            if (state_ == RuntimeState::Uninitialized)
            {
                static const char *default_argv[] = {"apec_m_solver"};
                const bool has_user_args = (argc > 0 && argv != nullptr);
                const int init_argc = has_user_args ? argc : 1;
                const char *const *init_argv = has_user_args
                                                   ? const_cast<const char *const *>(argv)
                                                   : default_argv;

                process_libmesh_init_ = std::make_shared<libMesh::LibMeshInit>(init_argc, init_argv, mpi_comm);
                apec::logging::initialize(mpi_comm, "apec_m");
                logging_initialized_ = true;
                state_ = RuntimeState::Initialized;
            }

            ++active_instances_;
            return process_libmesh_init_;
        }

        void release_instance()
        {
            std::lock_guard<std::mutex> guard(mutex_);
            if (active_instances_ > 0)
                --active_instances_;
        }

        void shutdown_runtime()
        {
            std::lock_guard<std::mutex> guard(mutex_);

            if (state_ == RuntimeState::Closed)
                return;

            if (active_instances_ != 0)
                apec::logging::error_and_abort("Cannot shutdown ApecMSolver runtime while solver instances are still active.");

            if (logging_initialized_)
            {
                apec::logging::shutdown();
                logging_initialized_ = false;
            }

            process_libmesh_init_.reset();
            state_ = RuntimeState::Closed;
        }

        bool consume_banner_needed()
        {
            std::lock_guard<std::mutex> guard(mutex_);
            if (banner_printed_)
                return false;
            banner_printed_ = true;
            return true;
        }

    private:
        RuntimeManager() = default;

        std::mutex mutex_;
        RuntimeState state_ = RuntimeState::Uninitialized;
        std::shared_ptr<libMesh::LibMeshInit> process_libmesh_init_;
        int active_instances_ = 0;
        bool logging_initialized_ = false;
        bool banner_printed_ = false;
    };

#define APEC_M_STR_IMPL(x) #x
#define APEC_M_STR(x) APEC_M_STR_IMPL(x)

    const char *apec_m_version_string()
    {
#ifdef APEC_M_VERSION
        return APEC_M_STR(APEC_M_VERSION);
#else
        return "dev";
#endif
    }

    const char *build_mode_string()
    {
#ifdef APEC_M_BUILD_TYPE
        return APEC_M_BUILD_TYPE;
#elif defined(NDEBUG)
        return "Release";
#else
        return "Debug";
#endif
    }

    void print_solver_banner_once()
    {
        if (!RuntimeManager::instance().consume_banner_needed())
            return;

        apec::logging::table_root(apec::logging::format_solver_banner(
            apec_m_version_string(),
            build_mode_string(),
            std::string(__DATE__) + " " + __TIME__));
    }

    py::object active_vtk_grid_from_apec(const py::object &apec_object)
    {
        if (!py::hasattr(apec_object, "mesh"))
            apec::logging::error_and_abort("Expected an APEC-like Python object exposing a 'mesh' property.");

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

    MPI_Comm mpi_comm_from_python(const py::object &mpi_comm)
    {
        if (mpi_comm.is_none())
            return MPI_COMM_WORLD;

        if (!py::hasattr(mpi_comm, "py2f"))
            apec::logging::error_and_abort("Expected an mpi4py communicator with py2f().");

        const auto f_comm = static_cast<MPI_Fint>(py::cast<int>(mpi_comm.attr("py2f")()));
        return MPI_Comm_f2c(f_comm);
    }

} // namespace

void ApecMSolver::ensure_libmesh_initialized(MPI_Comm mpi_comm)
{
    libmesh_init = RuntimeManager::instance().acquire(mpi_comm, argc, argv);
    print_solver_banner_once();
}

ApecMSolver::ApecMSolver(py::object mpi_comm)
{
    mpi_comm_object = std::move(mpi_comm);
    const MPI_Comm active_comm = mpi_comm_from_python(mpi_comm_object);
    ensure_libmesh_initialized(active_comm);
    closed_ = false;
    apec::logging::info("Solver created successfully");
}

ApecMSolver::~ApecMSolver()
{
    close();
}

void ApecMSolver::close()
{
    if (closed_)
        return;

    // Destroy all libMesh and PETSc objects before MPI finalisation.
    // Order matters: equation systems must be cleared before the mesh is reset.
    if (equation_system)
        equation_system->clear();

    coupling_functor.reset();
    equation_system.reset();

    if (mesh)
        mesh->clear();

    mesh.reset();

    active_vtk_grid = nullptr;
    apec_object = py::none();
    mpi_comm_object = py::none();
    loaded = false;
    closed_ = true;

    libmesh_init.reset();
    RuntimeManager::instance().release_instance();
}

void ApecMSolver::shutdown_runtime()
{
    RuntimeManager::instance().shutdown_runtime();
}

void ApecMSolver::load_model(py::object apec_object)
{
    if (closed_)
        apec::logging::error_and_abort("Solver instance has been closed and cannot be reused.");

    apec::logging::info("Loading model into solver");
    this->apec_object = apec_object;

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
    equation_system->load_model_data(this->apec_object);
    equation_system->configure_systems();
    equation_system->initialize_systems();

    apec::logging::info("Preparing to distribute Mesh and Degrees-Of-Freedom across processors");
    auto &main_system = equation_system->get_system<libMesh::System>(
        ApecMSystemBase::system_name_from_type(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    coupling_functor = std::make_unique<AugmentSparsityOnNodes>(*mesh,
                                                                grid_ptr);
    main_system.get_dof_map().add_coupling_functor(*coupling_functor);

    apec::logging::info("Finding elements' neighbors and distributing mesh across processors");
    mesh->prepare_for_use(true);
    apec::logging::table_root(apec::logging::format_mesh_summary_table(*mesh));

    apec::logging::info("Distributing Degrees-Of-Freedom");
    equation_system->initialize_physics();

    loaded = true;
}

void ApecMSolver::solve()
{
    if (!loaded || !mesh)
        apec::logging::error_and_abort("Solver has no loaded model. Call load_model() first.");

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

void ApecMSolver::write_system_to_vtk_grid(const std::string &system_name)
{
    if (!loaded || !mesh || !equation_system || !active_vtk_grid)
        apec::logging::error_and_abort("Solver is not fully initialized. Call load_model() first.");

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

vtkUnstructuredGrid *ApecMSolver::get_active_vtk_grid() const
{
    return active_vtk_grid;
}

const AugmentSparsityOnNodes *ApecMSolver::get_coupling_functor() const
{
    if (!coupling_functor)
        apec::logging::error_and_abort("Solver coupling functor is not initialized.");
    return coupling_functor.get();
}

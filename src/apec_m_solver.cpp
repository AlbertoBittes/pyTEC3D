#include "apec_m_solver.h"

#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "augment_sparsity_on_nodes.h"
#include "element_types_factory.h"
#include "vtkio_updated.h"
#include "apec_m_logging.h"

#include <libmesh/distributed_mesh.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_modification.h>

#include <mpi.h>

#include <vtkUnstructuredGrid.h>

#include <pybind11/stl.h>

#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;

namespace
{
    std::shared_ptr<libMesh::LibMeshInit> process_libmesh_init;
    bool g_solver_banner_printed = false;

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
#if defined(NDEBUG)
        return "Release";
#else
        return "Debug";
#endif
    }

    void print_solver_banner_once()
    {
        if (g_solver_banner_printed)
            return;

        g_solver_banner_printed = true;

        constexpr int inner_width = 76;
        const std::string horizontal_rule = "+" + std::string(inner_width, '-') + "+";
        const auto centered = [](const std::string &text) -> std::string
        {
            constexpr int width = inner_width;
            const std::string clipped = (static_cast<int>(text.size()) > width) ? text.substr(0, width) : text;
            const int pad_total = width - static_cast<int>(clipped.size());
            const int pad_left = pad_total / 2;
            const int pad_right = pad_total - pad_left;
            return "|" + std::string(pad_left, ' ') + clipped + std::string(pad_right, ' ') + "|";
        };
        const auto kv = [](const std::string &key, const std::string &value) -> std::string
        {
            constexpr int width = inner_width;
            constexpr int key_width = 20;
            std::ostringstream out;
            out << "| " << std::left << std::setw(key_width) << key << " : ";
            std::string row = out.str() + value;
            if (static_cast<int>(row.size()) > width + 1)
                row = row.substr(0, width + 1);
            if (static_cast<int>(row.size()) < width + 1)
                row += std::string(width + 1 - row.size(), ' ');
            row += "|";
            return row;
        };

        std::ostringstream banner;
        banner << horizontal_rule << "\n";
        banner << centered("A P E C - M  S O L V E R") << "\n";
        banner << centered("Low-Frequency Electromagnetics Finite Elements solver") << "\n";
        banner << horizontal_rule << "\n";
        banner << kv("Version", apec_m_version_string()) << "\n";
        banner << kv("Build mode", build_mode_string()) << "\n";
        banner << kv("Build stamp", std::string(__DATE__) + " " + __TIME__) << "\n";
        banner << horizontal_rule;
        apec::logging::table_root(banner.str());
    }

    std::string join_elem_dimensions(const std::set<unsigned char> &dims)
    {
        std::ostringstream out;
        bool first = true;
        for (const auto dim : dims)
        {
            if (!first)
                out << ",";
            out << static_cast<unsigned int>(dim);
            first = false;
        }
        return out.str();
    }

    std::string join_elem_orders(const std::set<Order> &orders)
    {
        std::ostringstream out;
        bool first = true;
        for (const auto order : orders)
        {
            if (!first)
                out << ",";
            out << Utility::enum_to_string<Order>(order);
            first = false;
        }
        return out.str();
    }

    std::string join_counts_by_rank(const std::vector<dof_id_type> &counts)
    {
        std::ostringstream out;
        for (std::size_t rank = 0; rank < counts.size(); ++rank)
        {
            if (rank != 0)
                out << ", ";
            out << "p" << rank << "=" << counts[rank];
        }
        return out.str();
    }

    std::string format_mesh_summary_table(const libMesh::MeshBase &mesh)
    {
        constexpr int inner_width = 76;
        const std::string horizontal_rule = "+" + std::string(inner_width, '-') + "+";
        const auto kv = [](const std::string &key, const std::string &value) -> std::string
        {
            constexpr int width = inner_width;
            constexpr int key_width = 26;
            std::ostringstream out;
            out << "| " << std::left << std::setw(key_width) << key << " : ";
            std::string row = out.str() + value;
            if (static_cast<int>(row.size()) > width + 1)
                row = row.substr(0, width + 1);
            if (static_cast<int>(row.size()) < width + 1)
                row += std::string(width + 1 - row.size(), ' ');
            row += "|";
            return row;
        };

        std::vector<dof_id_type> local_nodes_by_rank;
        std::vector<dof_id_type> local_elems_by_rank;
        mesh.comm().allgather(mesh.n_local_nodes(), local_nodes_by_rank);
        mesh.comm().allgather(mesh.n_local_elem(), local_elems_by_rank);

        std::ostringstream out;
        out << horizontal_rule << "\n";
        out << "| " << std::left << std::setw(inner_width) << "Mesh Summary" << "|\n";
        out << horizontal_rule << "\n";
        out << kv("elem dimensions", "{" + join_elem_dimensions(mesh.elem_dimensions()) + "}") << "\n";
        out << kv("elem default orders", "{" + join_elem_orders(mesh.elem_default_orders()) + "}") << "\n";
        out << kv("supported nodal order", Utility::enum_to_string<Order>(mesh.supported_nodal_order())) << "\n";
        out << kv("spatial dimension", std::to_string(mesh.spatial_dimension())) << "\n";
        out << kv("total nodes", std::to_string(mesh.n_nodes())) << "\n";
        out << kv("total elements", std::to_string(mesh.n_elem())) << "\n";
        out << kv("nodes by rank", join_counts_by_rank(local_nodes_by_rank)) << "\n";
        out << kv("elements by rank", join_counts_by_rank(local_elems_by_rank)) << "\n";
        out << horizontal_rule;
        return out.str();
    }

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
    {
        print_solver_banner_once();
        return;
    }

    if (process_libmesh_init)
    {
        libmesh_init = process_libmesh_init;
        apec::logging::initialize(mpi_comm, "apec_m");
        print_solver_banner_once();
        return;
    }

    // libMesh expects a valid argv array. Passing nullptr may segfault in some builds.
    static const char *default_argv[] = {"apec_m_solver"};
    const bool has_user_args = (argc > 0 && argv != nullptr);
    const int init_argc = has_user_args ? argc : 1;
    const char *const *init_argv = has_user_args
                                       ? const_cast<const char *const *>(argv)
                                       : default_argv;

    process_libmesh_init = std::make_shared<libMesh::LibMeshInit>(init_argc, init_argv, mpi_comm);
    libmesh_init = process_libmesh_init;
    apec::logging::initialize(mpi_comm, "apec_m");
    print_solver_banner_once();
}

ApecMSolver::ApecMSolver(py::object mpi_comm)
{
    mpi_comm_object = std::move(mpi_comm);
    const MPI_Comm active_comm = mpi_comm_from_python(mpi_comm_object);
    ensure_libmesh_initialized(active_comm);
    apec::logging::info("Solver created successfully");
}

ApecMSolver::~ApecMSolver()
{
    close();
}

void ApecMSolver::close()
{
    // Ensure libMesh/PETSc objects are destroyed before MPI finalization.
    if (equation_system)
        equation_system->clear();

    coupling_functor.reset();
    equation_system.reset();

    if (mesh)
        mesh->clear();

    element_types_factory.reset();
    mesh.reset();

    active_vtk_grid = nullptr;
    apec_object = py::none();
    mpi_comm_object = py::none();
    loaded = false;

    libmesh_init.reset();
    apec::logging::shutdown();

    if (process_libmesh_init && process_libmesh_init.use_count() == 1)
        process_libmesh_init.reset();
}

void ApecMSolver::shutdown_runtime()
{
    apec::logging::shutdown();
    if (process_libmesh_init && process_libmesh_init.use_count() == 1)
        process_libmesh_init.reset();
}

void ApecMSolver::load_model(py::object apec_object)
{
    apec::logging::info("Loading model into solver");
    this->apec_object = apec_object;
    element_types_factory = std::make_shared<ElementTypesFactory>(apec_object.attr("element_types"));

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
    equation_system->configure_systems();
    equation_system->initialize_systems();

    apec::logging::info("Preparing to distribute Mesh and Degrees-Of-Freedom across processors");
    apec::logging::info("Identifying coupled nodes for mesh distribution");
    auto &main_system = equation_system->get_system<ApecMLinearSystem>(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    coupling_functor = std::make_unique<AugmentSparsityOnNodes>(*mesh,
                                                                grid_ptr);
    main_system.get_dof_map().add_coupling_functor(*coupling_functor);

    apec::logging::info("Finding elements' neighbors and distributing mesh across processors");
    mesh->prepare_for_use(true);
    apec::logging::table_root(format_mesh_summary_table(*mesh));

    apec::logging::info("Distributing Degrees-Of-Freedom");
    equation_system->initialize_physics();

    loaded = true;
}

void ApecMSolver::solve()
{
    if (!loaded || !mesh)
        libmesh_error_msg("Solver has no loaded model. Call load_model() first.");

    apec::logging::info("Starting system solve");

    equation_system->solve_equation_systems(
        py::cast<py::dict>(apec_object.attr("time_step")));

    write_system_to_vtk_grid(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    write_system_to_vtk_grid(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::NODAL_SOLUTION));
    write_system_to_vtk_grid(
        ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::ELEMENT_SOLUTION));

    apec::logging::info("Solve complete and solution fields written back to VTK grid");
}

void ApecMSolver::write_system_to_vtk_grid(const std::string &system_name)
{
    if (!loaded || !mesh || !equation_system || !active_vtk_grid)
        libmesh_error_msg("Solver is not fully initialized. Call load_model() first.");

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

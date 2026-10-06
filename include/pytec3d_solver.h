/**
 * @file pytec3d_solver.h
 * @brief Top-level solver class and pybind11-facing API for the pyTEC3D engine.
 */
#ifndef pyTEC3D_SOLVER_H
#define pyTEC3D_SOLVER_H

#include <memory>
#include <string>

#include <pybind11/pybind11.h>

class pyTEC3DModel;

namespace pytec3d
{
    class ElementTypeParser;
}

namespace libMesh
{
    class LibMeshInit;
}

class pyTEC3DSolver
{
public:
    /**
     * @brief Constructs the solver and initialises the libMesh / MPI runtime.
     * @param num_threads_per_mpi_process Number of threads to use per MPI rank.
     */
    explicit pyTEC3DSolver(int num_threads_per_mpi_process = 1);

    /** @brief Calls close() to free all libMesh and PETSc resources. */
    ~pyTEC3DSolver();

    /**
     * @brief Imports the model from a Python pyTEC3D model object, distributes mesh and DOFs.
     * @param solver_settings dict with the PETSc settings of the Python solver, read at every solve.
     */
    void load_model(pybind11::object pytec3d_object, pybind11::object solver_settings);

    /** @brief Releases all libMesh, PETSc, and VTK resources.  Idempotent. */
    void close();

    /**
     * @brief Shuts down the process-level libMesh / spdlog runtime.
     *
     * Must be called after all pyTEC3DSolver instances have been closed. After
     * shutdown the runtime cannot be re-initialised in the same process.
     */
    static void shutdown_runtime();

    /** @brief Executes one time-step solve (dict with time, theta, transient_effects) and writes results back to the VTK grid. */
    void solve(const pybind11::dict &time_step);

    /** @brief Localises a named system's solution onto the active VTK grid. */
    void write_system_to_vtk_grid(const std::string &system_name);

    /** @brief Returns true if load_model() has completed successfully. */
    bool is_loaded() const;

    /** @brief The registries compiled into the module, as a dict (domains, fields, element types, properties). */
    static pybind11::dict registries();

private:
    std::shared_ptr<const pytec3d::ElementTypeParser> element_type_parser_;
    int argc = 0;
    char **argv = nullptr;
    std::shared_ptr<libMesh::LibMeshInit> libmesh_init;
    int num_threads_per_mpi_process_ = 1;
    std::shared_ptr<pyTEC3DModel> model_;
    bool loaded = false;
    bool closed_ = false;

    void ensure_libmesh_initialized();
};

#endif // pyTEC3D_SOLVER_H

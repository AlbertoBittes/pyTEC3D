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
    class ElementTypeRegistry;
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
     */
    void load_model(pybind11::object pytec3d_object);

    /** @brief Releases all libMesh, PETSc, and VTK resources.  Idempotent. */
    void close();

    /**
     * @brief Shuts down the process-level libMesh / spdlog runtime.
     *
     * Must be called after all pyTEC3DSolver instances have been closed. After
     * shutdown the runtime cannot be re-initialised in the same process.
     */
    static void shutdown_runtime();

    /** @brief Executes one time-step solve and writes results back to the VTK grid. */
    void solve();

    /** @brief Localises a named system's solution onto the active VTK grid. */
    void write_system_to_vtk_grid(const std::string &system_name);

    /** @brief Returns true if load_model() has completed successfully. */
    bool is_loaded() const;

    /** @brief Base element types, fields and domains loaded from the JSON definitions. */
    const pytec3d::ElementTypeRegistry &element_type_registry() const { return *element_type_registry_; }

private:
    std::shared_ptr<const pytec3d::ElementTypeRegistry> element_type_registry_;
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

/**
 * @file apec_m_solver.h
 * @brief Top-level solver class and pybind11-facing API for the APEC-M engine.
 */
#ifndef APEC_M_SOLVER_H
#define APEC_M_SOLVER_H

#include <memory>
#include <string>

#include <pybind11/pybind11.h>

#include "apec_m_visibility.h"

class ApecMModel;

namespace libMesh
{
    class LibMeshInit;
}

class APEC_M_HIDDEN ApecMSolver
{
public:
    /**
     * @brief Constructs the solver and initialises the libMesh / MPI runtime.
     * @param num_threads_per_mpi_process Number of threads to use per MPI rank.
     */
    explicit ApecMSolver(int num_threads_per_mpi_process = 1);

    /** @brief Calls close() to free all libMesh and PETSc resources. */
    ~ApecMSolver();

    /**
     * @brief Imports the model from a Python APEC model object, distributes mesh and DOFs.
     */
    void load_model(pybind11::object apec_object);

    /** @brief Releases all libMesh, PETSc, and VTK resources.  Idempotent. */
    void close();

    /**
     * @brief Shuts down the process-level libMesh / spdlog runtime.
     *
     * Must be called after all ApecMSolver instances have been closed. After
     * shutdown the runtime cannot be re-initialised in the same process.
     */
    static void shutdown_runtime();

    /** @brief Executes one time-step solve and writes results back to the VTK grid. */
    void solve();

    /** @brief Localises a named system's solution onto the active VTK grid. */
    void write_system_to_vtk_grid(const std::string &system_name);

    /** @brief Returns true if load_model() has completed successfully. */
    bool is_loaded() const;

private:
    int argc = 0;
    char **argv = nullptr;
    std::shared_ptr<libMesh::LibMeshInit> libmesh_init;
    int num_threads_per_mpi_process_ = 1;
    std::shared_ptr<ApecMModel> model_;
    bool loaded = false;
    bool closed_ = false;

    void ensure_libmesh_initialized();
};

#endif // APEC_M_SOLVER_H

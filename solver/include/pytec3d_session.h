/**
 * @file pytec3d_session.h
 * @brief A solver session: the process-wide libMesh / MPI runtime and one model loaded into it.
 *
 * The session is plain C++, with no Python or VTK dependency: the model is handed to it as views of the caller's
 * arrays and the model data as JSON (model_input.h), the time step and the solver settings with each solve. A Python
 * API (python_bindings.cpp) or any C++ program can drive it.
 */
#pragma once

#include <memory>
#include <string>

#include "model_input.h"
#include "solver_configuration.h"
#include "time_step.h"

class pyTEC3DModel;

namespace libMesh
{
    class LibMeshInit;
}

class pyTEC3DSession
{
public:
    /**
     * @brief Creates a session, initialising the libMesh / MPI runtime on first use in the process.
     * @param num_threads_per_mpi_process Threads per MPI rank; only the first session of the process sets it.
     */
    explicit pyTEC3DSession(int num_threads_per_mpi_process = 1);

    /** @brief Calls close(). */
    ~pyTEC3DSession();

    pyTEC3DSession(const pyTEC3DSession &) = delete;
    pyTEC3DSession &operator=(const pyTEC3DSession &) = delete;

    /**
     * @brief Builds the mesh from the views of @p input, reads the model data, distributes mesh and DOFs.
     *
     * The session keeps the array views (not the arrays: they must stay alive, see model_input.h) and copies the model
     * data. The mesh, the TYPE / MATERIAL / CONSTANT / COORDINATE_SYSTEM ids, the CP_<k> couplings and the rotations
     * are read here: changing them later takes a new load_model().
     */
    void load_model(const pytec3d::ModelInput &input);

    /**
     * @brief Replaces the model arrays of the loaded model by @p arrays (same sizes and fields): to hand the solver
     * arrays the caller has reallocated. Loads, Dirichlet values and outputs are read from them at the next solve.
     */
    void set_arrays(const pytec3d::ModelArrays &arrays);

    /** @brief Releases the model and its libMesh and PETSc resources; the runtime stays up. Idempotent. */
    void close();

    /**
     * @brief Shuts down the process-wide libMesh / spdlog runtime.
     *
     * Must be called after all sessions have been closed. After shutdown the runtime cannot be re-initialised
     * in the same process (MPI cannot be initialised twice).
     */
    static void shutdown_runtime();

    /**
     * @brief Solves one time step with @p settings and writes the solution into the output arrays.
     *
     * The settings may differ from one solve to the next; a change of the linear solver settings refactorizes the
     * system matrix.
     */
    void solve(const TimeStepInput &time_step, const SolverSettings &settings);

    /** @brief Writes the solution of the system @p system_name into the output arrays. */
    void write_system_output(const std::string &system_name);

    /** @brief What the last solve() did. */
    const StepReport &last_step_report() const;

    /** @brief True once load_model() has completed and until close(). */
    bool is_loaded() const;

private:
    void require_loaded() const;

    std::shared_ptr<libMesh::LibMeshInit> libmesh_init_; ///< Shared by all sessions of the process.
    std::unique_ptr<pyTEC3DModel> model_;
    bool closed_ = false;
};

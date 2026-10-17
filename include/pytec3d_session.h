/**
 * @file pytec3d_session.h
 * @brief A solver session: the process-wide libMesh / MPI runtime and one model loaded into it.
 *
 * The session is plain C++: python_bindings.cpp converts the Python objects (the model, its grid, the
 * time_step and solver_settings dicts) into the types below and exposes the session to Python.
 */
#pragma once

#include <memory>
#include <string>

#include "coordinate_rotations.h"
#include "solver_configuration.h"
#include "time_step.h"

class pyTEC3DModel;
class vtkUnstructuredGrid;

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
     * @brief Imports the grid and the model data, distributes mesh and DOFs.
     *
     * The session keeps a reference to @p grid (it writes the solution into it) and copies the model data and the
     * solver settings: later changes to them only take effect when the model is loaded again.
     * @param model_json The model data as JSON text (model_data.h).
     * @param rotations Local-to-global rotations of the coordinate systems; only used during the call.
     */
    void load_model(vtkUnstructuredGrid *grid,
                    const std::string &model_json,
                    const SolverSettings &solver_settings,
                    const CoordinateRotations &rotations);

    /** @brief Releases the model and its libMesh, PETSc and VTK resources; the runtime stays up. Idempotent. */
    void close();

    /**
     * @brief Shuts down the process-wide libMesh / spdlog runtime.
     *
     * Must be called after all sessions have been closed. After shutdown the runtime cannot be re-initialised
     * in the same process (MPI cannot be initialised twice).
     */
    static void shutdown_runtime();

    /** @brief Solves one time step and writes the solution into the grid. */
    void solve(const TimeStepInput &time_step);

    /** @brief Writes the solution of the system @p system_name into the grid. */
    void write_system_to_vtk_grid(const std::string &system_name);

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

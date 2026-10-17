/**
 * @file pytec3d_model.h
 * @brief A model loaded into the solver: the grid, the libMesh mesh and the equation systems built from them.
 */
#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include "coordinate_rotations.h"
#include "solver_configuration.h"
#include "time_step.h"

class AugmentSparsityOnNodes;
class pyTEC3DEquationSystems;

namespace libMesh
{
    class DistributedMesh;
    namespace Parallel
    {
        class Communicator;
    }
}

class pyTEC3DModel
{
public:
    pyTEC3DModel();
    ~pyTEC3DModel();

    pyTEC3DModel(const pyTEC3DModel &) = delete;
    pyTEC3DModel &operator=(const pyTEC3DModel &) = delete;

    /**
     * @brief Imports @p grid and the model data, distributes mesh and DOFs.
     *
     * Keeps a reference to @p grid: the solution is written into it after each solve.
     */
    void load(const libMesh::Parallel::Communicator &comm,
              vtkUnstructuredGrid *grid,
              const std::string &model_json,
              const SolverSettings &solver_settings,
              const CoordinateRotations &rotations);

    /** @brief Solves one time step and writes the solution of every system into the grid. */
    void solve(const TimeStepInput &time_step);

    void write_system_to_vtk_grid(const std::string &system_name);

    /** @brief What the last solve() did. */
    const StepReport &last_step_report() const;

    /** @brief Releases the equation systems and the mesh. */
    void clear();

    bool is_loaded() const;
    vtkUnstructuredGrid *get_active_vtk_grid() const;
    const AugmentSparsityOnNodes *get_coupling_functor() const;

private:
    void require_loaded() const;

    // Destroyed in reverse order: the equation systems before the coupling functor they use, the mesh last.
    vtkSmartPointer<vtkUnstructuredGrid> grid_;
    std::unique_ptr<libMesh::DistributedMesh> mesh_;
    std::unique_ptr<AugmentSparsityOnNodes> coupling_functor_;
    std::unique_ptr<pyTEC3DEquationSystems> equation_systems_;
    std::size_t flipped_elements_ = 0; ///< Reoriented when the mesh was imported.
    bool loaded_ = false;
};

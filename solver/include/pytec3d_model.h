/**
 * @file pytec3d_model.h
 * @brief A model loaded into the solver: the views of the caller's arrays, the libMesh mesh and the equation systems.
 */
#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "array_io.h"
#include "model_input.h"
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
     * @brief Builds the mesh from the views of @p input, reads the model data, distributes mesh and DOFs.
     *
     * Keeps the array views of @p input: the loads are read from them and the solution written into them.
     */
    void load(const libMesh::Parallel::Communicator &comm, const pytec3d::ModelInput &input);

    /** @brief Replaces the array views (same sizes as at the load). */
    void set_arrays(const pytec3d::ModelArrays &arrays);

    /** @brief Solves one time step and writes the solution of every system into the output arrays. */
    void solve(const TimeStepInput &time_step, const SolverSettings &settings);

    /** @brief Writes the solution of the system @p system_name into the output arrays. */
    void write_system_output(const std::string &system_name);

    /** @brief What the last solve() did. */
    const StepReport &last_step_report() const;

    /** @brief Releases the equation systems and the mesh. */
    void clear();

    bool is_loaded() const;

    /** @brief The model arrays (model_input.h). */
    const pytec3d::ModelArrays &arrays() const { return arrays_; }

    const AugmentSparsityOnNodes *get_coupling_functor() const;

private:
    void require_loaded() const;
    /** @brief Aborts when the arrays miss an output field of the systems (a field name the systems and the
     * element types do not agree on). */
    void check_outputs() const;

    // Destroyed in reverse order: the equation systems before the coupling functor they use, the mesh last.
    pytec3d::ModelArrays arrays_;
    std::vector<pytec3d::FieldArrayRequirement> required_fields_; ///< Of the element types of the model data.
    std::size_t n_points_ = 0;
    std::size_t n_cells_ = 0;
    std::unique_ptr<libMesh::DistributedMesh> mesh_;
    std::unique_ptr<AugmentSparsityOnNodes> coupling_functor_;
    std::unique_ptr<pyTEC3DEquationSystems> equation_systems_;
    std::size_t flipped_elements_ = 0; ///< Reoriented when the mesh was imported.
    bool loaded_ = false;
};

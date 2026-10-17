/**
 * @file pytec3d_equation_systems.h
 * @brief pyTEC3D extension of libMesh::EquationSystems.
 *
 * Owns the full set of solution systems (DOF Solution, Nodal Solution, Element
 * Solution), orchestrates variable initialisation, physics initialisation, and
 * the per-time-step solve sequence.
 */
#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include <libmesh/equation_systems.h>
#include <libmesh/tensor_value.h>

#include "coordinate_rotations.h"
#include "model_data.h"
#include "solver_configuration.h"
#include "time_step.h"

class ModelElementTypes;
class pyTEC3DSystemBase;
class pyTEC3DModel;
class pyTEC3DDOFSystemBase;

/**
 * @brief Manages all libMesh solution systems for one pyTEC3D simulation.
 *
 * Derives from libMesh::EquationSystems; belongs to a pyTEC3DModel (model()).
 * Responsible for:
 *   - Registering and configuring DOF/nodal/element solution systems.
 *   - Driving variable and physics initialisation (variable adds, mesh prepare, DOF distribute).
 *   - Executing the per-time-step solve via solve_equation_systems().
 */
class pyTEC3DEquationSystems
    : public libMesh::EquationSystems

{
public:
    using RotationMap = std::unordered_map<libMesh::dof_id_type, libMesh::TensorValue<libMesh::Real>>;

    pyTEC3DEquationSystems(libMesh::MeshBase &mesh, pyTEC3DModel &model);

    /** @brief The model the systems belong to: its grid and coupling functor. */
    pyTEC3DModel &model() const { return model_; }
    ~pyTEC3DEquationSystems() override;

    /** @brief Parses the model data JSON (model.data) and builds the element types it selects. */
    void load_model_data(const std::string &model_json, const SolverSettings &solver_settings);

    /**
     * @brief Dry run of the coefficients on one element per (TYPE, MATERIAL, CONSTANT) combination of the mesh.
     *
     * Aborts, naming the element type, domain and material or constants, when a term active in the domain
     * needs a material property or real constant the model data leaves undefined.
     */
    void validate_model_data() const;

    /** @brief The model data of the loaded model. */
    const pytec3d::ModelData &get_model_data() const;

    /** @brief The linear solver settings, fixed when the model is loaded. */
    const SolverSettings &get_solver_settings() const;

    /** @brief Updates time-step state and decides whether this solve advances time. */
    void prepare_time_step_solve(const TimeStepInput &time_step);

    /** @brief Executes a single time-step solve: prepares constraints, assembles, solves, and updates post-processing fields. */
    void solve_equation_systems(const TimeStepInput &time_step);

    /** @brief Registers each system type with the EquationSystems container. */
    void configure_systems();

    /** @brief Adds physical variables to each registered system. */
    void initialize_systems();

    /** @brief Performs any post-init work that requires variables to already exist. */
    void post_initialize_systems();

    /** @brief Calls libMesh::EquationSystems::init() to distribute DOFs across the mesh. */
    void initialize_physics();

    /** @brief Returns true when the DOF system is configured as nonlinear. */
    bool uses_nonlinear_system() const;

    /** @brief Returns the active DOF solution system (linear or nonlinear). */
    pyTEC3DDOFSystemBase &get_dof_solution_system();

    /** @brief Returns the active DOF solution system (linear or nonlinear). */
    const pyTEC3DDOFSystemBase &get_dof_solution_system() const;

    /** @brief What the current (or last) time-step solve did; filled by the systems. */
    StepReport &step_report() { return step_report_; }
    const StepReport &step_report() const { return step_report_; }

    /** @brief Returns the current time-step state. */
    const TimeStep *current_time_step() const;

    /** @brief Returns the previous (already accepted) time-step state. */
    const TimeStep *previous_time_step() const;

    /** @brief True when a law is selected on an element of the grid: the DOF system is solved with Newton. */
    bool has_nonlinear_elements() const;

    /** @brief Highest time order of the element types of the model data: 0 static, 1, or 2. */
    unsigned int time_order() const;

    /** @brief Scheme of the current time step for the time order of the model (time_step.h). */
    TimeScheme time_scheme() const;

    /** @brief The element types of the model, by TYPE id; aborts before load_model_data(). */
    const ModelElementTypes &element_types() const;

    /** @brief Returns precomputed element-center local-to-global rotations keyed by element id. */
    const RotationMap &get_element_local_to_global_rotations() const;

    /** @brief Returns precomputed node local-to-global rotations keyed by node id. */
    const RotationMap &get_node_local_to_global_rotations() const;

    /** @brief Precomputes the local-to-global rotations of the local elements (at their centres) and of the nodes. */
    void build_rotation_caches(const CoordinateRotations &rotations);

private:
    /** @brief The largest input of each nonlinear law per id, into the step report; warns on extrapolation. */
    void compute_law_ranges();

    pyTEC3DModel &model_;
    bool use_nonlinear_system_ = false;
    std::unique_ptr<TimeStep> time_step_;
    StepReport step_report_;
    std::unique_ptr<TimeStep> old_time_step_;
    std::unique_ptr<ModelElementTypes> element_types_;
    pytec3d::ModelData model_data_;
    SolverSettings solver_settings_;
    RotationMap element_local_to_global_rotations_;
    RotationMap node_local_to_global_rotations_;
};


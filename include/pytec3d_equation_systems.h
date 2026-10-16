/**
 * @file pytec3d_equation_systems.h
 * @brief pyTEC3D extension of libMesh::EquationSystems.
 *
 * Owns the full set of solution systems (DOF Solution, Nodal Solution, Element
 * Solution), orchestrates variable initialisation, physics initialisation, and
 * the per-time-step solve sequence.
 */
#ifndef pyTEC3D_EQUATION_SYSTEMS_H
#define pyTEC3D_EQUATION_SYSTEMS_H

#include <unordered_map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <array>

#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/tensor_value.h>
#include <pybind11/pybind11.h>

#include "fe_parent_object.h"
#include "solver_configuration.h"
#include "time_step.h"
#include "model_data.h"

class ElementTypesFactory;
class pyTEC3DSystemBase;
class pyTEC3DModel;
class pyTEC3DDOFSystemBase;

/**
 * @brief Manages all libMesh solution systems for one pyTEC3D simulation.
 *
 * Derives from libMesh::EquationSystems and FEParentObjectBase<pyTEC3DModel>.
 * Responsible for:
 *   - Registering and configuring DOF/nodal/element solution systems.
 *   - Driving variable and physics initialisation (variable adds, mesh prepare, DOF distribute).
 *   - Executing the per-time-step solve via solve_equation_systems().
 */
class pyTEC3DEquationSystems
    : public libMesh::EquationSystems,
      public FEParentObjectBase<pyTEC3DModel>

{
public:
    using RotationMap = std::unordered_map<libMesh::dof_id_type, libMesh::TensorValue<libMesh::Real>>;

    /**
     * @brief Parses the model data JSON (model.data) and builds the element types it selects.
     * @p solver_settings is the dict of the Python solver, read at every solve.
     */
    void load_model_data(const std::string &model_json, pybind11::object solver_settings);

    /**
     * @brief Dry run of the coefficients on one element per (TYPE, MATERIAL, CONSTANT) combination of the mesh.
     *
     * Aborts, naming the element type, domain and material or constants, when a term active in the domain
     * needs a material property or real constant the model data leaves undefined.
     */
    void validate_model_data() const;

    /** @brief The model data of the loaded model. */
    const pytec3d::ModelData &get_model_data() const;

    /** @brief The solver_settings dict of the Python solver (none when it passed none). */
    const pybind11::object &get_solver_settings() const;

    /** @brief Updates time-step state and decides whether this solve advances time. */
    void prepare_time_step_solve(const pybind11::dict &time_step);

    /** @brief Executes a single time-step solve: prepares constraints, assembles, solves, and updates post-processing fields. */
    void solve_equation_systems(const pybind11::dict &time_step);

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

    /** @brief Returns the element-types factory built from the model description. */
    const std::shared_ptr<ElementTypesFactory> &get_element_types_factory() const;

    /** @brief Returns precomputed element-center local-to-global rotations keyed by element id. */
    const RotationMap &get_element_local_to_global_rotations() const;

    /** @brief Returns precomputed node local-to-global rotations keyed by node id. */
    const RotationMap &get_node_local_to_global_rotations() const;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    pyTEC3DEquationSystems(Args &&...args)
        : libMesh::EquationSystems(std::forward<Args>(args)...)
    {
    }

    void build_rotation_caches(pybind11::object pytec3dmodel_object);

private:
    bool use_nonlinear_system_ = false;
    std::unique_ptr<TimeStep> time_step_;
    std::unique_ptr<TimeStep> old_time_step_;
    std::shared_ptr<ElementTypesFactory> element_types_factory_;
    pytec3d::ModelData model_data_;
    pybind11::object solver_settings_;
    RotationMap element_local_to_global_rotations_;
    RotationMap node_local_to_global_rotations_;
};

#endif // pyTEC3D_EQUATION_SYSTEMS_H

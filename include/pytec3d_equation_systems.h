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
#include "bh_curve.h"

class ElementTypesFactory;
class pyTEC3DSystemBase;
class pyTEC3DModel;
class pyTEC3DDOFSystemBase;

namespace pytec3d
{
    class ElementTypeParser;
}

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
    using PropertyRowsMap = std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::Real>>;
    using RotationMap = std::unordered_map<libMesh::dof_id_type, libMesh::TensorValue<libMesh::Real>>;

    /** @brief Parses and caches model metadata shared across systems/builders; @p parser must outlive the element types. */
    void load_model_data(pybind11::object pytec3d_object,
                         const std::shared_ptr<const pytec3d::ElementTypeParser> &parser);

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

    /** @brief Returns true when material BH curves are available. */
    bool has_bh_curves() const;

    /** @brief Returns parsed BH curves indexed by material and component. */
    const BHCurveMaterialsMap &get_bh_curves() const;

    /** @brief Returns the element-types factory built from the model description. */
    const std::shared_ptr<ElementTypesFactory> &get_element_types_factory() const;

    /** @brief Returns cached material rows indexed by material ID. */
    const PropertyRowsMap &get_material_rows() const;

    /** @brief Returns cached real-constant rows indexed by constants ID. */
    const PropertyRowsMap &get_constant_rows() const;

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
    BHCurveMaterialsMap bh_curves_;
    PropertyRowsMap materials_rows_;
    PropertyRowsMap constants_rows_;
    RotationMap element_local_to_global_rotations_;
    RotationMap node_local_to_global_rotations_;
};

#endif // pyTEC3D_EQUATION_SYSTEMS_H

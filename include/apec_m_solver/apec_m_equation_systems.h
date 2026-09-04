/**
 * @file apec_m_equation_systems.h
 * @brief APEC-M extension of libMesh::EquationSystems.
 *
 * Owns the full set of solution systems (DOF Solution, Nodal Solution, Element
 * Solution), orchestrates variable initialisation, physics initialisation, and
 * the per-time-step solve sequence.
 */
#ifndef APEC_M_EQUATION_SYSTEMS_H
#define APEC_M_EQUATION_SYSTEMS_H

#include <unordered_map>
#include <memory>
#include <tuple>
#include <vector>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <pybind11/pybind11.h>

#include "apec_m_solver/fe_parent_object.h"
#include "apec_m_solver/solver_configuration.h"
#include "apec_m_solver/apec_m_visibility.h"
#include "apec_m_solver/apec_m_time_step.h"
#include "apec_m_solver/bh_curve.h"

class ElementTypesFactory;
class ApecMSystemBase;
class ApecMModel;
class ApecMDOFSystemBase;
namespace apec_m
{
    class CoordinateSystem;
}

/**
 * @brief Manages all libMesh solution systems for one APEC-M simulation.
 *
 * Derives from libMesh::EquationSystems and FEParentObjectBase<ApecMModel>.
 * Responsible for:
 *   - Registering and configuring DOF/nodal/element solution systems.
 *   - Driving variable and physics initialisation (variable adds, mesh prepare, DOF distribute).
 *   - Executing the per-time-step solve via solve_equation_systems().
 */
class APEC_M_HIDDEN ApecMEquationSystems
    : public libMesh::EquationSystems,
      public FEParentObjectBase<ApecMModel>

{
public:
    using PropertyRowsMap = std::unordered_map<libMesh::dof_id_type, std::vector<libMesh::Real>>;
    using CoordinateSystemsMap = std::unordered_map<libMesh::dof_id_type, std::shared_ptr<apec_m::CoordinateSystem>>;

    /** @brief Parses and caches model metadata shared across systems/builders. */
    void load_model_data(pybind11::object apec_object);

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
    ApecMDOFSystemBase &get_dof_solution_system();

    /** @brief Returns the active DOF solution system (linear or nonlinear). */
    const ApecMDOFSystemBase &get_dof_solution_system() const;

    /** @brief Returns the current time-step state. */
    const ApecMTimeStep *current_time_step() const;

    /** @brief Returns the previous (already accepted) time-step state. */
    const ApecMTimeStep *previous_time_step() const;

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

    /** @brief Returns cached coordinate systems indexed by positive ID. */
    const CoordinateSystemsMap &get_coordinate_systems() const;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : libMesh::EquationSystems(std::forward<Args>(args)...)
    {
    }

private:
    bool use_nonlinear_system_ = false;
    std::unique_ptr<ApecMTimeStep> time_step_;
    std::unique_ptr<ApecMTimeStep> old_time_step_;
    std::shared_ptr<ElementTypesFactory> element_types_factory_;
    BHCurveMaterialsMap bh_curves_;
    PropertyRowsMap materials_rows_;
    PropertyRowsMap constants_rows_;
    CoordinateSystemsMap coordinate_systems_;
};

#endif // APEC_M_EQUATION_SYSTEMS_H

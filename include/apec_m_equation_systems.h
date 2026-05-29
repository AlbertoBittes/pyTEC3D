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

#include <memory>
#include <tuple>
#include <vector>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <pybind11/pybind11.h>

#include "fe_parent_object.h"
#include "solver_configuration.h"
#include "apec_m_visibility.h"

class ElementTypesFactory;
class ApecMSystemBase;
class ApecMSolver;

/**
 * @brief Manages all libMesh solution systems for one APEC-M simulation.
 *
 * Derives from libMesh::EquationSystems and FEParentObjectBase<ApecMSolver>.
 * Responsible for:
 *   - Registering and configuring DOF/nodal/element solution systems.
 *   - Driving variable and physics initialisation (variable adds, mesh prepare, DOF distribute).
 *   - Executing the per-time-step solve via solve_equation_systems().
 */
class APEC_M_HIDDEN ApecMEquationSystems
    : public libMesh::EquationSystems,
      public FEParentObjectBase<ApecMSolver>

{
public:
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

    /** @brief Zeroes all solution vectors in all systems. */
    void zero_all_systems();

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : libMesh::EquationSystems(std::forward<Args>(args)...)
    {
    }
};

#endif // APEC_M_EQUATION_SYSTEMS_H

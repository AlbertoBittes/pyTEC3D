/**
 * @file apec_m_system.h
 * @brief APEC-M solution system types: DOF, nodal, and element solution.
 *
 * Provides three libMesh-based system classes:
 *   - ApecMLinearSystem  — transient linear-implicit system for the primary DOF unknowns
 *                         (A_x, A_y, A_z / Az, VOLT, EMF).  Owns the assembler, constraint
 *                         handler, solver configuration, and time-integration state.
 *   - ApecMExplicitSystem — lightweight explicit system used for nodal and element post-
 *                          processing results (B, H, Be, He, Jt).
 *   - ApecMSystem<SysType> — CRTP mixin that adds variable initialisation and
 *                            system-type metadata to any libMesh System subclass.
 */
#ifndef APEC_M_SYSTEM_H
#define APEC_M_SYSTEM_H

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits> // Add this include for std::is_base_of
#include <vector>

#include <libmesh/equation_systems.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/system_subset.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>
#include <libmesh/petsc_matrix.h>
#include <pybind11/pybind11.h>

#include "assembler.h"
#include "global_constraints_object.h"
#include "solver_configuration.h"
#include "apec_m_logging.h"
#include "apec_m_solver.h"
#include "apec_m_visibility.h"
#include "fe_parent_object.h"

class ElementTypesFactory;

/**
 * @brief Base interface shared by all APEC-M system types.
 *
 * Exposes the compile-time system-type enum, the name-to-type map, and
 * pure-virtual housekeeping methods implemented by each concrete system.
 */
class ApecMSystemBase
{
public:
    enum SystemType
    {
        DOF_SOLUTION_SYSTEM = 0,
        NODAL_SOLUTION,
        ELEMENT_SOLUTION,
        NUM_SYSTEMS
    };

    virtual ~ApecMSystemBase() = default;
    virtual void zero_all() = 0;
    virtual bool is_nodal() const = 0;

    static const std::map<SystemType, std::string> SYSTEM_INFO_MAP;
    static const std::string &system_name_from_type(SystemType system_type);
    static SystemType system_type_from_name(const std::string &system_name);
};

#include "variables_map.h"
#include "apec_m_equation_systems.h"
#include "element_types.h"

/**
 * @brief CRTP mixin that adds APEC-M variable initialisation to any libMesh System subclass.
 *
 * Reads the element-type map via the solver parent and adds the appropriate
 * physical DOF or result variables to the system for each subdomain.
 *
 * @tparam SysLike  A libMesh System subclass (e.g. TransientSystem<LinearImplicitSystem>).
 */
template <typename SysLike>
class ApecMSystem
    : public SysLike,
      public ApecMSystemBase
{
public:
    static_assert(std::is_base_of<libMesh::System, SysLike>::value,
                  "SysLike must be System or its subclass");

    ApecMSystemBase::SystemType system_type = ApecMSystemBase::NUM_SYSTEMS;
    std::string system_name;

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...)
    {
        std::string sys_name = this->name();
        const auto sys_type = ApecMSystemBase::system_type_from_name(sys_name);
        system_type = sys_type;
        system_name = ApecMSystemBase::system_name_from_type(sys_type);
    }

    void init_variables();

    inline bool is_nodal() const override
    {
        return system_type != ApecMSystemBase::ELEMENT_SOLUTION;
    }

    virtual void post_init() {};

    void zero_all() override
    {
        for (auto it = this->vectors_begin(); it != this->vectors_end(); ++it)
            if (it->second)
                it->second->zero();

        this->solution->zero();
        this->update();
    }

protected:
    ApecMEquationSystems &typed_equation_systems()
    {
        return static_cast<ApecMEquationSystems &>(this->get_equation_systems());
    }

    const ApecMEquationSystems &typed_equation_systems() const
    {
        return static_cast<const ApecMEquationSystems &>(this->get_equation_systems());
    }
};

/**
 * @brief Primary DOF solution system: transient linear-implicit system for A / VOLT / EMF unknowns.
 *
 * Owns the time-step state, assembler, constraint handler, and PETSc solver configuration.
 * Implements the generalised-theta time-integration scheme and optional subset-solve for
 * constrained DOF removal.
 */
class APEC_M_HIDDEN ApecMLinearSystem : public ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>
{
public:
    /** @brief Time-integration state for one time step. */
    struct TimeStep
    {
        libMesh::Real time = 0.0;       ///< Current simulation time.
        libMesh::Real theta = 1.0;      ///< Generalised-theta parameter (1=backward Euler, 0.5=Crank-Nicolson).
        bool transient_effects = false; ///< When false, the capacitance matrix contribution is skipped.
        unsigned int time_step_idx = 0; ///< Monotonically increasing time-step index.
        libMesh::Real dt = 0.0;         ///< Time-step size (current_time - previous_time).
        bool updated = true;            ///< True when dt or theta changed since the last step.
    };

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>(
              std::forward<Args>(args)...)
    {
        global_constraints_object = std::make_shared<GlobalConstraintsObject>(this);
        this->attach_constraint_object(*global_constraints_object);
        this->add_vector("solution_time_derivative", true, libMesh::GHOSTED);
        assembler = std::make_shared<AssemblerLinear>(&this->typed_equation_systems());
        this->attach_assemble_object(*assembler);
    }

    /** @brief Releases the active subset-solve restriction from the DofMap. */
    void release_subset_solve();

    /**
     * @brief Advances or updates the time-step state from the Python time-step dict.
     * @param time_step_dict  Dict with keys @c time, @c theta, @c transient_effects.
     */
    void prepare_time_step_solve(const pybind11::dict &time_step_dict);

    /** @brief Computes the solution time-derivative vector using the theta scheme. */
    void compute_solution_time_derivative();

    /** @brief Restricts the solve to unconstrained DOFs by removing constrained rows and columns. */
    void delete_rows_cols_constrained_dofs();

    /** @brief Zeroes the system matrix and all solution vectors. */
    void zero_all() override;

    std::shared_ptr<AssemblerLinear> get_assembler() const
    {
        return assembler;
    }

    std::shared_ptr<GlobalConstraintsObject> get_global_constraints_object() const
    {
        return global_constraints_object;
    }

    ApecMEquationSystems &get_apec_equation_systems()
    {
        return this->typed_equation_systems();
    }

    const ApecMEquationSystems &get_apec_equation_systems() const
    {
        return this->typed_equation_systems();
    }

    void add_matrices() override
    {
        libMesh::ImplicitSystem::add_matrices();
        this->get_system_matrix().use_hash_table(true);
    }

    libMesh::PetscLinearSolver<libMesh::Number> *get_petsc_solver() const
    {
        if (!this->linear_solver)
            apec::logging::error_and_abort("ApecMLinearSystem: linear solver is not initialized.");
        return static_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(this->linear_solver.get());
    }

    std::unique_ptr<TimeStep> time_step;     ///< State for the current time step.
    std::unique_ptr<TimeStep> old_time_step; ///< State for the previous time step (used for dt computation).

    /** @brief Called by ApecMEquationSystems after DOF distribution to finish solver setup. */
    void post_init() override
    {
        apec::logging::info("Initializing system matrix");
        this->zero_out_matrix_and_rhs = false;
        this->get_system_matrix().zero();
        apec::logging::info("Setting linear solver configurations");
        auto &equation_systems = this->typed_equation_systems();
        auto *solver_parent = equation_systems.get_parent_object();
        auto *petsc_solver = get_petsc_solver();
        if (solver_parent && petsc_solver)
        {
            solver_config_copy_ = std::make_unique<PetscSolverConfiguration>(
                solver_parent->get_apec_object().attr("solver_settings"),
                petsc_solver);
            solver_config_copy_->configure_solver();
        }
        else
        {
            apec::logging::error_and_abort("ApecMLinearSystem: solver parent or PETSc linear solver not available.");
        }
    }

private:
    std::shared_ptr<GlobalConstraintsObject> global_constraints_object;
    std::shared_ptr<AssemblerLinear> assembler;
    std::unique_ptr<PetscSolverConfiguration> solver_config_copy_;
    std::shared_ptr<libMesh::SystemSubset> subset_solve_holder_;
};

/**
 * @brief Explicit system for post-processing result fields (nodal B/H or element Be/He/Jt).
 *
 * Does not own an assembler or constraint object.  Values are written directly
 * by AssemblerLinear::write_nodal_solution() and write_element_solution().
 */
class ApecMExplicitSystem : public ApecMSystem<libMesh::ExplicitSystem>
{
public:
    // Variadic template constructor: forwards all args to ExplicitSystem
    template <typename... Args>
    ApecMExplicitSystem(Args &&...args)
        : ApecMSystem<libMesh::ExplicitSystem>(
              std::forward<Args>(args)...)
    {
        this->remove_vector("RHS");
    }
};

#endif // APEC_M_SYSTEM_H

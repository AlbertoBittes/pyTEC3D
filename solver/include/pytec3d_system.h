/**
 * @file pytec3d_system.h
 * @brief pyTEC3D solution system hierarchy.
 *
 * pyTEC3DSystem<SysLike> adds the pyTEC3DSystemBase interface to a libMesh system. The DOF system
 * (pyTEC3DDOFSystem, linear or nonlinear) solves the time step; the nodal and element solution systems are explicit
 * systems holding the fields computed from it. The templates are instantiated explicitly in pytec3d_system.cpp
 * for the libMesh systems used here.
 */
#pragma once

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <libmesh/dense_vector.h>
#include <libmesh/explicit_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/linear_solver.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/transient_system.h>
#include <libmesh/wrapped_petsc.h>
#include <petscksp.h>

#include "constraints.h"
#include "element_evaluator.h"
#include "pytec3d_system_base.h"
#include "reduced_linear_solver.h"
#include "solver_configuration.h"
#include "time_step.h"

class pyTEC3DDOFSystemBase
{
public:
    virtual ~pyTEC3DDOFSystemBase() noexcept = default;

    virtual void init_time_step(const TimeStepInput &time_step) = 0;
    virtual void prepare_time_step_solve(const TimeStepInput &time_step) = 0;
    /** @brief Computes the time derivatives of the solution at the end of the step from the time scheme. */
    virtual void update_time_derivatives() = 0;
    virtual bool check_update_system_matrix() = 0;
    virtual bool check_update_element_matrices() = 0;
    /** @brief Computes and stores the element matrices of the local elements, once per model. */
    virtual void write_element_matrices() = 0;
    virtual libMesh::LinearSolver<libMesh::Number> *get_linear_solver() const = 0;
    /** @brief The solver settings of the next solve; a change of the linear solver settings updates the system matrix. */
    virtual void apply_solver_settings(const SolverSettings &settings) = 0;
    virtual const NodalLoadInjectorElemByVarMap *ensure_nodal_load_injector_map() = 0;

    /** @brief The operators of the static nodal / element fields of each active local element (in the order of the
     * mesh's local elements), kept in memory from the assembly of the operators: {F u, F du/dt, F d2u/dt2}. */
    using FieldOperators = std::vector<std::vector<libMesh::DenseMatrix<libMesh::Number>>>;
    virtual const FieldOperators &nodal_field_operators() const = 0;
    virtual const FieldOperators &element_field_operators() const = 0;
};

template <typename SysLike>
class pyTEC3DSystem
    : public SysLike,
      public pyTEC3DSystemBase
{
public:
    pyTEC3DSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number);

    void init_variables() override;
    void solve() override;
    bool is_nodal() const override;
    void post_init() override;
    void assemble_system() override = 0;

    pyTEC3DEquationSystems &get_pytec3d_equation_systems() override;
    const pyTEC3DEquationSystems &get_pytec3d_equation_systems() const override;
    const TimeStep *current_time_step() const override;
    const TimeStep *previous_time_step() const override;

protected:
    pyTEC3DEquationSystems &typed_equation_systems();
    const pyTEC3DEquationSystems &typed_equation_systems() const;

private:
    static void assemble_dispatch(libMesh::EquationSystems &es, const std::string &name);

    pyTEC3DSystemBase::SystemType system_type_ = pyTEC3DSystemBase::NUM_SYSTEMS;
};

template <typename ImplicitSys>
class pyTEC3DDOFSystem
    : public pyTEC3DSystem<libMesh::TransientSystem<ImplicitSys>>,
      public pyTEC3DDOFSystemBase
{
public:
    pyTEC3DDOFSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number);
    ~pyTEC3DDOFSystem() noexcept override;

    void prepare_time_step_solve(const TimeStepInput &time_step) override;
    void init_time_step(const TimeStepInput &time_step) override;
    void update_time_derivatives() override;

    void solve() override;
    void post_init() override;

    bool check_update_system_matrix() override;
    bool check_update_element_matrices() override;
    void write_element_matrices() override;

    GlobalConstraintsObject<ImplicitSys> *get_global_constraints_object() const;
    libMesh::LinearSolver<libMesh::Number> *get_linear_solver() const override;
    void apply_solver_settings(const SolverSettings &settings) override;

    /** @brief The constraints of the DofMap to the linear solver: P and g (@p rows_changed: other rows, a new P). */
    void apply_constraints(bool rows_changed);
    /** @brief The solver of the linear systems: of the system, or of the Newton iterations. */
    virtual ReducedPetscLinearSolver &reduced_linear_solver() = 0;

    const NodalLoadInjectorElemByVarMap *ensure_nodal_load_injector_map() override;

protected:
    /** @brief The matrix holding K_eff = K + c_v C + c_a M: the system matrix (linear) or linear_stiffness (Newton). */
    virtual libMesh::SparseMatrix<libMesh::Number> &effective_matrix() = 0;
    /** @brief K_eff of the scheme of the current step, from the K_eff of the previous scheme (+ dc_v C + dc_a M). */
    void update_effective_matrix();
    /** @brief The right-hand side F - C h_v - M h_a of the step: the element right-hand sides (loads, magnetization,
     * nodal loads) and the history of the step as global products. */
    void assemble_rhs(libMesh::NumericVector<libMesh::Number> &F);
    /** @brief The active local elements, in order (the element loops run over them on the threads). */
    const std::vector<const libMesh::Elem *> &local_elements();

    std::unique_ptr<GlobalConstraintsObject<ImplicitSys>> global_constraints_object_;
    std::unique_ptr<PetscSolverConfiguration> solver_config_;
    std::unique_ptr<NodalLoadInjectorElemByVarMap> nodal_load_injector_map_;
    std::vector<const libMesh::Elem *> local_elements_;
    FieldOperators nodal_operators_, element_operators_;
    /// The rhs bases of the local elements (ElementEvaluator::rhs_basis()), computed with the operators; an element
    /// without one (rhs_cached_ false) has its rhs built with its FEs at every step.
    std::vector<pytec3d::ElementEvaluator::RhsBasis> rhs_bases_;
    std::vector<char> rhs_cached_;

public:
    const FieldOperators &nodal_field_operators() const override { return nodal_operators_; }
    const FieldOperators &element_field_operators() const override { return element_operators_; }

protected:
    bool element_matrices_written_ = false; ///< K, C, M assembled, the field operators stored (once per model).
    bool solver_settings_changed_ = false;  ///< The linear solver settings changed: update the system matrix.
    bool has_damping_ = false, has_mass_ = false; ///< C (M) has entries.
    libMesh::Real keff_c_v_ = 0.0, keff_c_a_ = 0.0; ///< The coefficients of C and M in effective_matrix().
};

class pyTEC3DLinearSystem : public pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>
{
public:
    pyTEC3DLinearSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number);
    ~pyTEC3DLinearSystem() noexcept override;

    void add_matrices() override;
    void assemble_system() override;
    void solve() override;
    void post_init() override;
    ReducedPetscLinearSolver &reduced_linear_solver() override;

protected:
    libMesh::SparseMatrix<libMesh::Number> &effective_matrix() override { return this->get_system_matrix(); }
};

class pyTEC3DNonLinearSystem : public pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>
{
public:
    pyTEC3DNonLinearSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number);
    ~pyTEC3DNonLinearSystem() noexcept override;

    void add_matrices() override;
    void assemble_system() override;
    void assemble_jacobian_only();
    void assemble_residual_only();

    bool check_update_system_matrix() override;
    void solve() override;
    void post_init() override;
    void apply_solver_settings(const SolverSettings &settings) override;
    ReducedPetscLinearSolver &reduced_linear_solver() override;

    /** @brief As pyTEC3DDOFSystem, and the law-output rates of the new step become the old ones. */
    void init_time_step(const TimeStepInput &time_step) override;
    /** @brief As pyTEC3DDOFSystem, and the rates of the law outputs under dt() (law_rates_). */
    void update_time_derivatives() override;

protected:
    libMesh::SparseMatrix<libMesh::Number> &effective_matrix() override { return this->get_matrix("linear_stiffness"); }

private:
    /** @brief The first and second time derivatives of the law outputs under dt() of one element (q' and q'' of
     * ElementBuilder::law_output), for the terms with 1 and 2 dt() (index 0 and 1): the history of the time scheme. */
    struct LawRates
    {
        std::array<libMesh::DenseVector<libMesh::Number>, 2> v; ///< q'
        std::array<libMesh::DenseVector<libMesh::Number>, 2> a; ///< q''
    };
    std::unordered_map<libMesh::dof_id_type, LawRates> law_rates_;     ///< At the end of the step, by element id.
    std::unordered_map<libMesh::dof_id_type, LawRates> old_law_rates_; ///< At the start of the step.

    /** @brief Newton callbacks of libMesh: assemble the Jacobian and / or the residual at the current solution. */
    static void matvec_callback(const libMesh::NumericVector<libMesh::Number> &X,
                                libMesh::NumericVector<libMesh::Number> *R,
                                libMesh::SparseMatrix<libMesh::Number> *J,
                                libMesh::NonlinearImplicitSystem &S);
    static void residual_callback(const libMesh::NumericVector<libMesh::Number> &X,
                                  libMesh::NumericVector<libMesh::Number> &R,
                                  libMesh::NonlinearImplicitSystem &S);

    void assemble_system_impl(bool assemble_jacobian, bool assemble_residual);
    /** @brief The local elements with laws selected: the only ones evaluated at each Newton iteration. */
    const std::vector<const libMesh::Elem *> &law_elements();
    std::vector<const libMesh::Elem *> law_elements_;
    bool law_elements_ready_ = false;
    /** @brief Adds the current iterate, of residual norm @p residual, to the Newton history of the step report. */
    void record_newton_iterate(libMesh::Real residual);

    bool residual_ready_from_matvec_ = false;
    std::unique_ptr<libMesh::NumericVector<libMesh::Number>> previous_iterate_; ///< Of the Newton history.
    libMesh::Real full_newton_step_ = -1.0; ///< Norm of the last full Newton step, set by the solver's monitor.
};

class pyTEC3DExplicitSystem : public pyTEC3DSystem<libMesh::ExplicitSystem>
{
public:
    pyTEC3DExplicitSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number);
};

class pyTEC3DNodalSolutionSystem : public pyTEC3DExplicitSystem
{
public:
    using pyTEC3DExplicitSystem::pyTEC3DExplicitSystem;

    void assemble_system() override;
};

class pyTEC3DElementSolutionSystem : public pyTEC3DExplicitSystem
{
public:
    using pyTEC3DExplicitSystem::pyTEC3DExplicitSystem;

    void assemble_system() override;
};

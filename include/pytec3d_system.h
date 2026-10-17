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

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <libmesh/explicit_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/linear_solver.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/system_subset.h>
#include <libmesh/transient_system.h>
#include <libmesh/wrapped_petsc.h>
#include <petscksp.h>

#include "constraints.h"
#include "pytec3d_system_base.h"
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
    virtual const NodalLoadInjectorElemByVarMap *ensure_nodal_load_injector_map() = 0;
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

    virtual void delete_rows_cols_constrained_dofs() = 0;

    const NodalLoadInjectorElemByVarMap *ensure_nodal_load_injector_map() override;

protected:
    std::unique_ptr<GlobalConstraintsObject<ImplicitSys>> global_constraints_object_;
    std::unique_ptr<libMesh::SolverConfiguration> solver_config_;
    std::unique_ptr<libMesh::SystemSubset> subset_solve_holder_;
    std::unique_ptr<NodalLoadInjectorElemByVarMap> nodal_load_injector_map_;
    bool element_matrices_written_ = false; ///< The element matrices of this model are in the emat files.
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
    void delete_rows_cols_constrained_dofs() override;

    void release_subset_solve();

private:
    void clear_reduced_solve_cache();

    std::vector<unsigned int> reduced_subset_dofs_cache_;
    libMesh::WrappedPetsc<IS> reduced_subset_is_;
    libMesh::WrappedPetsc<VecScatter> reduced_scatter_;
    libMesh::WrappedPetsc<Mat> reduced_submat_;
    libMesh::WrappedPetsc<Mat> reduced_subprecond_;
    libMesh::WrappedPetsc<Vec> reduced_subrhs_;
    libMesh::WrappedPetsc<Vec> reduced_subsolution_;

    Mat reduced_parent_mat_cache_ = nullptr;
    Mat reduced_parent_precond_mat_cache_ = nullptr;
    Vec reduced_parent_rhs_vec_cache_ = nullptr;
    bool reduced_cache_ready_ = false;
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
    void delete_rows_cols_constrained_dofs() override;

private:
    /** @brief Newton callbacks of libMesh: assemble the Jacobian and / or the residual at the current solution. */
    static void matvec_callback(const libMesh::NumericVector<libMesh::Number> &X,
                                libMesh::NumericVector<libMesh::Number> *R,
                                libMesh::SparseMatrix<libMesh::Number> *J,
                                libMesh::NonlinearImplicitSystem &S);
    static void residual_callback(const libMesh::NumericVector<libMesh::Number> &X,
                                  libMesh::NumericVector<libMesh::Number> &R,
                                  libMesh::NonlinearImplicitSystem &S);

    void assemble_system_impl(bool assemble_jacobian, bool assemble_residual);
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

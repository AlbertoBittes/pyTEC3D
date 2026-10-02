/**
 * @file pytec3d_system.h
 * @brief pyTEC3D solution system hierarchy.
 */
#ifndef pyTEC3D_SYSTEM_H
#define pyTEC3D_SYSTEM_H

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <libmesh/equation_systems.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/nonlinear_implicit_system.h>
#include <libmesh/implicit_system.h>
#include <libmesh/explicit_system.h>
#include <libmesh/system_subset.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/newton_solver.h>
#include <libmesh/nonlinear_solver.h>
#include <libmesh/wrapped_petsc.h>
#include <libmesh/linear_solver.h>
#include <petscksp.h>
#include <pybind11/pybind11.h>

#include "constraints.h"
#include "solver_configuration.h"
#include "logging.h"
#include "pytec3d_solver.h"
#include "time_step.h"

class ElementTypesFactory;
class pyTEC3DEquationSystems;

using NodalLoadInjectorElemByVarMap =
    std::unordered_map<std::string,
                       std::unordered_map<libMesh::dof_id_type, libMesh::dof_id_type>>;

class pyTEC3DSystemBase
{
public:
    enum SystemType
    {
        DOF_SOLUTION_SYSTEM = 0,
        NODAL_SOLUTION,
        ELEMENT_SOLUTION,
        NUM_SYSTEMS
    };

    virtual ~pyTEC3DSystemBase() = default;
    virtual bool is_nodal() const = 0;
    virtual void assemble_system() = 0;
    virtual pyTEC3DEquationSystems &get_pytec3d_equation_systems() = 0;
    virtual const pyTEC3DEquationSystems &get_pytec3d_equation_systems() const = 0;
    virtual const TimeStep *current_time_step() const = 0;
    virtual const TimeStep *previous_time_step() const = 0;
    virtual void solve() = 0;
    virtual void init_variables() = 0;
    virtual void post_init() = 0;

    static const std::map<SystemType, std::string> SYSTEM_INFO_MAP;
    static const std::string &system_name_from_type(SystemType system_type);
    static SystemType system_type_from_name(const std::string &system_name);
};

class pyTEC3DDOFSystemBase
{
public:
    virtual ~pyTEC3DDOFSystemBase() noexcept = default;

    virtual void init_time_step(const pybind11::dict &time_step_dict) = 0;
    virtual void prepare_time_step_solve(const pybind11::dict &time_step_dict) = 0;
    virtual void compute_solution_time_derivative() = 0;
    virtual bool check_update_system_matrix() = 0;
    virtual bool check_update_element_matrices() = 0;
    virtual void write_element_matrices() = 0;
    virtual libMesh::LinearSolver<libMesh::Number> *get_linear_solver() const = 0;
    virtual const NodalLoadInjectorElemByVarMap *ensure_nodal_load_injector_map() = 0;
};

#include "variables_map.h"
#include "pytec3d_equation_systems.h"
#include "element_types.h"

template <typename SysLike>
class pyTEC3DSystem
    : public SysLike,
      public pyTEC3DSystemBase
{
public:
    static_assert(std::is_base_of<libMesh::System, SysLike>::value,
                  "SysLike must be System or its subclass");

    pyTEC3DSystemBase::SystemType system_type = pyTEC3DSystemBase::NUM_SYSTEMS;
    std::string system_name;

    template <typename... Args>
    pyTEC3DSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...)
    {
        std::string sys_name = this->name();
        const auto sys_type = pyTEC3DSystemBase::system_type_from_name(sys_name);
        system_type = sys_type;
        system_name = pyTEC3DSystemBase::system_name_from_type(sys_type);
        this->attach_assemble_function(&pyTEC3DSystem<SysLike>::assemble_dispatch);
    }

    void init_variables() override;

    void solve() override
    {
        SysLike::solve();
    }

    inline bool is_nodal() const override
    {
        return system_type != pyTEC3DSystemBase::ELEMENT_SOLUTION;
    }

    virtual void post_init() override {}

    virtual void assemble_system() override = 0;

    pyTEC3DEquationSystems &get_pytec3d_equation_systems() override
    {
        return this->typed_equation_systems();
    }

    const pyTEC3DEquationSystems &get_pytec3d_equation_systems() const override
    {
        return this->typed_equation_systems();
    }

    const TimeStep *current_time_step() const override
    {
        return this->typed_equation_systems().current_time_step();
    }

    const TimeStep *previous_time_step() const override
    {
        return this->typed_equation_systems().previous_time_step();
    }

protected:
    pyTEC3DEquationSystems &typed_equation_systems()
    {
        return static_cast<pyTEC3DEquationSystems &>(this->get_equation_systems());
    }

    const pyTEC3DEquationSystems &typed_equation_systems() const
    {
        return static_cast<const pyTEC3DEquationSystems &>(this->get_equation_systems());
    }

private:
    static void assemble_dispatch(libMesh::EquationSystems &es, const std::string &name)
    {
        auto &sys = es.get_system(name);
        auto *pytec3d_sys = dynamic_cast<pyTEC3DSystemBase *>(&sys);
        if (!pytec3d_sys)
            pytec3d::logging::error_and_abort("pyTEC3DSystem::assemble_dispatch: system is not an pyTEC3DSystemBase.");
        pytec3d_sys->assemble_system();
    }
};

template <typename ImplicitSys>
class pyTEC3DDOFSystem
    : public pyTEC3DSystem<libMesh::TransientSystem<ImplicitSys>>,
      public pyTEC3DDOFSystemBase
{
public:
    ~pyTEC3DDOFSystem() noexcept override = default;

    template <typename... Args>
    pyTEC3DDOFSystem(Args &&...args)
        : pyTEC3DSystem<libMesh::TransientSystem<ImplicitSys>>(
              std::forward<Args>(args)...)
    {
        global_constraints_object = std::make_shared<GlobalConstraintsObject<ImplicitSys>>(this);
        this->attach_constraint_object(*global_constraints_object);
        this->add_vector("solution_time_derivative", true, libMesh::GHOSTED);
    }

    void prepare_time_step_solve(const pybind11::dict &time_step_dict) override;
    void init_time_step(const pybind11::dict &time_step_dict) override;
    void compute_solution_time_derivative() override;

    void solve() override
    {
        libMesh::TransientSystem<ImplicitSys>::solve();
    }

    void init_variables() override
    {
        pyTEC3DSystem<libMesh::TransientSystem<ImplicitSys>>::init_variables();
    }

    void post_init() override;

    virtual bool check_update_system_matrix() override;
    bool check_update_element_matrices() override;
    void write_element_matrices() override;

    std::shared_ptr<GlobalConstraintsObject<ImplicitSys>> get_global_constraints_object() const
    {
        return global_constraints_object;
    }

    libMesh::LinearSolver<libMesh::Number> *get_linear_solver() const override
    {
        if (!this->linear_solver)
            return nullptr;
        return this->linear_solver.get();
    }

    virtual void delete_rows_cols_constrained_dofs() = 0;

    const NodalLoadInjectorElemByVarMap *ensure_nodal_load_injector_map() override;

protected:
    std::shared_ptr<GlobalConstraintsObject<ImplicitSys>> global_constraints_object;
    std::unique_ptr<libMesh::SolverConfiguration> solver_config_copy_;
    std::shared_ptr<libMesh::SystemSubset> subset_solve_holder_;
    std::shared_ptr<NodalLoadInjectorElemByVarMap> nodal_load_injector_map_;
};

class pyTEC3DLinearSystem : public pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>
{
public:
    ~pyTEC3DLinearSystem() noexcept override = default;

    template <typename... Args>
    pyTEC3DLinearSystem(Args &&...args)
        : pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>(std::forward<Args>(args)...)
    {
    }

    void add_matrices() override
    {
        libMesh::ImplicitSystem::add_matrices();
        this->get_system_matrix().use_hash_table(true);
    }

    void assemble_system() override;

    void solve() override;

    void post_init() override;

    void release_subset_solve();

    void delete_rows_cols_constrained_dofs() override;

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
    ~pyTEC3DNonLinearSystem() noexcept override = default;

    template <typename... Args>
    pyTEC3DNonLinearSystem(Args &&...args)
        : pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>(std::forward<Args>(args)...)
    {
        diff_solver = std::make_unique<libMesh::NewtonSolver>(*this);
        diff_solver->verbose = true;
        auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(diff_solver.get());
        if (!newton_solver)
            pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem ctor: failed to create Newton solver.");
        newton_solver->require_residual_reduction = true;
        auto nonlinear_solver = this->nonlinear_solver.get();
        nonlinear_solver->matvec = &pyTEC3DNonLinearSystem::matvec_callback;
        nonlinear_solver->residual = &pyTEC3DNonLinearSystem::residual_callback;
    }

    void add_matrices() override
    {
        libMesh::ImplicitSystem::add_matrices();
        this->get_system_matrix().use_hash_table(true);
        this->add_matrix("linear_stiffness");
        this->get_matrix("linear_stiffness").use_hash_table(true);
    }

    static void matvec_callback(const libMesh::NumericVector<libMesh::Number> &X,
                                libMesh::NumericVector<libMesh::Number> *R,
                                libMesh::SparseMatrix<libMesh::Number> *J,
                                libMesh::NonlinearImplicitSystem &S)
    {
        (void)X;
        auto *sys = dynamic_cast<pyTEC3DNonLinearSystem *>(&S);
        if (!sys)
            pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::matvec_callback: system is not pyTEC3DNonLinearSystem.");
        if (R && J)
        {
            sys->assemble_system();
            sys->residual_ready_from_matvec_ = true;
        }
        else if (J)
        {
            sys->assemble_jacobian_only();
            sys->residual_ready_from_matvec_ = false;
        }
        else if (R)
        {
            sys->assemble_residual_only();
            sys->residual_ready_from_matvec_ = false;
        }
    }

    static void residual_callback(const libMesh::NumericVector<libMesh::Number> &X,
                                  libMesh::NumericVector<libMesh::Number> &R,
                                  libMesh::NonlinearImplicitSystem &S)
    {
        (void)X;
        (void)R;
        auto *sys = dynamic_cast<pyTEC3DNonLinearSystem *>(&S);
        if (!sys)
            pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::residual_callback: system is not pyTEC3DNonLinearSystem.");
        if (sys->residual_ready_from_matvec_)
        {
            sys->residual_ready_from_matvec_ = false;
            return;
        }
        sys->assemble_residual_only();
    }

    void assemble_system() override;
    void assemble_jacobian_only();
    void assemble_residual_only();

    bool check_update_system_matrix() override;

    void solve() override;

    void post_init() override;

    void delete_rows_cols_constrained_dofs() override;

protected:
    void assemble_system_impl(bool assemble_jacobian, bool assemble_residual);
    bool residual_ready_from_matvec_ = false;
    std::unique_ptr<libMesh::SolverConfiguration> nonlinear_solver_config_copy_;
};

class pyTEC3DExplicitSystem : public pyTEC3DSystem<libMesh::ExplicitSystem>
{
public:
    template <typename... Args>
    pyTEC3DExplicitSystem(Args &&...args)
        : pyTEC3DSystem<libMesh::ExplicitSystem>(
              std::forward<Args>(args)...)
    {
        this->remove_vector("RHS");
    }
};

class pyTEC3DNodalSolutionSystem : public pyTEC3DExplicitSystem
{
public:
    template <typename... Args>
    pyTEC3DNodalSolutionSystem(Args &&...args)
        : pyTEC3DExplicitSystem(std::forward<Args>(args)...)
    {
    }

    void assemble_system() override;
};

class pyTEC3DElementSolutionSystem : public pyTEC3DExplicitSystem
{
public:
    template <typename... Args>
    pyTEC3DElementSolutionSystem(Args &&...args)
        : pyTEC3DExplicitSystem(std::forward<Args>(args)...)
    {
    }

    void assemble_system() override;
};

#endif // pyTEC3D_SYSTEM_H

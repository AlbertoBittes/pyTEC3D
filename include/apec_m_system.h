/**
 * @file apec_m_system.h
 * @brief APEC-M solution system hierarchy.
 */
#ifndef APEC_M_SYSTEM_H
#define APEC_M_SYSTEM_H

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
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
#include <pybind11/pybind11.h>

#include "global_constraints_object.h"
#include "solver_configuration.h"
#include "apec_m_logging.h"
#include "apec_m_solver.h"
#include "apec_m_visibility.h"

class ElementTypesFactory;
class ApecMEquationSystems;

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
    virtual bool is_nodal() const = 0;
    virtual void assemble_system() = 0;

    static const std::map<SystemType, std::string> SYSTEM_INFO_MAP;
    static const std::string &system_name_from_type(SystemType system_type);
    static SystemType system_type_from_name(const std::string &system_name);
};

class ApecMDOFSystemBase
{
public:
    struct TimeStep
    {
        libMesh::Real time = 0.0;
        libMesh::Real theta = 1.0;
        bool transient_effects = false;
        unsigned int time_step_idx = 0;
        libMesh::Real dt = 0.0;
        bool updated = true;
    };

    virtual ~ApecMDOFSystemBase() = default;

    virtual void init_time_step(const pybind11::dict &time_step_dict) = 0;
    virtual void prepare_time_step_solve(const pybind11::dict &time_step_dict) = 0;
    virtual void compute_solution_time_derivative() = 0;
    virtual void solve() = 0;
    virtual void init_variables() = 0;
    virtual void post_init() = 0;
    virtual void write_element_data() = 0;
    virtual void write_element_matrices() = 0;
    virtual void write_element_tangent_stiffness_residual() = 0;
    virtual void write_element_rhs() = 0;
    virtual bool check_update_system_matrix() = 0;
    virtual bool check_update_element_matrices() = 0;
    virtual std::shared_ptr<GlobalConstraintsObject> get_global_constraints_object() const = 0;
    virtual libMesh::PetscLinearSolver<libMesh::Number> *get_petsc_solver() const = 0;
    virtual ApecMEquationSystems &get_apec_equation_systems() = 0;
    virtual const ApecMEquationSystems &get_apec_equation_systems() const = 0;
    virtual const TimeStep *current_time_step() const = 0;
    virtual const TimeStep *previous_time_step() const = 0;
};

#include "variables_map.h"
#include "apec_m_equation_systems.h"
#include "element_types.h"

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

    template <typename... Args>
    ApecMSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...)
    {
        std::string sys_name = this->name();
        const auto sys_type = ApecMSystemBase::system_type_from_name(sys_name);
        system_type = sys_type;
        system_name = ApecMSystemBase::system_name_from_type(sys_type);
        this->attach_assemble_function(&ApecMSystem<SysLike>::assemble_dispatch);
    }

    void init_variables();

    inline bool is_nodal() const override
    {
        return system_type != ApecMSystemBase::ELEMENT_SOLUTION;
    }

    virtual void post_init() {}

    virtual void assemble_system() override = 0;

protected:
    ApecMEquationSystems &typed_equation_systems()
    {
        return static_cast<ApecMEquationSystems &>(this->get_equation_systems());
    }

    const ApecMEquationSystems &typed_equation_systems() const
    {
        return static_cast<const ApecMEquationSystems &>(this->get_equation_systems());
    }

private:
    static void assemble_dispatch(libMesh::EquationSystems &es, const std::string &name)
    {
        auto &sys = es.get_system(name);
        auto *apec_sys = dynamic_cast<ApecMSystemBase *>(&sys);
        if (!apec_sys)
            apec::logging::error_and_abort("ApecMSystem::assemble_dispatch: system is not an ApecMSystemBase.");
        apec_sys->assemble_system();
    }
};

template <typename ImplicitSys>
class ApecMDOFSystem
    : public ApecMSystem<libMesh::TransientSystem<ImplicitSys>>,
      public ApecMDOFSystemBase
{
public:
    template <typename... Args>
    ApecMDOFSystem(Args &&...args)
        : ApecMSystem<libMesh::TransientSystem<ImplicitSys>>(
              std::forward<Args>(args)...)
    {
        global_constraints_object = std::make_shared<GlobalConstraintsObject>(this);
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
        ApecMSystem<libMesh::TransientSystem<ImplicitSys>>::init_variables();
    }

    void post_init() override
    {
        ApecMSystem<libMesh::TransientSystem<ImplicitSys>>::post_init();
    }

    void write_element_data() override;
    void write_element_matrices() override;
    void write_element_tangent_stiffness_residual() override;
    void write_element_rhs() override;

    bool check_update_system_matrix() override;
    bool check_update_element_matrices() override;

    std::shared_ptr<GlobalConstraintsObject> get_global_constraints_object() const override
    {
        return global_constraints_object;
    }

    libMesh::PetscLinearSolver<libMesh::Number> *get_petsc_solver() const override
    {
        if (!this->linear_solver)
            return nullptr;
        return dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(this->linear_solver.get());
    }

    ApecMEquationSystems &get_apec_equation_systems() override
    {
        return this->typed_equation_systems();
    }

    const ApecMEquationSystems &get_apec_equation_systems() const override
    {
        return this->typed_equation_systems();
    }

    const TimeStep *current_time_step() const override
    {
        return time_step.get();
    }

    const TimeStep *previous_time_step() const override
    {
        return old_time_step.get();
    }

    std::unique_ptr<TimeStep> time_step;
    std::unique_ptr<TimeStep> old_time_step;

protected:
    std::shared_ptr<GlobalConstraintsObject> global_constraints_object;
    std::unique_ptr<PetscSolverConfiguration> solver_config_copy_;
};

class APEC_M_HIDDEN ApecMLinearSystem : public ApecMDOFSystem<libMesh::LinearImplicitSystem>
{
public:
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : ApecMDOFSystem<libMesh::LinearImplicitSystem>(std::forward<Args>(args)...)
    {
    }

    void add_matrices() override
    {
        libMesh::ImplicitSystem::add_matrices();
        this->get_system_matrix().use_hash_table(true);
    }

    void assemble_system() override;

    void solve() override;

    void release_subset_solve();

    void delete_rows_cols_constrained_dofs();

    void post_init() override;

private:
    std::shared_ptr<libMesh::SystemSubset> subset_solve_holder_;
};

class APEC_M_HIDDEN ApecMNonLinearSystem : public ApecMDOFSystem<libMesh::NonlinearImplicitSystem>
{
public:
    template <typename... Args>
    ApecMNonLinearSystem(Args &&...args)
        : ApecMDOFSystem<libMesh::NonlinearImplicitSystem>(std::forward<Args>(args)...)
    {
        newton_solver = std::make_unique<libMesh::NewtonSolver>(*this);
    }

    void add_matrices() override
    {
        libMesh::ImplicitSystem::add_matrices();
        this->get_system_matrix().use_hash_table(true);
        this->add_matrix("linear_stiffness");
        this->get_matrix("linear_stiffness").use_hash_table(true);
    }

    void assemble_system() override;

    void solve() override;

    void post_init() override;

protected:
    std::unique_ptr<libMesh::NewtonSolver> newton_solver = nullptr;
    unsigned int max_iterations = 50;
    libMesh::Real convergence_tolerance = 1e-6;
};

class ApecMExplicitSystem : public ApecMSystem<libMesh::ExplicitSystem>
{
public:
    template <typename... Args>
    ApecMExplicitSystem(Args &&...args)
        : ApecMSystem<libMesh::ExplicitSystem>(
              std::forward<Args>(args)...)
    {
        this->remove_vector("RHS");
    }
};

class ApecMNodalSolutionSystem : public ApecMExplicitSystem
{
public:
    template <typename... Args>
    ApecMNodalSolutionSystem(Args &&...args)
        : ApecMExplicitSystem(std::forward<Args>(args)...)
    {
    }

    void assemble_system() override;
};

class ApecMElementSolutionSystem : public ApecMExplicitSystem
{
public:
    template <typename... Args>
    ApecMElementSolutionSystem(Args &&...args)
        : ApecMExplicitSystem(std::forward<Args>(args)...)
    {
    }

    void assemble_system() override;
};

#endif // APEC_M_SYSTEM_H

#ifndef APEC_M_SYSTEM_H
#define APEC_M_SYSTEM_H

#include <cmath>
#include <limits>
#include <map>
#include <memory>
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
#include "apec_m_solver.h"
#include "fe_parent_object.h"

using namespace libMesh;

class ElementTypesFactory;

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

    static const std::map<SystemType, std::string> SYSTEM_INFO_MAP;
    static const std::string &system_name_from_type(SystemType system_type);
};

#include "variables_map.h"
#include "apec_m_equation_systems.h"
#include "element_types.h"

// The basic libMesh class (System, EquationSystems, see libMesh guide) are extendend to meet APEC-M code structure requirements
template <typename SysLike>
class ApecMSystem
    : public SysLike,
      public ApecMSystemBase,
{
public:
    using ContextView = TypedContextView<ElementTypesFactory, PetscSolverConfiguration>;

    static_assert(std::is_base_of<System, SysLike>::value,
                  "SysLike must be System or its subclass");

    ApecMSystemBase::SystemType system_type = ApecMSystemBase::NUM_SYSTEMS;
    std::string system_name;

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...)
    {
        // Fetch system type from Globals using the system name
        std::string sys_name = this->name();
        auto sys_type = ApecMSystemBase::NUM_SYSTEMS;
        for (const auto &pair : ApecMSystemBase::SYSTEM_INFO_MAP)
        {
            if (pair.second == sys_name)
            {
                sys_type = pair.first;
                break;
            }
        }
        if (sys_type == ApecMSystemBase::NUM_SYSTEMS)
            libmesh_error_msg("System type not found for system name: " + sys_name);
        auto it = ApecMSystemBase::SYSTEM_INFO_MAP.find(sys_type);
        if (it == ApecMSystemBase::SYSTEM_INFO_MAP.end())
            libmesh_error_msg("System name not found for system type.");
        system_type = sys_type;
        system_name = it->second;
    }

    void init_variables();

    virtual void post_init() {};

    void zero_all() override
    {
        for (auto it = this->vectors_begin(); it != this->vectors_end(); ++it)
            if (it->second)
                it->second->zero();

        this->solution->zero();
        this->update();
    }
};

class ApecMLinearSystem : public ApecMSystem<TransientSystem<LinearImplicitSystem>>
{
public:
    struct TimeStep
    {
        pybind11::dict time_step_dict;
        unsigned int time_step_idx = 0;
        Real dt = 0.0;
        bool updated = true;
    };

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : ApecMSystem<TransientSystem<LinearImplicitSystem>>(
              std::forward<Args>(args)...)
    {
        global_constraints_object = std::make_shared<GlobalConstraintsObject>(this);
        this->attach_constraint_object(*global_constraints_object);
        this->add_vector("solution_time_derivative", true, GHOSTED);
        auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&this->get_equation_systems());
        assembler = std::make_shared<AssemblerLinear>(equation_systems);
        this->attach_assemble_object(*assembler);
    }

    void release_subset_solve();

    void prepare_time_step_solve(const pybind11::dict &time_step_dict);
    void compute_solution_time_derivative(const pybind11::dict &time_step);
    void delete_rows_cols_constrained_dofs();
    void zero_all() override;

    std::shared_ptr<AssemblerLinear> get_assembler() const
    {
        return assembler;
    }

    std::shared_ptr<GlobalConstraintsObject> get_global_constraints_object() const
    {
        return global_constraints_object;
    }

    void add_matrices() override
    {
        ImplicitSystem::add_matrices();
        this->get_system_matrix().use_hash_table(true);
    }

    std::unique_ptr<TimeStep> time_step;
    std::unique_ptr<TimeStep> old_time_step;

    void post_init() override
    {
        this->zero_out_matrix_and_rhs = false;
        this->get_system_matrix().zero();
        // These are maybe not needed
        // this->rhs->zero();
        // this->current_local_solution->zero();
        // this->current_local_solution->close();
        auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&this->get_equation_systems());
        auto *solver_parent = equation_systems ? equation_systems->get_parent_object() : nullptr;
        auto *petsc_solver = dynamic_cast<PetscLinearSolver<Number> *>(this->linear_solver.get());
        if (solver_parent && petsc_solver)
        {
            solver_config_copy_.reset(new PetscSolverConfiguration(
                solver_parent->get_apec_object().attr("solver_settings"),
                petsc_solver));
            solver_config_copy_->configure_solver();
        }
        else
        {
            libmesh_error_msg("ApecMLinearSystem: solver parent or PETSc linear solver not available.");
        }
    }

private:
    std::shared_ptr<GlobalConstraintsObject> global_constraints_object;
    std::shared_ptr<AssemblerLinear> assembler;
    std::unique_ptr<PetscSolverConfiguration> solver_config_copy_;
};

class ApecMExplicitSystem : public ApecMSystem<ExplicitSystem>
{
public:
    // Variadic template constructor: forwards all args to ExplicitSystem
    template <typename... Args>
    ApecMExplicitSystem(Args &&...args)
        : ApecMSystem<ExplicitSystem>(
              std::forward<Args>(args)...)
    {
        this->remove_vector("RHS");
    }
};

#endif // APEC_M_SYSTEM_H

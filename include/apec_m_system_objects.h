#ifndef APEC_M_SYSTEM_OBJECTS_H
#define APEC_M_SYSTEM_OBJECTS_H

#include <memory>
#include <type_traits> // Add this include for std::is_base_of

#include <libmesh/equation_systems.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/system_subset.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/periodic_boundaries.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>
#include <libmesh/petsc_matrix.h>

#include "fem_service.h"
#include "element_types.h"
#include "element_types_factory.h"
#include "assembler.h"
#include "globals.h"
#include "global_constraints_object.h"
#include "solver_configuration.h"
#include "apec_m_eq_system_object.h"

using namespace libMesh;

class ApecMSystemBase
{
public:
    virtual ~ApecMSystemBase() = default;
    virtual void zero_all() = 0;
};

// The basic libMesh class (System, EquationSystems, see libMesh guide) are extendend to meet APEC-M code structure requirements
template <typename SysLike>
class ApecMSystem
    : public SysLike,
      public ApecMSystemBase,
      public FEMService<ElementTypesFactory, PetscSolverConfiguration>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory, PetscSolverConfiguration>;

    static_assert(std::is_base_of<System, SysLike>::value,
                  "SysLike must be System or its subclass");

    Globals::SystemInfo system_info; // Attribute to hold SystemInfo

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...)
    {
        // Fetch SystemInfo from Globals using the system name
        std::string sys_name = this->name();
        auto sys_type = Globals::get_key_by_predicate_one(
            Globals::SYSTEM_INFO_MAP,
            [sys_name](const Globals::SystemInfo &info)
            {
                return info.name == sys_name;
            });
        system_info = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, sys_type);
    }

    static void ramp_vectors(const Substep &time_substep, NumericVector<Number> &old_vector, NumericVector<Number> &current_vector, NumericVector<Number> &ramped_vector);

    void init_variables();

    virtual void post_init() {};

    void zero_all() override
    {
        for (auto it = this->vectors_begin(); it != this->vectors_end(); ++it)
            if (it->second)
                it->second->zero();

        auto zero_optional = [](NumericVector<Number> *vec)
        {
            if (vec)
                vec->zero();
        };

        zero_optional(this->solution.get());
        zero_optional(this->current_local_solution.get());
        // zero_optional(this->old_local_solution.get());
        // zero_optional(this->older_local_solution.get());
        // zero_optional(this->rhs.get());
    }

private:
};

class ApecMLinearSystem : public ApecMSystem<TransientSystem<LinearImplicitSystem>>
{
public:
    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : ApecMSystem<TransientSystem<LinearImplicitSystem>>(
              std::forward<Args>(args)...)
    {
        global_constraints_object = std::make_shared<GlobalConstraintsObject>(this);
        this->attach_constraint_object(*global_constraints_object);
        this->add_vector("node_dof_constraints");
        this->add_vector("old_node_dof_constraints");
        this->add_vector("ramped_node_dof_constraints");
        this->add_vector("solution_time_derivative", true, GHOSTED);
        assembler = std::make_shared<AssemblerLinear>(static_cast<ApecMEquationSystems *>(&this->get_equation_systems()));
        this->attach_assemble_object(*assembler);
    }

    void release_subset_solve();

    void compute_solution_time_derivative(const Substep &time_substep);
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

        if (this->get_dof_map().get_periodic_boundaries()->size() > 0)
        {
            this->get_system_matrix().use_hash_table(true);
        }
    }

    void post_init() override
    {
        this->zero_out_matrix_and_rhs = false;
        this->get_system_matrix().zero();
        this->rhs->zero();
        this->current_local_solution->zero();
        this->current_local_solution->close();
        auto solver_config = this->context<PetscSolverConfiguration>();

        if (solver_config)
        {
            solver_config->set_petsc_linear_solver(
                dynamic_cast<PetscLinearSolver<Number> *>(this->linear_solver.get()));
            solver_config->configure_solver();
        }
        else
        {
            libmesh_error_msg("ApecMLinearSystem: No PetscSolverConfiguration found in system context.");
        }
    }

private:
    std::shared_ptr<GlobalConstraintsObject> global_constraints_object;
    std::shared_ptr<AssemblerLinear> assembler;
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
        if (system_info.type == Globals::NODAL_LOADS || system_info.type == Globals::ELEMENT_LOADS)
        {
            this->add_vector("old_loads");
            this->add_vector("ramped_loads");
        }
    }
};

#endif // APEC_M_SYSTEM_OBJECTS_H

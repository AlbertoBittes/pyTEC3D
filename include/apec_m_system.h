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

#include "fem_service.h"
#include "assembler.h"
#include "global_constraints_object.h"
#include "solver_configuration.h"
#include "coordinate_system.h"
#include "time_step.h"
#include "fe_parent_object.h"

using namespace libMesh;

class ElementTypesFactory;

class ApecMSystemBase
{
public:
    enum SystemType
    {
        DOF_SOLUTION_SYSTEM = 0,
        NODAL_LOADS,
        ELEMENT_LOADS,
        NODAL_SOLUTION,
        ELEMENT_SOLUTION,
        SUPPORT_SOLUTION,
        NODAL_B,
        NUM_SYSTEMS
    };

    virtual ~ApecMSystemBase() = default;
    virtual void zero_all() = 0;
    virtual bool global_to_local_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                            const std::string &variable_label,
                                                            const std::vector<dof_id_type> &selected_dof_ids = {}) = 0;
    virtual bool local_to_global_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                            const std::string &variable_label,
                                                            const std::vector<dof_id_type> &selected_dof_ids = {}) = 0;
    virtual void ramp_vectors(const Substep &time_substep,
                              NumericVector<Number> &old_vector,
                              NumericVector<Number> &current_vector,
                              NumericVector<Number> &ramped_vector) = 0;

    static const std::map<SystemType, std::string> SYSTEM_INFO_MAP;
    static bool is_nodal_system(SystemType system_type);
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
      public FEMService<ElementTypesFactory, PetscSolverConfiguration>
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

    void ramp_vectors(const Substep &time_substep,
                      NumericVector<Number> &old_vector,
                      NumericVector<Number> &current_vector,
                      NumericVector<Number> &ramped_vector) override;

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

    bool global_to_local_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                    const std::string &variable_label,
                                                    const std::vector<dof_id_type> &selected_dof_ids = {}) override;
    bool local_to_global_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                    const std::string &variable_label,
                                                    const std::vector<dof_id_type> &selected_dof_ids = {}) override;
    void dof_numbers_vars(const std::string &variable_name,
                          bool select_elements,
                          const std::vector<dof_id_type> &object_ids,
                          std::vector<dof_id_type> &dof_ids,
                          bool global);
};

template <typename SysLike>
void ApecMSystem<SysLike>::ramp_vectors(const Substep &time_substep,
                                        NumericVector<Number> &old_vector,
                                        NumericVector<Number> &current_vector,
                                        NumericVector<Number> &ramped_vector)
{
    if (time_substep.get_ramped_load())
    {
        Real substep_number = static_cast<Real>(time_substep.get_index()) + 1;
        Real number_of_substeps = static_cast<Real>(time_substep.get_parent_object()->get_n_substeps());
        for (numeric_index_type k = ramped_vector.first_local_index(); k < ramped_vector.last_local_index(); ++k)
        {
            Real old_val = old_vector(k);
            Real curr_val = current_vector(k);

            if (std::isnan(curr_val))
            {
                ramped_vector.set(k, std::numeric_limits<Real>::quiet_NaN());
            }
            else
            {
                Real start_val = std::isnan(old_val) ? 0.0 : old_val;
                Real ramped_val = start_val + (curr_val - start_val) / number_of_substeps * substep_number;
                ramped_vector.set(k, ramped_val);
            }
        }
    }
    else
    {
        for (numeric_index_type k = ramped_vector.first_local_index(); k < ramped_vector.last_local_index(); ++k)
        {
            Real curr_val = current_vector(k);
            ramped_vector.set(k, curr_val);
        }
    }
    ramped_vector.close();
    this->update();
}

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
        auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&this->get_equation_systems());
        assembler = std::make_shared<AssemblerLinear>(equation_systems);
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
        this->get_system_matrix().use_hash_table(true);
    }

    void post_init() override
    {
        this->zero_out_matrix_and_rhs = false;
        this->get_system_matrix().zero();
        // These are maybe not needed
        // this->rhs->zero();
        // this->current_local_solution->zero();
        // this->current_local_solution->close();
        auto solver_config = this->context<PetscSolverConfiguration>();

        if (solver_config)
        {
            solver_config_copy_.reset(new PetscSolverConfiguration(*solver_config));
            solver_config_copy_->set_petsc_linear_solver(
                dynamic_cast<PetscLinearSolver<Number> *>(this->linear_solver.get()));
            solver_config_copy_->configure_solver();
        }
        else
        {
            libmesh_error_msg("ApecMLinearSystem: No PetscSolverConfiguration found in system context.");
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
        if (system_type == ApecMSystemBase::NODAL_LOADS || system_type == ApecMSystemBase::ELEMENT_LOADS)
        {
            this->add_vector("old_loads");
            this->add_vector("new_loads");
        }
    }
};

#endif // APEC_M_SYSTEM_H

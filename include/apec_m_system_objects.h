#ifndef APEC_M_SYSTEM_OBJECTS_H
#define APEC_M_SYSTEM_OBJECTS_H

#include <memory>
#include <type_traits> // Add this include for std::is_base_of

#include <libmesh/equation_systems.h>
#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/system_subset.h>

#include "fem_service_abstract.h"
#include "element_types.h"
#include "element_types_factory.h"
#include "assembler.h"
#include "globals.h"
#include "global_constraints_object.h"

using namespace libMesh;

class ElementBuilder;
// class GlobalConstraintsObject;

// The basic libMesh class (System, EquationSystems, see libMesh guide) are extendend to meet APEC-M code structure requirements
template <typename SysLike, typename AssemblerLike>
class ApecMSystem
    : public SysLike,
      public FEMServiceAbstract<ElementTypesFactory>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory>;

    static_assert(std::is_base_of<System, SysLike>::value,
                  "SysLike must be System or its subclass");

    static_assert(std::is_base_of<System::Assembly, AssemblerLike>::value,
                  "AssemblerLike must be System::Assembly or its subclass");

    Globals::SystemInfo system_info; // Attribute to hold SystemInfo

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...),
          FEMServiceAbstract<ElementTypesFactory>(nullptr)
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
        this->set_assembler();
    }

    void set_assembler()
    {
        assembler = std::make_shared<AssemblerLike>(this);
        this->attach_assemble_object(*assembler);
    }

    std::shared_ptr<AssemblerLike> get_assembler() const
    {
        return assembler;
    }

private:
    std::shared_ptr<AssemblerLike> assembler;
};

class ApecMLinearSystem : public ApecMSystem<TransientSystem<LinearImplicitSystem>, AssemblerLinear>
{
public:
    using ContextView = TypedContextView<ElementTypesFactory>;

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : ApecMSystem<TransientSystem<LinearImplicitSystem>, AssemblerLinear>(
              std::forward<Args>(args)...)
    {
        global_constraints_object = std::make_shared<GlobalConstraintsObject>(this);
        this->attach_constraint_object(*global_constraints_object);
        this->add_vector("node_dof_constraints");
        this->add_vector("old_node_dof_constraints");
    }

    std::shared_ptr<GlobalConstraintsObject> get_global_constraints_object() const
    {
        return global_constraints_object;
    }

    void release_subset_solve();

    void compute_solution_time_derivative(const Substep &time_substep);

private:
    std::shared_ptr<GlobalConstraintsObject> global_constraints_object;
};

typedef ApecMSystem<ExplicitSystem, AssemblerSolution> ApecMSolutionSystem;
typedef ApecMSystem<TransientSystem<ExplicitSystem>, AssemblerSolution> ApecMLoadsSystem;

#endif // APEC_M_SYSTEM_OBJECTS_H
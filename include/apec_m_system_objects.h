#ifndef APEC_M_SYSTEM_OBJECTS_H
#define APEC_M_SYSTEM_OBJECTS_H

#include <memory>
#include <type_traits> // Add this include for std::is_base_of

#include <libmesh/transient_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/system_subset.h>
#include "fem_service_abstract.h"
#include "element_types.h"
#include "boundary_conditions_object.h"
#include "assembler.h"
#include "system_type.h"

using namespace libMesh;

class ElementBuilder;

// The basic libMesh class (System, EquationSystems, see libMesh guide) are extendend to meet APEC-M code structure requirements
template <typename SysLike, typename AssemblerLike>
class ApecMSystem
    : public SysLike,
      public FEMServiceAbstract<ElementTypeFactory>
{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    static_assert(std::is_base_of<System, SysLike>::value,
                  "SysLike must be System or its subclass");

    static_assert(std::is_base_of<System::Assembly, AssemblerLike>::value,
                  "AssemblerLike must be System::Assembly or its subclass");

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMSystem(Args &&...args)
        : SysLike(std::forward<Args>(args)...),
          FEMServiceAbstract<ElementTypeFactory>(nullptr) // or whatever default you want
    {
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
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : ApecMSystem<TransientSystem<LinearImplicitSystem>, AssemblerLinear>(
              std::forward<Args>(args)...)
    {
        boundary_conditions_object = std::make_shared<BoundaryConditionsObject>(this);
        this->attach_constraint_object(*boundary_conditions_object);
    }

    std::shared_ptr<BoundaryConditionsObject> get_boundary_conditions_object() const
    {
        return boundary_conditions_object;
    }

    void release_subset_solve();

    void compute_solution_time_derivative(const Substep &time_substep);

private:
    std::shared_ptr<BoundaryConditionsObject> boundary_conditions_object;
};

typedef ApecMSystem<ExplicitSystem, AssemblerSolution> ApecMSolutionSystem;

class ApecMEquationSystems
    : public EquationSystems,
      public FEMServiceAbstract<ElementTypeFactory>

{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : EquationSystems(std::forward<Args>(args)...), FEMServiceAbstract<ElementTypeFactory>(nullptr)
    {
        element_builder = std::make_shared<ElementBuilder>(nullptr, this);
    }

    std::shared_ptr<ElementBuilder> get_element_builder() const
    {
        return element_builder;
    }

private:
    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};

#endif // APEC_M_SYSTEM_OBJECTS_H
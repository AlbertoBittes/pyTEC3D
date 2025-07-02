#ifndef APEC_M_MODEL_OBJECTS_H
#define APEC_M_MODEL_OBJECTS_H

#include <memory>
#include <type_traits> // Add this include for std::is_base_of

#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "fem_service_abstract.h"
#include "element_builder.h"
#include "assembler.h"
#include "boundary_conditions_object.h"

using namespace libMesh;

enum SystemType
{
    MAIN = 0,
    NODAL_EXPLICIT,
    ELEMENT_EXPLICIT,
    ELEMENT_LOADS,
    // NODAL_BCS,
    NUM_SYSTEMS
};

// Map SystemType to system name
inline const char *get_system_name(SystemType type)
{
    switch (type)
    {
    case MAIN:
        return "ApecMLinearSystem";
    case NODAL_EXPLICIT:
        return "NodalSolutionSystem";
    case ELEMENT_EXPLICIT:
        return "ElementSolutionSystem";
    case ELEMENT_LOADS:
        return "ElementLoadsSystem";
    // case NODAL_BCS:
    //     return "NodalBCsSystem";
    default:
        return "";
    }
}

// Map SystemType to VariableCategory
inline Globals::VariableCategory get_variable_category(SystemType type)
{
    switch (type)
    {
    case MAIN:
        return Globals::VariableCategory::DOFVariables;
    case NODAL_EXPLICIT:
        return Globals::VariableCategory::NodalSolutionVariables;
    case ELEMENT_EXPLICIT:
        return Globals::VariableCategory::ElementSolutionVariables;
    case ELEMENT_LOADS:
        return Globals::VariableCategory::LoadsVariables;
    // case NODAL_BCS:
    //     return NodalBCsVariables;
    default:
        return Globals::VariableCategory::DOFVariables; // fallback
    }
}

// Map VariableCategory to SystemType (if needed)
inline SystemType get_system_type_from_category(Globals::VariableCategory cat)
{
    switch (cat)
    {
    case Globals::VariableCategory::DOFVariables:
        return MAIN;
    case Globals::VariableCategory::NodalSolutionVariables:
        return NODAL_EXPLICIT;
    case Globals::VariableCategory::ElementSolutionVariables:
        return ELEMENT_EXPLICIT;
    case Globals::VariableCategory::LoadsVariables:
        return ELEMENT_LOADS;
    // case NodalBCsVariables:
    //     return NODAL_BCS;
    default:
        return NUM_SYSTEMS;
    }
}

// Map system name to SystemType
inline SystemType get_system_type_from_name(const std::string &name)
{
    if (name == "ApecMLinearSystem")
        return MAIN;
    if (name == "NodalSolutionSystem")
        return NODAL_EXPLICIT;
    if (name == "ElementSolutionSystem")
        return ELEMENT_EXPLICIT;
    if (name == "ElementLoadsSystem")
        return ELEMENT_LOADS;
    // if (name == "NodalBCsSystem")
    //     return NODAL_BCS;
    return NUM_SYSTEMS;
}

template <typename SysLike, typename AssemblerLike>
class ApecMSystem
    : public FEMServiceAbstract<ElementTypeFactory>,
      public SysLike
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

    // Setter for system_type
    void set_system_type(SystemType type)
    {
        system_type = type;
    }

    // Getter for system_type
    SystemType get_system_type() const
    {
        return system_type;
    }

private:
    std::shared_ptr<AssemblerLike> assembler;
    SystemType system_type;
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

    void release_subset_solve()
    {
        if (this->_subset != nullptr)
        {
            delete this->_subset;
            this->_subset = nullptr;
        }
    }

private:
    std::shared_ptr<BoundaryConditionsObject> boundary_conditions_object;
};

typedef ApecMSystem<ExplicitSystem, AssemblerSolution> ApecMSolutionSystem;

class ApecMEquationSystems
    : public FEMServiceAbstract<ElementTypeFactory>,
      public EquationSystems
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

#endif // APEC_M_MODEL_OBJECTS_H
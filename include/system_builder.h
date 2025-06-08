#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include <string>
#include "fe_data_object.h"
#include "fem_service_abstract.h"
#include "globals.h"
#include "element_types.h" // for VariableCategory

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
inline VariableCategory get_variable_category(SystemType type)
{
    switch (type)
    {
    case MAIN:
        return VariableCategory::DOFVariables;
    case NODAL_EXPLICIT:
        return VariableCategory::NodalSolutionVariables;
    case ELEMENT_EXPLICIT:
        return VariableCategory::ElementSolutionVariables;
    case ELEMENT_LOADS:
        return VariableCategory::LoadsVariables;
    // case NODAL_BCS:
    //     return NodalBCsVariables;
    default:
        return VariableCategory::DOFVariables; // fallback
    }
}

// Map VariableCategory to SystemType (if needed)
inline SystemType get_system_type_from_category(VariableCategory cat)
{
    switch (cat)
    {
    case VariableCategory::DOFVariables:
        return MAIN;
    case VariableCategory::NodalSolutionVariables:
        return NODAL_EXPLICIT;
    case VariableCategory::ElementSolutionVariables:
        return ELEMENT_EXPLICIT;
    case VariableCategory::LoadsVariables:
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

// SystemBuilder as a service using shared_ptr to a TypedContextView of ElementTypeFactory
class SystemBuilder
    : public System::Initialization,
      public FEMServiceAbstract<ElementTypeFactory>,
      public FEParentObjectBase<System>
{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Constructor: pass context view as shared_ptr to FEMServiceAbstract
    SystemBuilder(System *sys, std::shared_ptr<ContextView> ctx_view = nullptr)
        : FEParentObjectBase<System>(sys), FEMServiceAbstract<ElementTypeFactory>(ctx_view) {}

    // Override the initialize method
    void initialize() override;

    // Example: access the ElementTypeFactory from the context view
    std::shared_ptr<ElementTypeFactory> element_type_factory() { return this->context<ElementTypeFactory>(); }
};

#endif // SYSTEM_BUILDER_H
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
        return DOFVariables;
    case NODAL_EXPLICIT:
        return NodalSolutionVariables;
    case ELEMENT_EXPLICIT:
        return ElementSolutionVariables;
    case ELEMENT_LOADS:
        return LoadsVariables;
    // case NODAL_BCS:
    //     return NodalBCsVariables;
    default:
        return DOFVariables; // fallback
    }
}

// Map VariableCategory to SystemType (if needed)
inline SystemType get_system_type_from_category(VariableCategory cat)
{
    switch (cat)
    {
    case DOFVariables:
        return MAIN;
    case NodalSolutionVariables:
        return NODAL_EXPLICIT;
    case ElementSolutionVariables:
        return ELEMENT_EXPLICIT;
    case LoadsVariables:
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

// Map system name to VariableCategory
inline VariableCategory get_variable_category_from_name(const std::string &name)
{
    return get_variable_category(get_system_type_from_name(name));
}

class SystemBuilder : public System::Initialization, public FEDataObjectBase<System>, public FEMServiceAbstract
{
public:
    // Constructor: pass context to FEMServiceAbstract
    SystemBuilder(System *sys, std::shared_ptr<FEMDataContext> context = nullptr)
        : FEDataObjectBase<System>(sys), FEMServiceAbstract(context) {}

    // Override the initialize method
    void initialize() override;
};

#endif // SYSTEM_BUILDER_H
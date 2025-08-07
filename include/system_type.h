#ifndef SYSTEM_TYPE_H
#define SYSTEM_TYPE_H

#include <string>

namespace Globals
{
    enum class VariableCategory;
}
// Here some useful enums and helper functions for the systems of equations

// The main system is represented by the magnetic vector potential DOFs (A_x, A_y, etc.). It is
// the system which matrix and rhs are built by Assembler and ElementBuilder classes.

// It is solved to get the nodal degrees of freedom

// The NodalSolutionSystem will hold those quantities derived from magnetic vector potential and other
// DOF variables, associated to mesh nodes, of the main system. Equivalent in APDL: Nodal Solutions in POST1 (B, etc.)

// The ElementSolutionSystem system will hold those quantities derived from magnetic vector potential and other
// DOF variables, associated to mesh elements, of the main system. Equivalent in APDL: Element Solutions in POST1 (J, etc.)

// The ElementLoadsSystem holds the body loads (Js in APDL)
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

#endif // SYSTEM_TYPE_H
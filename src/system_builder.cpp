#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"

using namespace libMesh;

void SystemBuilder::createSystem()
{
    if (transient)
    {
        if (ferromagnetic)
        {
            system_name = "NonlinearTransient";
            // Add TransientSystem<NonlinearImplicitSystem> as System and return its pointer
            equation_systems.add_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
        }
        else
        {
            system_name = "LinearTransient";
            // Add TransientSystem<LinearImplicitSystem> as System and return its pointer
            equation_systems.add_system<TransientSystem<LinearImplicitSystem>>(system_name);
        }
    }
    else
    {
        if (ferromagnetic)
        {
            system_name = "Nonlinear";
            // Add NonlinearImplicitSystem as System and return its pointer
            equation_systems.add_system<NonlinearImplicitSystem>(system_name);
        }
        else
        {
            system_name = "Linear";
            // Add LinearImplicitSystem as System and return its pointer
            equation_systems.add_system<LinearImplicitSystem>(system_name);
        }
        // Assume the rest of the method continues to handle other cases based on the attributes
    }
};

void SystemBuilder::setVariables()
{
    // Assuming "SystemName" is the name of the system you've added in createSystem
    auto &system = equation_systems.get_system<System>(0);

    // Example combinations - replace with actual logic and variable names
    if (fe_family == "family1" && dim == 2 && formulation == "AV")
    {
        system.add_variable("variable1", FIRST, LAGRANGE);
    }
    else if (fe_family == "family1" && dim == 2 && formulation == "form2")
    {
        system.add_variable("variable2", FIRST, LAGRANGE);
    }
    else if (fe_family == "family1" && dim == 3 && formulation == "form1")
    {
        system.add_variable("variable3", SECOND, LAGRANGE);
    }
    else if (fe_family == "family1" && dim == 3 && formulation == "form2")
    {
        system.add_variable("variable4", SECOND, LAGRANGE);
    }
    else if (fe_family == "family2" && dim == 2 && formulation == "form1")
    {
        system.add_variable("variable5", FIRST, MONOMIAL);
    }
    else if (fe_family == "family2" && dim == 3 && formulation == "form2")
    {
        system.add_variable("variable6", SECOND, MONOMIAL);
    }
    // Add more conditions as needed for other combinations
}

template <typename SystemType>
SystemType &SystemBuilder::getSystemByType()
{
    // Ensure SystemType is derived from System
    static_assert(std::is_base_of<System, SystemType>::value, "SystemType must be derived from System");

    if (system_name == "NonlinearTransient")
    {
        return equation_systems.get_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
    }
    else if (system_name == "LinearTransient")
    {
        return equation_systems.get_system<TransientSystem<LinearImplicitSystem>>(system_name);
    }
    else if (system_name == "Nonlinear")
    {
        return equation_systems.get_system<NonlinearImplicitSystem>(system_name);
    }
    else if (system_name == "Linear")
    {
        return equation_systems.get_system<LinearImplicitSystem>(system_name);
    }
    return nullptr;
}
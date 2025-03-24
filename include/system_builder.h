#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include "libmesh/equation_systems.h"
#include "../include/json.hpp"

using namespace libMesh;

class SystemBuilder : public System::Initialization
{
public:
    // Modified constructor to accept additional parameters
    SystemBuilder(EquationSystems &eq_systems)
        : equation_systems(eq_systems)
    {
        createSystem(); // Call createSystem at the end of the constructor
        setVariables(); // Call setVariables at the end of the constructor
    }

    void createSystem();

    void setVariables();

private:
    EquationSystems &equation_systems;
    bool transient;
    bool ferromagnetic;
    // Assuming var1, var2, var3 are placeholders; replace with actual types and names
    std::string system_name;
};

#endif // SYSTEM_BUILDER_H
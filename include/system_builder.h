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
        : equation_systems(eq_systems) {}

    void initialize() override;

private:
    EquationSystems &equation_systems;
};

#endif // SYSTEM_BUILDER_H
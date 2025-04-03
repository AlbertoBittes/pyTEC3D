#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include "libmesh/equation_systems.h"
#include "../include/json.hpp"

using namespace libMesh;

class SystemBuilder : public System::Initialization
{
public:
    // Modified constructor to accept additional parameters
    SystemBuilder(EquationSystems &eq_systems, std::string &system_name)
        : equation_systems(eq_systems)
    {
        this->system_name = system_name;
    }

    void initialize() override;

private:
    EquationSystems &equation_systems;
    std::string system_name;
};

#endif // SYSTEM_BUILDER_H
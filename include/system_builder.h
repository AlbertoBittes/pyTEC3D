#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include "libmesh/equation_systems.h"
#include "../include/json.hpp"

using namespace libMesh;

class SystemBuilder
{
public:
    EquationSystems &equation_systems;
    bool transient;
    bool ferromagnetic;
    // Assuming var1, var2, var3 are placeholders; replace with actual types and names
    std::string fe_family;
    int dim;
    std::string formulation;
    int type_num;
    std::string system_name;
    std::string order;
    const nlohmann::json &config;
    std::shared_ptr<System::Assembly> assembler; // Add AssemblerBase attribute

    // Modified constructor to accept additional parameters
    SystemBuilder(EquationSystems &eq_systems, const nlohmann::json &config)
        : equation_systems(eq_systems), config(config)
    {
        createSystem(); // Call createSystem at the end of the constructor
        setVariables(); // Call setVariables at the end of the constructor
    }

    void createSystem();

    void setVariables();
};

#endif // SYSTEM_BUILDER_H
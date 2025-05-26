#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include <string>
#include "libmesh/equation_systems.h"
#include "../include/fe_data_object.h"
#include "../include/globals.h"
using namespace libMesh;

class SystemBuilder : public System::Initialization, public FEDataObject<Globals::DummyType, EquationSystems>
{
public:
    // Default constructor
    SystemBuilder()
        : FEDataObject<Globals::DummyType, EquationSystems>(nullptr),
          system_name("") {}
    // Constructor with EquationSystems reference and system name
    SystemBuilder(EquationSystems *eq_systems, const std::string &system_name)
        : FEDataObject<Globals::DummyType, EquationSystems>(eq_systems),
          system_name(system_name) {}

    // Setters
    void set_system_name(const std::string &name) { system_name = name; }

    // Getters
    const std::string &get_system_name() const { return system_name; }

    // Override the initialize method
    void initialize() override;

private:
    std::string system_name; // Name of the system
};

#endif // SYSTEM_BUILDER_H
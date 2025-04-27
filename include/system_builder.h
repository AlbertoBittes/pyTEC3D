#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include "libmesh/equation_systems.h"
#include <string>

using namespace libMesh;

class SystemBuilder : public System::Initialization
{
public:
    // Default constructor
    SystemBuilder()
        : equation_systems(nullptr), system_name("") {}

    // Constructor with EquationSystems reference and system name
    SystemBuilder(EquationSystems &eq_systems, const std::string &system_name)
        : equation_systems(&eq_systems), system_name(system_name) {}

    // Setters
    void set_equation_systems(EquationSystems &eq_systems) { equation_systems = &eq_systems; }
    void set_system_name(const std::string &name) { system_name = name; }

    // Getters
    EquationSystems *get_equation_systems() const { return equation_systems; }
    const std::string &get_system_name() const { return system_name; }

    // Override the initialize method
    void initialize() override;

private:
    EquationSystems *equation_systems; // Pointer to the EquationSystems object
    std::string system_name;           // Name of the system
};

#endif // SYSTEM_BUILDER_H
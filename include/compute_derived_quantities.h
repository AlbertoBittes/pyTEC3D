#ifndef COMPUTE_DERIVED_QUANTITIES_H
#define COMPUTE_DERIVED_QUANTITIES_H

#include "libmesh/equation_systems.h"
#include "../include/assembler.h"

class DerivedQuantitiesCalculator
{
public:
    // Constructor that accepts an AssemblerLinear object and an EquationSystems object
    DerivedQuantitiesCalculator(std::shared_ptr<AssemblerLinear> assembler, EquationSystems &equation_systems);
    void compute_derived_quantities(std::string sys_name);

private:
    // Member variables to store the AssemblerLinear and EquationSystems objects
    std::shared_ptr<AssemblerLinear> assembler;
    EquationSystems &equation_systems;
};

#endif // COMPUTE_DERIVED_QUANTITIES_H
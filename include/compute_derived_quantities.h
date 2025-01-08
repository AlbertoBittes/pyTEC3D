#ifndef COMPUTE_MAGNETIC_FIELD_H
#define COMPUTE_MAGNETIC_FIELD_H

#include "libmesh/equation_systems.h"
#include "../include/assembler.h"

class AssemblerDerivedQuantitiesLinear : public AssemblerLinear
{
public:
    void compute_derived_quantities(libMesh::EquationSystems &es);
    void compute_magnetic_derived_quantities_2D(libMesh::EquationSystems &es);
    void compute_magnetic_derived_quantities_3D(libMesh::EquationSystems &es);
};

#endif // COMPUTE_MAGNETIC_FIELD_H
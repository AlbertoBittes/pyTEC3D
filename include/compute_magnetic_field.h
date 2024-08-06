#ifndef COMPUTE_MAGNETIC_FIELD_H
#define COMPUTE_MAGNETIC_FIELD_H

#include "libmesh/equation_systems.h"

void compute_magnetic_field_static(libMesh::EquationSystems &es);

void compute_magnetic_field_static_2D(libMesh::EquationSystems &es);

#endif // COMPUTE_MAGNETIC_FIELD_H
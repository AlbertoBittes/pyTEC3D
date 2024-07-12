#ifndef ASSEMBLE_MATRIX_H
#define ASSEMBLE_MATRIX_H

#include "libmesh/equation_systems.h"

using namespace libMesh;

void assemble_poisson(EquationSystems &es, const std::string &system_name);

void assemble_poisson_3d(EquationSystems &es,
                         const std::string &system_name);

#endif // ASSEMBLE_MATRIX_H
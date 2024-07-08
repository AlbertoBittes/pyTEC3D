#ifndef ASSEMBLE_MATRIX_H
#define ASSEMBLE_MATRIX_H

#include "libmesh/equation_systems.h"

void assemble_poisson(libmesh::EquationSystems &es, const std::string &system_name);

void assemble_poisson_3d(libmesh::EquationSystems &es,
                         const std::string &system_name);

#endif // ASSEMBLE_MATRIX_H
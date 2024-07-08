#ifndef POISSON_H
#define POISSON_H

#include "libmesh/equation_systems.h"

using namespace libMesh;
// Remove the incomplete type 'void'
void assemble_poisson(libMesh::EquationSystems &es, const std::string &system_name);

void assemble_poisson_3d(libMesh::EquationSystems &es,
                         const std::string &libmesh_dbg_var(system_name))

#endif // POISSON_H
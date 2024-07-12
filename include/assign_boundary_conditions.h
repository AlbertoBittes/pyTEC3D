#ifndef ASSIGN_BC_H
#define ASSIGN_BC_H

#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"

using namespace libMesh;

void assign_boundary_conditions(LinearImplicitSystem &system, unsigned int a_var);

void assign_boundary_ids(Mesh &mesh);

#endif // ASSIGN_BC_H
#ifndef NODE_DOF_CONSTRAINTS_H
#define NODE_DOF_CONSTRAINTS_H

#include <libmesh/equation_systems.h>
#include <string>

// Forward declaration from libMesh
using namespace libMesh;

// Declaration of the node_dof_constraints function
void node_dof_constraints(EquationSystems &equation_systems, const std::string &system_name);

#endif // NODE_DOF_CONSTRAINTS_H
#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "loads_object.h"
#include "apec_m_system_objects.h"
#include "element_types_factory.h"

void ElementBuilder::build_infin112()
{
    libmesh_error_msg("build_infin112 is not implemented yet.");
}

void ElementBuilder::build_infin112_rhs()
{
    libmesh_error_msg("build_infin112_rhs is not implemented yet.");
}

void ElementBuilder::build_infin112_nodal_solution_vector()
{
    libmesh_error_msg("build_infin112_nodal_solution_vector is not implemented yet.");
}

void ElementBuilder::build_infin112_element_solution_vector()
{
    libmesh_error_msg("build_infin112_element_solution_vector is not implemented yet.");
}
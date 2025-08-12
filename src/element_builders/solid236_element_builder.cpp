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

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236()
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236_rhs()
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

void ElementBuilder::build_solid236_nodal_solution_vector()
{
    libmesh_error_msg("build_solid236_nodal_solution_vector is not implemented yet.");
}

void ElementBuilder::build_solid236_element_solution_vector()
{
    // This is a placeholder for the actual implementation
    libmesh_error_msg("build_solid236_element_solution_vector is not implemented yet.");
}
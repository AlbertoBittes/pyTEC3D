// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Introduction Example 1 - Creation of a Mesh Object</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the fourth example program.  It builds on
// the third example program by showing how to formulate
// the code in a dimension-independent way.  Very minor
// changes to the example will allow the problem to be
// solved in one, two or three dimensions.
//
// This example will also introduce the PerfLog class
// as a way to monitor your code's performance.  We will
// use it to instrument the matrix assembly code and look
// for bottlenecks where we should focus optimization efforts.
//
// This example also shows how to extend example 3 to run in
// parallel.  Notice how little has changed!

// C++ include files that we need
#include <algorithm>
#include <iostream>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"

#include "../include/assemble_matrix.h"
#include "../include/assign_boundary_conditions.h"
#include "../include/loads_bc_object.h"
#include "../include/session.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

using namespace libMesh;

// Begin the main program.
int main(int argc, char **argv)
{
    // Initialize libMesh and any dependent libraries, like in example 2.
    LibMeshInit init(argc, argv);

    Session session(argv[1]);

    // Declare a performance log for the main program
    // PerfLog perf_main("Main Program");

    // Check for proper calling arguments.
    libmesh_error_msg_if(argc < 2, "Usage:\n"
                                       << "\t " << argv[0] << " config.json");

    // Brief message to the user regarding the program name
    // and command line arguments.
    libMesh::out << "Running " << argv[0];

    for (int i = 1; i < argc; i++)
        libMesh::out << " " << argv[i];

    libMesh::out << std::endl
                 << std::endl;

    const std::string input_filename = session.config["Mesh Input Section"]["mesh_filepath"].get<std::string>();

// We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
    libmesh_example_requires(false, "--enable-dirichlet");
#endif

    // Read FE order from config file
    std::string order = session.config["Finite Elements Settings"]["approx_order"].get<std::string>();
    std::string family = session.config["Finite Elements Settings"]["fe_family"].get<std::string>();
    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());
    mesh.read(input_filename);

    std::string family_type = family;

    unsigned int dim = mesh.mesh_dimension();

    if (dim == 3 && family == "LAGRANGE")
    {
        family_type = "LAGRANGE_VEC";
    }
    // Print information about the mesh to the screen.
    mesh.print_info();
    Mesh bond_mesh(init.comm());

    // Create an equation systems object.

    mesh.get_boundary_info().print_info();
    mesh.get_boundary_info().print_summary();
    mesh.get_boundary_info().sync(bond_mesh);
    bond_mesh.print_info();

    EquationSystems equation_systems(mesh);

    System *system = nullptr;
    int formulation = 1;
    unsigned int a_var;
    // Declare the system and its variables.
    // Create a system named "Poisson"
    if (dim == 2)
    {

        // Add the variable "u" to "Poisson".  "u"
        // will be approximated using second-order approximation by default
        // Read FE Family from command line
        LinearImplicitSystem &system =
            equation_systems.add_system<LinearImplicitSystem>("Poisson");
        a_var = system.add_variable("AZ", Utility::string_to_enum<Order>(order),
                                    Utility::string_to_enum<FEFamily>(family));
        system.attach_assemble_function(assemble_poisson);
    }
    else if (dim == 3)
    {
        // Add the variable "u" to "Poisson".  "u"
        // will be approximated using second-order approximation by default
        a_var = system.add_variable(
            "A", Utility::string_to_enum<Order>(order),
            Utility::string_to_enum<FEFamily>(family)); /* code */
        system.attach_assemble_function(assemble_poisson_3d);
    }

    // Assign boundary conditions
    assign_boundary_conditions(system, a_var);

    // In general, when reusing a system-indexed exact solution, we want
    // to use the default system-ordering constructor for
    // DirichletBoundary, so we demonstrate that here.  In this case,
    // though, we have only one variable, so system- and local-
    // orderings are the same.

    // We must add the Dirichlet boundary condition _before_
    // we call equation_systems.init()

    // Initialize the data structures for the equation system.
    equation_systems.init();

    // Print information about the system to the screen.
    equation_systems.print_info();

    // Solve the system "Poisson", just like example 2.
    system.solve();

    GmshIO(mesh).write_equation_systems(
        session.config["Output"]["output_filepath"].get<std::string>(), equation_systems);

    // All done.
    return 0;
}

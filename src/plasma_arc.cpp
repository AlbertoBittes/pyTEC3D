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
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/gmsh_io.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/zero_function.h"
#include "../include/assemble_matrix.h"
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the EquationSystems object and the
// name of the system we are assembling as input.  From the
// EquationSystems object we have access to the Mesh and
// other objects we might need.

void assemble_poisson(EquationSystems &es,
                      const std::string &system_name);
void assign_boundary_conditions(Mesh *mesh);
// Exact solution function prototype.
// Real exact_solution(const Real x,
//                     const Real y,
//                     const Real z = 0.);

// // Define a wrapper for exact_solution that will be needed below
// void exact_solution_wrapper(DenseVector<Number> &output,
//                             const Point &p,
//                             const Real)
// {
//   output(0) = exact_solution(p(0),
//                              (LIBMESH_DIM > 1) ? p(1) : 0,
//                              (LIBMESH_DIM > 2) ? p(2) : 0);
// }

// Begin the main program.
int main(int argc, char **argv)
{
  int b;
  std::cin >> b;
  // Initialize libMesh and any dependent libraries, like in example 2.
  LibMeshInit init(argc, argv);

  // Declare a performance log for the main program
  // PerfLog perf_main("Main Program");

  // Check for proper calling arguments.
  libmesh_error_msg_if(argc < 4, "Usage:\n"
                                     << "\t " << argv[0] << " -d 2(3) in.mesh");

  // Brief message to the user regarding the program name
  // and command line arguments.
  libMesh::out << "Running " << argv[0];

  for (int i = 1; i < argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl
               << std::endl;

  const std::string input_filename = argv[3];
  // Read problem dimension from command line.  Use int
  // instead of unsigned since the GetPot overload is ambiguous
  // otherwise.
  const int dim = libMesh::command_line_next("-d", 2);

  // We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
  libmesh_example_requires(false, "--enable-dirichlet");
#endif

  // Read FE order from command line. For now 1st order fixed
  std::string order = "FIRST";
  order = libMesh::command_line_next("-o", order);
  order = libMesh::command_line_next("-Order", order);

  // Create a mesh, with dimension to be overridden later, distributed
  // across the default MPI communicator.
  Mesh mesh(init.comm());
  // GmshIO mesh_io(mesh);
  mesh.read(input_filename);

  // Print information about the mesh to the screen.
  mesh.print_info();
  Mesh bond_mesh(init.comm());

  // Create an equation systems object.
  EquationSystems equation_systems(mesh);

  // Declare the system and its variables.
  // Create a system named "Poisson"
  if (dim == 2)
  {
    LinearImplicitSystem &system =
        equation_systems.add_system<LinearImplicitSystem>("Poisson");
    // Add the variable "u" to "Poisson".  "u"
    // will be approximated using second-order approximation by default
    // Read FE Family from command line
    std::string family = "LAGRANGE";
    family = libMesh::command_line_next("-f", family);
    family = libMesh::command_line_next("-FEFamily", family);
    unsigned int a_var = system.add_variable("Az",
                                             Utility::string_to_enum<Order>(order),
                                             Utility::string_to_enum<FEFamily>(family));
  }
  else if (dim == 3)
  {
    LinearImplicitSystem &system =
        equation_systems.add_system<LinearImplicitSystem>("Poisson");
    // Add the variable "u" to "Poisson".  "u"
    // will be approximated using second-order approximation by default
    std::string family = "LAGRANGE_VEC";
    family = libMesh::command_line_next("-f", family);
    family = libMesh::command_line_next("-FEFamily", family);
    unsigned int a_var = system.add_variable("A",
                                             Utility::string_to_enum<Order>(order),
                                             Utility::string_to_enum<FEFamily>(family)); /* code */
  }

  mesh.get_boundary_info().print_info();
  mesh.get_boundary_info().print_summary();
  mesh.get_boundary_info().sync(bond_mesh);
  bond_mesh.print_info();
  ZeroFunction<Number> zero;

  // boundaryInfo.get
  system.attach_assemble_function(assemble_poisson);
  std::set<boundary_id_type> boundary_ids;
  boundary_ids.insert(0);

  DirichletBoundary dirichlet_bc(boundary_ids, {a_var}, zero);

  // In general, when reusing a system-indexed exact solution, we want
  // to use the default system-ordering constructor for
  // DirichletBoundary, so we demonstrate that here.  In this case,
  // though, we have only one variable, so system- and local-
  // orderings are the same.

  // We must add the Dirichlet boundary condition _before_
  // we call equation_systems.init()
  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();
  mesh.print_info();

  // Solve the system "Poisson", just like example 2.
  system.solve();

  GmshIO(mesh).write_equation_systems("/home/bitteal/PlasmaArc/simulations/out_2.e", equation_systems);

  // All done.
  return 0;
}

// Give the system a pointer to the matrix assembly
// function.
// system.attach_assemble_function(assemble_poisson);

// Construct a Dirichlet boundary condition object

//   // Indicate which boundary IDs we impose the BC on
//   // We either build a line, a square or a cube, and
//   // here we indicate the boundaries IDs in each case
//   std::set<boundary_id_type> boundary_ids;
//   // the dim==1 mesh has two boundaries with IDs 0 and 1
//   boundary_ids.insert(0);
//   boundary_ids.insert(1);
//   // the dim==2 mesh has four boundaries with IDs 0, 1, 2 and 3
//   if (dim >= 2)
//   {
//     boundary_ids.insert(2);
//     boundary_ids.insert(3);
//   }
//   // the dim==3 mesh has four boundaries with IDs 0, 1, 2, 3, 4 and 5
//   if (dim == 3)
//   {
//     boundary_ids.insert(4);
//     boundary_ids.insert(5);
//   }

// #ifdef LIBMESH_ENABLE_DIRICHLET
//   // In general, when reusing a system-indexed exact solution, we want
//   // to use the default system-ordering constructor for
//   // DirichletBoundary, so we demonstrate that here.  In this case,
//   // though, we have only one variable, so system- and local-
//   // orderings are the same.
//   // DirichletBoundary dirichlet_bc(boundary_ids, {u_var}, exact_solution_object);

//   // We must add the Dirichlet boundary condition _before_
//   // we call equation_systems.init()
//   // system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
// #endif

//   // Initialize the data structures for the equation system.
//   equation_systems.init();

//   // Print information about the system to the screen.
//   equation_systems.print_info();
//   mesh.print_info();

//   // Solve the system "Poisson", just like example 2.
//   system.solve();

// #ifdef LIBMESH_HAVE_EXODUS_API

//   ExodusII_IO(mesh).write_equation_systems((dim == 3) ? "out_3.e" : "out_2.e", equation_systems);

// #endif // #ifdef LIBMESH_HAVE_EXODUS_API

//   // All done.
//   return 0;
// }

// Define a new function that receives as input a pointer to mesh object and returns void
void assign_boundary_conditions(Mesh *mesh)
{
}
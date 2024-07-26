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
#include <filesystem>
// Basic include file needed for the mesh functionality.
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/node.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// To impose Dirichlet boundary conditions
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "../include/loads_bc_mat_object.h"
#include "../include/session.h"
#include "../include/system_builder.h"
#include "../include/global_resources.h"
#include "../include/node_dof_constraints.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/transient_system.h"
#include "libmesh/vtk_io.h"

// Include the necessary header files for vtkDataArray and vtkPointData
#include <vtkUnstructuredGrid.h>
#include "vtkDataArray.h"
#include "vtkPointData.h"

using namespace libMesh;

// Begin the main program.
int main(int argc, char **argv)
{
    // int a;
    // std::cin >> a;
    // Initialize libMesh and any dependent libraries, like in example 2.
    LibMeshInit init(argc, argv);

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
    libMesh::out << "Reading PlasmaArc configurations..." << std::endl;

    Session session(argv[1]);

    bool transient = session.config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
    bool ferromagnetic = session.config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
    std::string fe_family = session.config["Finite Elements Settings"]["fe_family"].get<std::string>();
    std::string formulation = session.config["Finite Elements Settings"]["formulation"].get<std::string>();
    std::string order = session.config["Finite Elements Settings"]["approx_order"].get<std::string>();
    std::string inp_name = session.config["Mesh Input Section"]["input_name"].get<std::string>();
    std::string loads_folderpath = session.config["Loads and BC Files Section"]["loads_folderpath"].get<std::string>();
    std::string bc_folderpath = session.config["Loads and BC Files Section"]["bc_folderpath"].get<std::string>();
    // Declare a performance log for the main program
    // PerfLog perf_main("Main Program");
    std::cout << "Reading mesh..." << std::endl;

// We use Dirichlet boundary conditions here
#ifndef LIBMESH_ENABLE_DIRICHLET
    libmesh_example_requires(false, "--enable-dirichlet");
#endif

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());
    VTKIO vtkio(mesh);
    // mesh.allow_renumbering(false);
    std::string filname = session.config["Mesh Input Section"]["mesh_filepath"].get<std::string>();
    vtkio.read(filname);
    vtkUnstructuredGrid *vtkGrid = vtkio.get_vtk_grid();
    vtkPointData *pointData = vtkGrid->GetPointData();

    // Step 4: Retrieve the Array
    vtkDataArray *ansysNodeNumArray = pointData->GetArray("ansys_node_num");

    if (ansysNodeNumArray)
    {
        // Process the array as needed
        std::cout << "Found array 'ansys_node_num' with "
                  << ansysNodeNumArray->GetNumberOfTuples() << " tuples." << std::endl;
    }
    else
    {
        std::cerr << "Array 'ansys_node_num' not found." << std::endl;
    }
    unsigned int dim = mesh.mesh_dimension();

    // Print information about the mesh to the screen.
    mesh.print_info();

    // Create an equation systems object.
    libMesh::out << "Boundary Information:" << std::endl;

    mesh.get_boundary_info().print_info();

    mesh.get_boundary_info().print_summary();

    EquationSystems equation_systems(mesh);

    unsigned int a_var;
    // Declare the system and its variables.
    // Create a system named "Poisson"

    SystemBuilder system_builder(equation_systems, transient, ferromagnetic,
                                 fe_family, dim, formulation, order);

    std::string system_name = "Linear";
    // Assign boundary  conditions
    System &system = equation_systems.get_system(system_name);
    if (transient)
    {
        if (ferromagnetic)
        {
            system_name = "NonlinearTransient";
            // Add TransientSystem<NonlinearImplicitSystem> as System and return its pointer
            TransientSystem<NonlinearImplicitSystem> *system = nullptr;
            system = &equation_systems.get_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
        }
        else
        {
            system_name = "LinearTransient";
            // Add TransientSystem<LinearImplicitSystem> as System and return its pointer
            TransientSystem<LinearImplicitSystem> *system = nullptr;
            system = &equation_systems.get_system<TransientSystem<LinearImplicitSystem>>(system_name);
        }
    }
    else
    {
        if (ferromagnetic)
        {
            system_name = "Nonlinear";
            // Add NonlinearImplicitSystem as System and return its pointer
            NonlinearImplicitSystem *system = nullptr;
            system = &equation_systems.get_system<NonlinearImplicitSystem>(system_name);
        }
        else
        {
            system_name = "Linear";
            // Add LinearImplicitSystem as System and return its pointer
            // LinearImplicitSystem *system = nullptr;
            // system = &equation_systems.get_system<LinearImplicitSystem>(system_name);
        }
        // Assume the rest of the method continues to handle other cases based on the attributes
    }

    LoadsBCMatObject loads_bc_object(loads_folderpath, bc_folderpath, inp_name);
    GlobalResources::setLoadsBCObject(loads_bc_object);

    // In general, when reusing a system-indexed exact solution, we want
    // to use the default system-ordering constructor for
    // DirichletBoundary, so we demonstrate that here.  In this case,
    // though, we have only one variable, so system- and local-
    // orderings are the same.

    // We must add the Dirichlet boundary condition _before_
    // we call equation_systems.init()
    system.attach_constraint_function(node_dof_constraints);
    // Initialize the data structures for the equation system.
    equation_systems.init();
    // Redo this with Dirichlet boundaries

    // equation_systems.get_system<LinearImplicitSystem>("Linear").get_dof_map().process_constraints(equation_systems.get_system<LinearImplicitSystem>("Linear").get_mesh());
    // equation_systems.get_system<LinearImplicitSystem>("Linear").get_dof_map().prepare_send_list();
    // Print information about the system to the screen.
    equation_systems.print_info();

    // Solve the system "Poisson", just like example 2.
    system.solve();

    // Write the results
    vtkio.write_equation_systems(session.config["Output"]["output_filepath"].get<std::string>(), equation_systems);

    // All done.
    return 0;
}

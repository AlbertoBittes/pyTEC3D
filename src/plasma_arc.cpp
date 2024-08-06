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
#include "libmesh/libmesh.h"
#include "libmesh/node.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h" // Add this line to include the BoundaryInfo header

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "../include/loads_bc_mat_object.h"
#include "../include/session.h"
#include "../include/system_builder.h"
#include "../include/global_resources.h"
#include "../include/node_dof_constraints.h"
#include "../include/compute_magnetic_field.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/transient_system.h"
#include "libmesh/vtk_io.h"

// Include the necessary header files for vtkDataArray and vtkPointData
#include <vtkUnstructuredGrid.h>
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include <vtkSmartPointer.h>

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

    std::string mesh_filepath = session.config["Mesh Input Section"]["mesh_filepath"].get<std::string>();
    std::string inp_name = std::filesystem::path(mesh_filepath).stem().string(); // Extract basename without extension
    // Declare a performance log for the main program
    // PerfLog perf_main("Main Program");
    libMesh::out << "Reading mesh..." << std::endl;

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());
    VTKIO vtkio(mesh);
    // mesh.allow_renumbering(false);
    vtkio.read(mesh_filepath);

    vtkUnstructuredGrid *vtkGrid = vtkio.get_vtk_grid();

    libMesh::out << "Preparing Mesh..." << std::endl;

    unsigned char dim = 0;

    double bounds[6];
    vtkGrid->GetBounds(bounds);

    // Check x-dimension
    if (bounds[1] - bounds[0] > 0)
    {
        dim++;
    }

    // Check y-dimension
    if (bounds[3] - bounds[2] > 0)
    {
        dim++;
    }

    // Check z-dimension
    if (bounds[5] - bounds[4] > 0)
    {
        dim++;
    }

    mesh.set_spatial_dimension(dim);
    mesh.prepare_for_use();

    // Print information about the mesh to the screen.
    mesh.print_info();

    libMesh::out << "Building equations systems..." << std::endl;

    EquationSystems equation_systems(mesh);

    // Declare the system and its variables.
    // Create a system named "Poisson"

    vtkCellData *cellData = vtkGrid->GetCellData();
    if (cellData == nullptr)
    {
        libmesh_error_msg("Error: vtkGrid does not contain Cell Data.");
    }
    SystemBuilder system_builder(equation_systems, session.config, cellData);

    std::string system_name = "Linear";
    // Assign boundary  conditions
    LinearImplicitSystem &system = equation_systems.get_system<LinearImplicitSystem>(system_name);

    // LoadsBCMatObject loads_bc_object(loads_folderpath, bc_folderpath, inp_name);
    // GlobalResources::setLoadsBCObject(loads_bc_object);

    // Assign boundary conditions
    libMesh::out << "Setting boundary conditions..." << std::endl;
    vtkPointData *pointData = vtkGrid->GetPointData();
    if (pointData == nullptr)
    {
        libmesh_error_msg("Error: vtkGrid does not contain Point Data.");
    }
    NodalBCConstraint node_dof_constraints(system, pointData);
    system.attach_constraint_object(node_dof_constraints);

    // Initialize the data structures for the equation system.
    equation_systems.init();

    // Print information about the system to the screen.
    equation_systems.print_info();

    // Build the K matrix, the RHS F and Solve the system Kx=F
    system.solve();

    compute_magnetic_field_static_2D(equation_systems);
    // Use filesystem to modify the basename
    std::filesystem::path output_path(mesh_filepath);
    output_path.replace_filename(output_path.stem().string() + "_out" + output_path.extension().string());
    // Write the results
    vtkio.write_equation_systems(output_path, equation_systems);

    // All done.
    return 0;
}

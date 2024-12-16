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

// C++ include files that we need
#include <algorithm>
#include <iostream>
#include <math.h>
#include <set>
#include <filesystem>
// Basic include file needed for the mesh functionality.
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h" // Add this line to include the BoundaryInfo header

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

#include "../include/session.h"
#include "../include/system_builder.h"
#include "../include/node_dof_constraints.h"
#include "../include/compute_magnetic_field.h"
#include "../include/materials_properties.h"
#include "../include/time_stepper.h" // Add this line to include the TimeStepper header
#include "libmesh/transient_system.h"
#include "libmesh/vtk_io.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/numeric_vector.h"

// Include the necessary header files for vtkDataArray and vtkPointData
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <iostream>

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

    std::string output_folderpath;
    if (session.config["Output Section"]["output_filepath"].is_null())
    {
        output_folderpath = std::filesystem::path(mesh_filepath).parent_path().string();
    }
    else
    {
        output_folderpath = session.config["Output Section"]["output_filepath"].get<std::string>();
    }

    std::string inp_name = std::filesystem::path(mesh_filepath).stem().string(); // Extract basename without extension
    // Declare a performance log for the main program
    // PerfLog perf_main("Main Program");
    libMesh::out << "Reading mesh..." << std::endl;

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());
    VTKIO vtkio(mesh);
    mesh.allow_renumbering(false);
    if (init.comm().rank() == 0)
        vtkio.read(mesh_filepath);

    MeshCommunication().broadcast(mesh);

    vtkUnstructuredGrid *vtkGrid = vtkio.get_vtk_grid();

    // This is needed because of wron libMesh function when reading dimensions
    // Ask if they plan to fix it
    // Plus, check because this is not working with MPI,
    //     Hmm... I think if you do it this way, the file is opened and read from simultaneously on all N cores, which is probably not great for the filesystem. Can you try wrapping the mesh.read() line in if (comm().rank() == 0), and then, on all processors call

    // MeshCommunication().broadcast(mesh);
    libMesh::out << "Checking mesh dimensions..." << std::endl;
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

    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh.prepare_for_use();

    libMesh::out << "Reading loads data..." << std::endl;
    vtkCellData *cell_data = vtkGrid->GetCellData();
    if (cell_data == nullptr)
    {
        libmesh_error_msg("Error: vtkGrid does not contain any Cell Data.");
    }

    // Here I read material properties and assign conductive and non-conductive subdomains
    libMesh::out << "Reading materials properties..." << std::endl;
    vtkDataArray *matArray = cell_data->GetArray("ansys_material_type");
    if (matArray == nullptr)
    {
        libmesh_error_msg("Error: 'ansys_material_type' array not found in vtkCellData.");
    }
    unsigned int matNum;
    MaterialsMap materials = readMaterialsProperties(session.config["Mesh Input Section"]["materials_filepath"].get<std::string>());

    libMesh::out << "Identifying conductive domain..." << std::endl;
    for (auto &elem : mesh.active_element_ptr_range())
    {
        matNum = matArray->GetTuple1(elem->id());
        auto it = materials.find(matNum);
        if (it == materials.end())
        {
            libmesh_error_msg("Error: Material number " + std::to_string(matNum) + " not found in material properties.");
        }
        if (std::isnan(it->second.rsvx))
        {
            elem->subdomain_id() = 0;
        }
        else
        {
            elem->subdomain_id() = 1;
        }
    }
    // Print information about the mesh to the screen.
    mesh.print_info();

    libMesh::out << "Building equations systems..." << std::endl;

    EquationSystems equation_systems(mesh);

    // Declare the system and its variables
    SystemBuilder system_builder(equation_systems, session.config, cell_data, materials);

    libMesh::out << "Distributing degrees of freedom..." << std::endl;

    // Get a reference to the system we want to solve
    TransientLinearImplicitSystem &system = equation_systems.get_system<TransientLinearImplicitSystem>("PlasmaArcLinearSystem");

    // Initialize the data structures for the equation system.
    equation_systems.init();

    // Print information about the system to the screen.
    equation_systems.print_info();

    equation_systems.parameters.set<MaterialsMap>("materials_properties") = materials;

    // Setup transient options
    libMesh::out << "Setting up transient analysis options..." << std::endl;
    // Initializing variables for transient analysis options
    bool compute_equilibrium_first_step = false;
    bool transient = session.config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
    TimeStepper timeStepper(session.config, system);
    equation_systems.parameters.set<TimeStepper>("time_stepper") = timeStepper;
    // Iterate through the time steps

    // Extract the time_steps array
    if (transient)
    {
        // Check if the user wants to compute the equilibrium first step
        compute_equilibrium_first_step = session.config["Transient Analysis Settings"]["equilibrium_step"].get<bool>();
        // Extract the time_steps array
        // The list of time steps doesn't include the initial time, t=0 always is assumed
        // A static initial step may be required by the user however to get equilibrium
        system.time = 0;
    }
    else
    {
        // For static analysis, we set the time to 0
        system.time = 0;
    }

    // Get Point Data from the vtkGrid, it contains the boundary conditions to be applied to mesh nodes
    libMesh::out << "Reading boundary conditions data..." << std::endl;
    vtkPointData *pointData = vtkGrid->GetPointData();
    if (pointData == nullptr)
    {
        libmesh_error_msg("Error: vtkGrid does not contain Point Data.");
    }
    // unsigned int first_idx;
    // unsigned int last_idx;
    // for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
    // {
    //     libMesh::out << "Array " << i << ": " << pointData->GetArrayName(i) << std::endl;
    //     system.add_vector(pointData->GetArrayName(i));
    //     NumericVector<Number> &new_vec = system.get_vector(pointData->GetArrayName(i));
    //     first_idx = new_vec.first_local_index();
    //     last_idx = new_vec.last_local_index();
    //     for(unsigned int j = first_idx; j <= last_idx; ++j)
    //     {
    //         new_vec.set(j, pointData->GetArray(i)->GetComponent(j, 0));
    //     }
    // }
    libMesh::out << "Setting output path..." << std::endl;
    // Use filesystem to modify the basename
    std::filesystem::path output_path(output_folderpath);
    std::filesystem::path mesh_path(mesh_filepath);
    std::filesystem::path final_output_path;
    std::string mesh_filename_without_extension = mesh_path.stem().string();
    std::string output_filename;

    // Assign boundary conditions
    libMesh::out << "Setting boundary conditions object..." << std::endl;
    NodalBCConstraint node_dof_constraints(system, pointData, 0);
    system.attach_constraint_object(node_dof_constraints);

    libMesh::out << "Check for equilibrium first step..." << std::endl;
    // Compute first step for equilibrium or static analysis
    if (!transient || compute_equilibrium_first_step)
    {
        if (compute_equilibrium_first_step)
        {
            libMesh::out << "Computing equilibrium first step..." << std::endl;
        }
        else
        {
            libMesh::out << "Computing static analysis..." << std::endl;
        }

        libMesh::out << "Reinitializing constrained DOFs..." << std::endl;
        // This will reinitialize the constraints for the system
        system.reinit_constraints();

        // Build the K matrix, the RHS F and Solve the system Kx=F
        libMesh::out << "Building matrix and solving system..." << std::endl;
        system.solve();

        libMesh::out << "Computing derived quantities..." << std::endl;
        // Call the function to compute the magnetic field
        if (dim == 2)
            compute_magnetic_field_static_2D(equation_systems);
        else
            compute_magnetic_field_static_3D(equation_systems);

        libMesh::out << "Writing results..." << std::endl;
        // Construct the output file path
        output_filename = mesh_filename_without_extension + "_out_0";
        final_output_path = output_path / output_filename;
        // Write the results
        vtkio.write_equation_systems(final_output_path, equation_systems);
    }

    for (auto &step : timeStepper)
    {
        libMesh::out << "Solving time step " << step.get_step() << std::endl;

        libMesh::out << "Time increment at the end of time step: " << step.get_delta_t() << std::endl;

        libMesh::out << "Number of time substeps: " << step.get_n_substeps() << std::endl;

        libMesh::out << "Setting boundary conditions..." << std::endl;

        node_dof_constraints.setTimeStep(step.get_step() + 1);
        system.reinit_constraints();
        // Iterate through the substeps of the current time step
        for (auto &substep : step)
        {
            libMesh::out << "  Substep: " << substep.get_index() << ", Time: " << substep.get_time()
                         << ", dt: " << substep.get_dt()
                         << ", Write Results: " << substep.get_write_results()
                         << ", Transient Effect: " << substep.get_transient_effect()
                         << ", Ramped Load: " << substep.get_ramped_load() << std::endl;
            // Assemble & solve the linear system
            libMesh::out << "Building matrix and solving system..." << std::endl;
            system.solve();
        }
    }
    for (unsigned int time_idx = 0; time_idx < time_steps.size(); ++time_idx)
    {
        // Increment the time counter, set the time and the
        // time step size as parameters in the EquationSystem.
        dt = time_steps[time_idx] - system.time;
        dt_substep = dt / n_substeps[time_idx];

        libMesh::out << "Solving time step " << time_idx + 1 << " of " << time_steps.size() << std::endl;

        libMesh::out << "Time increment at the end of time step: " << dt << std::endl;

        libMesh::out << "Number of time substeps: " << n_substeps[time_idx] << std::endl;

        libMesh::out << "Setting boundary conditions..." << std::endl;

        for (unsigned int substep = 0; substep < n_substeps[time_idx]; ++substep)
        {
            equation_systems.parameters.set<Real>("dt") = dt_substep;
            system.time += dt_substep;
            // A pretty update message
            libMesh::out << "Solving time substep " << substep + 1 << " of " << n_substeps[time_idx] << std::endl;
            libMesh::out << "Time increment: " << dt_substep << std::endl;
            libMesh::out << "Time: " << system.time << std::endl;
            // At this point we need to update the old
            // solution vector.  The old solution vector
            // will be the current solution vector from the
            // previous time step. We will do this by extracting the
            // system from the EquationSystems object and using
            // vector assignment. Since only TransientSystems
            // (and systems derived from them) contain old solutions
            // we need to specify the system type when we ask for it.
            libMesh::out << "Updating old solution..." << std::endl;
            *system.old_local_solution = *system.current_local_solution;

            // Here do 2 things: ramp the load and assemble the matrix only if first substep

            if (write_substeps_results[time_idx] || substep == n_substeps[time_idx] - 1)
            {
                libMesh::out << "Computing derived quantities..." << std::endl;
                // Call the function to compute the magnetic field
                if (dim == 2)
                    compute_magnetic_field_static_2D(equation_systems);
                else
                    compute_magnetic_field_static_3D(equation_systems);
                libMesh::out << "Time substep " << substep + 1 << " of " << n_substeps[time_idx] << " completed." << std::endl;
                libMesh::out << "Writing time substep results..." << std::endl;
                output_filename = mesh_filename_without_extension + "_out_" + std::to_string(cumulative_n_substeps + 1);
                final_output_path = output_path / output_filename;
                // Write the results
                vtkio.write_equation_systems(final_output_path, equation_systems);
            }
            else
            {
                libMesh::out << "Time substep " << substep + 1 << " of " << n_substeps[time_idx] << " completed." << std::endl;
            }
            cumulative_n_substeps += 1;
        }

        libMesh::out << "Time step " << time_idx + 1 << " of " << time_steps.size() << " completed." << std::endl;
    }
    // All done.
    return 0;
}

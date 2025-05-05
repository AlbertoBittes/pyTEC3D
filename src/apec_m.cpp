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
// #include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h" // Add this line to include the BoundaryInfo header

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

#include "../include/element_builder.h"
#include "../include/assembler.h"
#include "../include/session.h"
#include "../include/system_builder.h"
#include "../include/element_types.h" // Add this line to include the ElementTypesMap header
#include "../include/loads_object.h"
#include "../include/time_step.h" // Add this line to include the TimeStepper header
#include "libmesh/transient_system.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/communicator.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/dof_map.h"
#include "libmesh/solver_configuration.h"
#include "../include/node_bc_object.h" // Include the base class

#include "../include/solver_configuration.h"
#include "../include/vtkio_updated.h" // Ensure this header file defines VTKIOUpdated
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/periodic_boundary.h" // Add this line to include the PeriodicBoundary header
#include "../include/azimuthal_periodic_boundary.h"

using namespace libMesh;

EquationSystems *build_equation_systems(Mesh &mesh, Session &session)
{
    // libMesh::out << "Building equations systems..." << std::endl;
    EquationSystems *equation_systems = new EquationSystems(mesh);
    // equation_systems.parameters.set<std::set<subdomain_id_type> *>("conductive_subdomain") = &conductive_subdomain;
    equation_systems->parameters.set<Session *>("session") = &session;

    // equation_systems->parameters.set<std::shared_ptr<TimeStepper>>("time_stepper") = std::make_shared<TimeStepper>();

    return equation_systems;
}

// void read_parse_element_types(std::string element_types_filepath, EquationSystems &equation_systems)
// {
//     libMesh::out << "Reading element types..." << std::endl;
//     equation_systems.parameters.set<std::shared_ptr<ElementTypesCollector>>("element_types_collector") = std::make_shared<ElementTypesCollector>(element_types_filepath);
// }

void init_equation_systems(EquationSystems &equation_systems)
{
    // Declare the system and its variables
    libMesh::out << "Building system of equations..." << std::endl;

    std::string system_name = "ApecMLinearSystem";
    std::string der_system_name = "DerivedQuantitiesSystem";
    std::string der_elem_system_name = "DerivedQuantitiesSystemElem";
    // Add TransientLinearImplicitSystem as System and return its pointer
    libMesh::out << "Checking system options..." << std::endl;

    equation_systems.clear();

    // read_parse_element_types(equation_systems.parameters.get<Session *>("session")->config["Mesh Input Section"]["element_types_filepath"].get<std::string>(), equation_systems);

    equation_systems.add_system<TransientLinearImplicitSystem>(system_name);
    TransientLinearImplicitSystem &system = equation_systems.get_system<TransientLinearImplicitSystem>(system_name);
    ExplicitSystem &der_system = equation_systems.add_system<ExplicitSystem>(der_system_name);
    ExplicitSystem &der_elem_system = equation_systems.add_system<ExplicitSystem>(der_elem_system_name);

    SystemBuilder system_builder(equation_systems, system_name);
    system.attach_init_object(system_builder);

    SystemBuilder system_builder_der(equation_systems, der_system_name);
    der_system.attach_init_object(system_builder_der);

    SystemBuilder system_builder_der_elem(equation_systems, der_elem_system_name);
    der_elem_system.attach_init_object(system_builder_der_elem);

    // This will go in /solu
    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    equation_systems.init();
    equation_systems.print_info();

    system.zero_out_matrix_and_rhs = false;
    NodeBCObject *node_bc_object = new NodeBCObject();
    system.attach_constraint_object(*node_bc_object);
    // system.solution->zero();
    // system.solution->close();
    system.current_local_solution->zero();
    system.current_local_solution->close();
}

void read_and_prepare_mesh(Mesh &mesh, std::string mesh_filepath)
{
    libMesh::out << "Reading mesh..." << std::endl;
    // mesh.clear();
    VTKIO_Updated vtkio(mesh);
    if (mesh.comm().rank() == 0)
    {
        vtkio.read(mesh_filepath);
    }
    mesh.comm().broadcast(*(vtkio.get_element_types_map_csv()));
    std::shared_ptr<ElementTypesCollector> element_types_collector = std::make_shared<ElementTypesCollector>(*(vtkio.get_element_types_map_csv()));
    MeshCommunication().broadcast(mesh);
    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh.prepare_for_use();
    MeshTools::Modification::orient_elements(mesh);
    // Print information about the mesh to the screen.
    std::cout << "Mesh Information: " << std::endl;
    mesh.print_info();
}

void solve(EquationSystems &equation_systems, LoadStep &load_step)
{
    auto time_step = load_step.get_time_step();

    libMesh::out << "Solving time step " << time_step->get_step() + 1 << std::endl;

    libMesh::out << "Time increment at the end of time step: " << time_step->get_delta_t() << std::endl;

    libMesh::out << "Number of time substeps: " << time_step->get_n_substeps() << std::endl;

    libMesh::out << "Setting boundary conditions..." << std::endl;

    libMesh::out << "Reinitializing constrained DOFs..." << std::endl;

    auto session = equation_systems.parameters.get<Session *>("session");
    // read material properties

    // Read loads

    // Get the system
    TransientLinearImplicitSystem &system = equation_systems.get_system<TransientLinearImplicitSystem>("ApecMLinearSystem");
    // This will reinitialize the constraints for the system
    // Read boundary conditions

    // system.get_dof_map().print_dof_constraints();
    // Iterate through the substeps of the current time step
    for (auto &substep : *time_step)
    {
        libMesh::out << "Substep: " << substep.get_index() + 1 << ", Time: " << substep.get_time()
                     << ", dt: " << substep.get_dt()
                     << ", Write Results: " << substep.get_write_results()
                     << ", Transient Effect: " << substep.get_transient_effect() << std::endl;
        // update time in element data object
        // Assemble & solve the linear system
        libMesh::out << "Building matrix and solving system..." << std::endl;
        system.solve();

        system.get_linear_solver()->print_converged_reason();
        if (substep.get_write_results() || substep.get_index() == time_step->get_n_substeps() - 1)
        {
            libMesh::out << "Computing derived quantities..." << std::endl;
            // Call the function to compute the magnetic field

            // equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->compute_derived_quantities();
            // init.comm().barrier();
            // equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystem").solution->print(std::cout);
            // libMesh::out << "Writing time substep results..." << std::endl;
            // output_filename = mesh_filename_without_extension + "_out_" + std::to_string(timeStepper->get_cumulative_substep());
            // output_filename = output_filename + ".pvtu";
            // final_output_path = output_path / output_filename;
            std::set<std::string> system_names = {"ApecMLinearSystem", "DerivedQuantitiesSystem"};
            // // Write the results
            // VTKIO vtk_output(mesh);
            // exo_output.write_equation_systems(final_output_path, equation_systems, &system_names);
            libMesh::out << "Writing time substep results..." << std::endl;
            // session->set_output_filepath(time_step.get_time_stepper().get_cumulative_substep());
            // // Write the results
            VTKIO_Updated vtk_output(equation_systems.get_mesh(), &equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystemElem"));
            // vtk_output.write_equation_systems(session->final_output_path, equation_systems, &system_names);
            // change the extension of filepath to .exo
            // std::string exo_output_path = session.final_output_path;
            // exo_output_path.replace(exo_output_path.end() - 4, exo_output_path.end(), ".exo");
            // ExodusII_IO exo_output(mesh);
            // exo_output.write_timestep(exo_output_path, equation_systems, 1, substep.get_time(), &system_names);
            // exo_output.write_element_data(equation_systems);

            // Print to a file the rhs
            // std::ofstream rhs_file("../simulations/two_wires_3D/matrix/RHS_full_" + std::to_string(step.get_step()) + ".txt");
            // system.rhs->print_global(rhs_file);
            // system.rhs->zero();
            // system.rhs->close();
            // output_filename = mesh_filename_without_extension + "_out";
            // final_output_path = output_path / output_filename;
            // std::ostringstream file_name;
            // // // We write the file in the ExodusII format.
            // file_name << final_output_path.string() << "_"
            //           << std::setw(3)
            //           << std::setfill('0')
            //           << std::right
            //           << timeStepper->get_cumulative_substep();

            // // // Write the results
            // Nemesis_IO nem_output(mesh);
            // std::vector<std::string> var_names;
            // unsigned int var_num = system.n_vars();
            // for (unsigned int i = 0; i < var_num; ++i)
            // {
            //     var_names.push_back(system.variable_name(i));
            // }
            // ExplicitSystem &derived_system = equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystem");
            // var_num = derived_system.n_vars();
            // for (unsigned int i = 0; i < var_num; ++i)
            // {
            //     var_names.push_back(derived_system.variable_name(i));
            // }
            // nem_output.set_output_variables(var_names);

            // nem_output.write_timestep(file_name.str(), equation_systems, 1, substep.get_time());
            // // Empty var_names
            // nem_output.write_element_data(equation_systems);
        }

        libMesh::out << "Time substep " << substep.get_index() + 1 << " of " << time_step->get_n_substeps() << " completed." << std::endl;
    }
}

// Function to execute the batch mode logic
void run_batch_mode(const std::string &mesh_filename, libMesh::Parallel::Communicator &comm)
{
    libMesh::out << "Running in batch mode with mesh file: " << mesh_filename << std::endl;

    // Create a session with the mesh filename
    Session session(mesh_filename);

    Mesh mesh(comm);

    EquationSystems &equation_systems = *(build_equation_systems(mesh, session));

    // Read and prepare the mesh
    read_and_prepare_mesh(mesh, session.get_input_mesh_filepath());

    // Build and initialize equation systems
    init_equation_systems(equation_systems);

    std::shared_ptr<LoadStepper> load_stepper = std::make_shared<LoadStepper>();
    // Solve the system for each time step
    for (unsigned int n = 0; n < session.get_n_load_steps(); n++)
    {
        // Create a new load step
        LoadStep load_step(n, equation_systems, session);
        load_step.populate_time_step();

        // Solve the system for the current load step
        solve(equation_systems, load_step);
    }
    {
        solve(equation_systems, step);
    }

    // Clean up
    // delete &mesh;
    // delete &equation_systems;
}

// Main function
int main(int argc, char **argv)
{
    // Initialize libMesh and any dependent libraries
    LibMeshInit init(argc, argv);

    // Check for proper calling arguments
    if (argc < 2)
    {
        libmesh_error_msg("Usage:\n"
                          << "\t" << argv[0] << " -b mesh_filename\n"
                          << "\t" << argv[0] << " -i");
    }

    // Parse command-line arguments
    std::string mode = argv[1];
    if (mode == "-b")
    {
        if (argc < 3)
        {
            libmesh_error_msg("Error: Missing mesh filename. Usage:\n"
                              << "\t" << argv[0] << " -b mesh_filename");
        }

        // Get the mesh filename
        std::string mesh_filename = argv[2];

        // Run batch mode
        run_batch_mode(mesh_filename, init.comm());
    }
    else if (mode == "-i")
    {
        libMesh::out << "Running in interactive mode..." << std::endl;

        // Interactive mode logic
        libMesh::out << "Interactive mode is not yet implemented." << std::endl;
    }
    else
    {
        libmesh_error_msg("Error: Invalid mode. Usage:\n"
                          << "\t" << argv[0] << " -b mesh_filename\n"
                          << "\t" << argv[0] << " -i");
    }

    // All done
    return 0;
}

// This will be done much later

// for (auto &elem : mesh.active_local_element_ptr_range())
// {
//     for ( auto &side : elem->boundary_sides_begin)
// }
// BoundaryMesh boundary_mesh(mesh.comm());

// boundaryInfo.print_summary();
// // std::map<dof_id_type, dof_id_type> node_id_map;
// // std::map<dof_id_type, unsigned char> side_id_map;
// // boundary_mesh.print_info();
// // boundaryInfo.get_side_and_node_maps(boundary_mesh, node_id_map, side_id_map);
// // iterate over the boundary mesh
// for (auto &elem : boundary_mesh.active_local_element_ptr_range())
// {

// }
// Put all this part in a function, it must also add an array of bools called conductive to make it easier to identify conductive elements
// init.comm().set_union(elem_types);
// std::map<subdomain_id_type, subdomain_id_type> subdomain_map;
// subdomain_id_type subdomain_id = 0;
// for (const auto &elem_type : elem_types)
// {
//     subdomain_map[elem_type] = subdomain_id;
//     subdomain_id++;
// }
// BoundaryInfo &boundaryInfo = mesh.get_boundary_info(); // Retrieves the BoundaryInfo object
// // boundaryInfo.sync(boundary_mesh);
// {
//     Point center(0., 0., 0.);
//     Point axis(0., 0., 1.);
//     Real angle = -90 * 2 * libMesh::pi / 360;
//     AzimuthalPeriodicBoundary periodic_bc(center, axis, angle);
//     periodic_bc.myboundary = 7;
//     periodic_bc.pairedboundary = 8;
//     equation_systems.get_system<TransientLinearImplicitSystem>(0).get_dof_map().add_periodic_boundary(periodic_bc);
// }
// vtkDataArray *orig_id = partitionedGrid->GetPointData()->GetArray("VTKorigID");
// vtkDataArray *boundary_id = partitionedGrid->GetPointData()->GetArray("boundary_id");
// Print the names of the arrays in the vtkPointData
// for (int i = 0; i < partitionedGrid->GetPointData()->GetNumberOfArrays(); i++)
// {
//     libMesh::out << "Array " << i << ": " << partitionedGrid->GetPointData()->GetArrayName(i) << std::endl;
// }
// if (!boundary_id)
// {
//     libmesh_error_msg("Error: 'boundary_id' array not found in vtkPointData.");
// }

// vtkIdType size_vtk = orig_id->GetNumberOfTuples();
// dof_id_type orig_id_num;
// const Node *node_ptr = nullptr;
// unsigned int val;
// for (vtkIdType j = 0; j < size_vtk; j++)
// {
//     orig_id_num = static_cast<dof_id_type>(orig_id->GetTuple1(j));
//     node_ptr = mesh.node_ptr(orig_id_num);
//     val = boundary_id->GetTuple1(j);
//     if (val != 0)
//         boundaryInfo.add_node(node_ptr, boundary_id->GetTuple1(j));
// }
// boundaryInfo.build_side_list_from_node_list();
// boundaryInfo.clear_boundary_node_ids();
// boundaryInfo.print_info();
// boundaryInfo.print_summary();
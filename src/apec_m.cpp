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
#include "../include/node_dof_constraints.h"
#include "../include/materials_properties.h"
#include "../include/element_types.h" // Add this line to include the ElementTypesMap header
#include "../include/element_data_manager.h"
#include "../include/time_stepper.h" // Add this line to include the TimeStepper header
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

    equation_systems->parameters.set<std::shared_ptr<TimeStepper>>("time_stepper") = std::make_shared<TimeStepper>();

    return equation_systems;
}

void read_parse_element_types(std::string element_types_filepath, EquationSystems &equation_systems)
{
    libMesh::out << "Reading element types..." << std::endl;
    equation_systems.parameters.set<std::shared_ptr<ElementTypesCollector>>("element_types_collector") = std::make_shared<ElementTypesCollector>(element_types_filepath);
}

void init_equation_systems(EquationSystems &equation_systems)
{
    // Declare the system and its variables
    libMesh::out << "Building system of equations..." << std::endl;

    std::string system_name = "ApecMLinearSystem";
    std::string der_system_name = "DerivedQuantitiesSystem";
    std::string der_elem_system_name = "DerivedQuantitiesSystemElem";
    // Add TransientLinearImplicitSystem as System and return its pointer
    libMesh::out << "Checking system options..." << std::endl;
    // Here I will implement various checks about the system element types, the mesh and everything
    // ferromagnetic = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
    // fe_family = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["fe_family"].get<std::string>();
    // formulation = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["formulation"].get<std::string>();
    // order = equation_systems.parameters.get<Session *>("session")->config["Finite Elements Settings"]["approx_order"].get<std::string>();
    // // dim = equation_systems.get_mesh().mesh_dimension();
    // dim = 3;
    // // Check if dim is 2 or 3
    // if (dim != 2 && dim != 3)
    // {
    //     libmesh_error_msg("Error: dim must be 2 or 3.");
    // }
    // // Check if fe_family is not "LAGRANGE" and dim is not 3
    // if (fe_family != "LAGRANGE" && fe_family != "NEDELEC_ONE")
    // {
    //     libmesh_error_msg("Error: Element type must be 'LAGRANGE' (nodal formulation) or 'NEDELEC_ONE' (edge formulation).");
    // }
    // // Check if fe_family is not "LAGRANGE" and dim is not 3
    // if (fe_family != "LAGRANGE" && dim != 3)
    // {
    //     libmesh_error_msg("Error: Mesh must be 3D for edge elements.");
    // }
    // // Check if order is "FIRST" or "SECOND"
    // if (order != "FIRST" && order != "SECOND")
    // {
    //     libmesh_error_msg("Error: order must be 'FIRST' or 'SECOND'.");
    // }
    // // Check if formulation is "AV" or "Hphi"
    // if (formulation != "AV" && formulation != "Hphi")
    // {
    //     libmesh_error_msg("Error: formulation must be 'AV' or 'Hphi'.");
    // }

    // libMesh::out << "Creating system..." << std::endl;
    // if (ferromagnetic)
    // {
    //     system_name = "NonlinearTransient";
    //     // Add TransientSystem<NonlinearImplicitSystem> as System and return its pointer
    //     equation_systems.add_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
    // }
    // else
    // {
    equation_systems.clear();

    read_parse_element_types(equation_systems.parameters.get<Session *>("session")->config["Mesh Input Section"]["element_types_filepath"].get<std::string>(), equation_systems);

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

    equation_systems.parameters.set<std::shared_ptr<ElementDataManager>>("element_data_manager") = std::make_shared<ElementDataManager>(ElementDataManager());

    equation_systems.parameters.get<std::shared_ptr<ElementDataManager>>("element_data_manager")->build_load_data_vectors(equation_systems.get_mesh(), equation_systems.parameters.get<std::shared_ptr<ElementTypesCollector>>("element_types_collector")->getLoadsLabels());

    system.zero_out_matrix_and_rhs = false;
    system.solution->zero();
    system.solution->close();
    system.current_local_solution->zero();
    system.current_local_solution->close();

    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());
    // SolverType solver_type = Utility::string_to_enum<SolverType>("GMRES");
    // solver->set_solver_type(GMRES);
    PetscSolverConfiguration solver_configuration(*solver);
    solver_configuration.configure_solver();
    solver->set_solver_configuration(solver_configuration);
}

Mesh *read_and_prepare_mesh(std::string mesh_filepath, libMesh::Parallel::Communicator &comm)
{
    libMesh::out << "Reading mesh..." << std::endl;
    Mesh *mesh = new Mesh(comm);
    // mesh.clear();
    if (mesh->comm().rank() == 0)
    {
        VTKIO_Updated vtkio(*mesh);
        vtkio.read(mesh_filepath);
    }
    MeshCommunication().broadcast(*mesh);
    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh->prepare_for_use();
    MeshTools::Modification::orient_elements(*mesh);
    // Print information about the mesh to the screen.
    std::cout << "Mesh Information: " << std::endl;
    mesh->print_info();
    return mesh;
}

void solve(EquationSystems &equation_systems, TimeStep &time_step)
{
    libMesh::out << "Solving time step " << time_step.get_step() + 1 << std::endl;

    libMesh::out << "Time increment at the end of time step: " << time_step.get_delta_t() << std::endl;

    libMesh::out << "Number of time substeps: " << time_step.get_n_substeps() << std::endl;

    libMesh::out << "Setting boundary conditions..." << std::endl;

    libMesh::out << "Reinitializing constrained DOFs..." << std::endl;

    auto session = equation_systems.parameters.get<Session *>("session");
    // read material properties
    equation_systems.parameters.get<std::shared_ptr<ElementDataManager>>("element_data_manager")->materials = std::make_shared<MaterialsMap>(*(readMaterialsProperties(equation_systems.parameters.get<Session>("session").config["Mesh Input Section"]["materials_filepath"].get<std::string>())));

    // Read loads
    std::string loads_path = session->mesh_filename_without_extension + "_loads_" + std::to_string(time_step.get_step() + 1);

    equation_systems.parameters.get<std::shared_ptr<ElementDataManager>>("element_data_manager")->read_update_loads_vectors(loads_path);
    // Get the system
    TransientLinearImplicitSystem &system = equation_systems.get_system<TransientLinearImplicitSystem>("ApecMLinearSystem");
    // This will reinitialize the constraints for the system
    // Read boundary conditions
    NodalBCConstraint *node_dof_constraints = dynamic_cast<NodalBCConstraint *>(&system.get_constraint_object());
    std::string bc_path = session->mesh_filename_without_extension + "_bc_" + std::to_string(time_step.get_step() + 1);
    node_dof_constraints->read_update_boundary_conditions(bc_path, time_step);

    // system.get_dof_map().print_dof_constraints();
    // Iterate through the substeps of the current time step
    for (auto &substep : time_step)
    {
        libMesh::out << "Substep: " << substep.get_index() + 1 << ", Time: " << substep.get_time()
                     << ", dt: " << substep.get_dt()
                     << ", Write Results: " << substep.get_write_results()
                     << ", Transient Effect: " << substep.get_transient_effect() << std::endl;
        // update time in element data manager
        equation_systems.parameters.get<std::shared_ptr<ElementDataManager>>("element_data_manager")->update_time_step(substep);
        // Assemble & solve the linear system
        libMesh::out << "Building matrix and solving system..." << std::endl;
        system.solve();

        system.get_linear_solver()->print_converged_reason();
        if (substep.get_write_results() || substep.get_index() == time_step.get_n_substeps() - 1)
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
            session->set_output_filepath(time_step.get_time_stepper().get_cumulative_substep());
            // // Write the results
            VTKIO_Updated vtk_output(equation_systems.get_mesh(), &equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystemElem"));
            vtk_output.write_equation_systems(session->final_output_path, equation_systems, &system_names);
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

        libMesh::out << "Time substep " << substep.get_index() + 1 << " of " << time_step.get_n_substeps() << " completed." << std::endl;
    }
}

// Begin the main program.
int main(int argc, char **argv)
{
    // Initialize libMesh and any dependent libraries, like in example 2.
    LibMeshInit init(argc, argv);
    // {
    //     if (init.comm().rank() != 0)
    //     {
    //         int i = 0;
    //         while (0 == i)
    //             sleep(1);
    //     }
    // }
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

    libMesh::out << "Reading ApecM configurations..." << std::endl;

    Session session(argv[1]);
    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.

    // Declare a performance log for the main program
    // PerfLog perf_main("Main Program");
    Mesh &mesh = *(read_and_prepare_mesh(session.input_mesh_filepath, init.comm()));
    // mesh.allow_renumbering(false);

    // Read element types
    EquationSystems &equation_systems = *(build_equation_systems(mesh, session));

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

    // Use a pointer to the local variable

    // Setup transient options
    libMesh::out << "Setting up transient analysis options..." << std::endl;
    // Initializing variables for transient analysis options

    // Iterate through the time steps

    // Extract the time_steps array
    // Check if the user wants to compute the equilibrium first step
    // Extract the time_steps array
    // The list of time steps doesn't include the initial time, t=0 always is assumed
    // A static initial step may be required by the user however to get equilibrium

    // Get Point Data from the vtkGrid, it contains the boundary conditions to be applied to mesh nodes
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

    // Assign boundary conditions

    // Compute first step for equilibrium or static analysis
    // if (compute_equilibrium_first_step)
    // {
    //     libMesh::out << "Computing equilibrium first step..." << std::endl;
    // }
    // else
    // {
    //     libMesh::out << "Computing static analysis..." << std::endl;
    // }
    // Build the K matrix, the RHS F and Solve the system Kx=F
    // This should be already done

    // solver->set_preconditioner_type(Utility::string_to_enum<PreconditionerType>("ILU"));
    for (auto &step : *equation_systems.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper"))
    {
        solve(equation_systems, step);
    }
    // All done.
    return 0;
}

// std::vector<unsigned int> cell_data_array_sizes;
// Loop over all cell data arrays and add their names to the vector
// vtkCellData *cell_data = vtkGrid->GetCellData();
// for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
// {
//     const char *array_name = cell_data->GetArrayName(i);
//     if (array_name)
//     {
//         // Check if the array names contains "Js"
//         if (std::string(array_name).find("Js_") != std::string::npos)
//         {
//             mesh.add_elem_datum<Real *>(array_name);
//         }
//         else
//         {
//             mesh.add_elem_integer(array_name);
//         }
//     }

// Repeat for point data
// std::vector<std::string> point_data_int_array_names;
// std::string array_name;
// vtkPointData *point_data = vtkGrid->GetPointData();
// for (int i = 0; i < point_data->GetNumberOfArrays(); ++i)
// {
//     array_name = point_data->GetArrayName(i);
//     if (!(std::string(array_name).find("A_") != std::string::npos || std::string(array_name).find("V_") != std::string::npos))
//     {
//         point_data_int_array_names.push_back(array_name);
//     }
// }

// std::vector<std::string> cell_data_int_array_names;
// vtkCellData *cell_data = vtkGrid->GetCellData();
// for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
// {
//     array_name = cell_data->GetArrayName(i);
//     if (!(std::string(array_name).find("A_") != std::string::npos || std::string(array_name).find("V_") != std::string::npos))
//     {
//         cell_data_int_array_names.push_back(array_name);
//     }
// }
// mesh.add_elem_integers(point_data_int_array_names);
// mesh.add_elem_data<Real>(point_data_float_array_names);
// vtkDataArray *data_array;
// for (const auto &elem : mesh.active_local_element_ptr_range())
// {
//     for (unsigned int i = 0; i < cell_data_int_array_names.size(); ++i)
//     {
//         data_array = cell_data->GetArray(cell_data_int_array_names[i].c_str());
//         elem->set_extra_integer(i, data_array->GetTuple1(elem->id()));
//     }
//     for (unsigned int i = 0; i < cell_data_float_array_names.size(); ++i)
//     {
//         data_array = cell_data->GetArray(cell_data_float_array_names[i].c_str());

//         elem->set_extra_datum(i, data_array->GetTuple(elem->id()));
//     }
// }
// for (const auto &node : mesh.local_node_ptr_range())
// {
//     for (unsigned int i = 0; i < point_data_int_array_names.size(); ++i)
//     {
//         data_array = point_data->GetArray(point_data_int_array_names[i].c_str());
//         node->set_extra_integer(i, data_array->GetTuple1(node->id()));
//     }
//     for (unsigned int i = 0; i < point_data_float_array_names.size(); ++i)
//     {
//         data_array = point_data->GetArray(point_data_float_array_names[i].c_str());
//         node->set_extra_datum(i, data_array->GetTuple(node->id()));
//     }
// }
// Node &node = mesh.node_ref(18);
// Real *real_val = node.get_extra_datum<Real *>(5);
// vtkIntArray *int_array_vtk;

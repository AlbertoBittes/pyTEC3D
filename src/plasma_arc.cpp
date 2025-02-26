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
#include <libmesh/exodusII_io.h>
#include <libmesh/nemesis_io.h>
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
#include "../include/time_stepper.h" // Add this line to include the TimeStepper header
#include "libmesh/transient_system.h"
#include "libmesh/vtk_io.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/dof_map.h"
#include "libmesh/solver_configuration.h"
#include "../include/solver_configuration.h"
#include "../include/vtkio_updated.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/enum_solver_type.h"
#include "libmesh/enum_solver_package.h"

// Include the necessary header files for vtkDataArray and vtkPointData
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>
#include <vtkExtractCells.h>

using namespace libMesh;

// vtkSmartPointer<vtkUnstructuredGrid> read_mesh_vtk(Mesh &mesh, const std::string &mesh_filepath)
// {
//     vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

//     if (mesh.comm().rank() == 0)
//     {
//         VTKIO vtkio(mesh);
//         vtkio.read(mesh_filepath);
//         vtkGrid = vtkio.get_vtk_grid();
//     }
//     else
//     {
//         vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();

//         // Pass the filename along to the reader
//         reader->SetFileName(mesh_filepath.c_str());

//         // Force reading
//         reader->Update();

//         // read in the grid
//         vtkGrid = reader->GetOutput();
//     }
//     return vtkGrid;
// }

// vtkSmartPointer<vtkUnstructuredGrid> partition_vtk_grid(Mesh &mesh, vtkSmartPointer<vtkUnstructuredGrid> vtkGrid, std::map<subdomain_id_type, subdomain_id_type> *subdomain_map_pointer)
// {
//     std::set<subdomain_id_type> elem_types;

//     vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
//     vtkSmartPointer<vtkExtractCells> extractCells = vtkSmartPointer<vtkExtractCells>::New();
//     for (const auto &elem : mesh.active_local_element_ptr_range())
//     {
//         elem_types.insert(static_cast<subdomain_id_type>(elem->type()));
//         cellIds->InsertNextId(elem->id());
//     }
//     extractCells->SetInputData(vtkGrid);
//     extractCells->SetCellList(cellIds);
//     extractCells->Update();
//     vtkSmartPointer<vtkUnstructuredGrid> partitionedGrid = extractCells->GetOutput();

//     mesh.comm().set_union(elem_types);
//     subdomain_id_type subdomain_id = 0;
//     for (const auto &elem_type : elem_types)
//     {
//         (*subdomain_map_pointer)[elem_type] = subdomain_id;
//         subdomain_id++;
//     }
//     return partitionedGrid;
// }

// std::set<subdomain_id_type> extract_conductive_subdomain(Mesh &mesh, vtkCellData *cell_data, std::map<subdomain_id_type, subdomain_id_type> &subdomain_map, MaterialsMap &materials)
// {
//     // Here I read material properties and assign conductive and non-conductive subdomains
//     libMesh::out << "Reading materials properties..." << std::endl;
//     vtkDataArray *matArray = cell_data->GetArray("ansys_material_type");
//     if (matArray == nullptr)
//     {
//         libmesh_error_msg("Error: 'ansys_material_type' array not found in vtkCellData.");
//     }
//     unsigned int matNum;

//     // Use a pointer to the local variable
//     MaterialsMap *materials_ptr = &materials;

//     // Also here you can set the subdomain_id for the conductive elements, just choose
//     libMesh::out << "Identifying conductive domain..." << std::endl;
//     unsigned int k = 0;
//     std::set<subdomain_id_type> conductive_subdomain;
//     subdomain_id_type elem_type;
//     // Do this before partitioning the mesh
//     for (auto &elem : mesh.active_local_element_ptr_range())
//     {
//         elem_type = static_cast<int>(elem->type());
//         matNum = matArray->GetTuple1(k);
//         auto it = materials_ptr->find(matNum);
//         if (it == materials_ptr->end())
//         {
//             libmesh_error_msg("Error: Material number " + std::to_string(matNum) + " not found in material properties.");
//         }
//         if (std::isnan(it->second.rsvx))
//         {
//             elem->subdomain_id() = subdomain_map[elem_type];
//         }
//         else
//         {
//             elem->subdomain_id() = subdomain_map[elem_type] + 100;
//             conductive_subdomain.insert(subdomain_map[elem_type] + 100);
//         }
//         k++;
//     }
//     mesh.comm().set_union(conductive_subdomain);
//     return conductive_subdomain;
// }
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
    libMesh::out << "Reading PlasmaArc configurations..." << std::endl;

    Session session(argv[1]);

    // Declare a performance log for the main program
    // PerfLog perf_main("Main Program");
    libMesh::out << "Reading mesh..." << std::endl;

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());
    mesh.allow_renumbering(false);

    vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (init.comm().rank() == 0)
    {
        VTKIO vtkio(mesh);
        vtkio.read(session.mesh_filepath);
        vtkGrid = vtkio.get_vtk_grid();
    }
    else
    {
        vtkSmartPointer<vtkXMLPUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();

        // Pass the filename along to the reader
        reader->SetFileName(session.mesh_filepath.c_str());

        // Force reading
        reader->Update();

        // read in the grid
        vtkGrid = reader->GetOutput();
    }
    MeshCommunication().broadcast(mesh);

    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh.prepare_for_use();
    MeshTools::Modification::orient_elements(mesh);

    // Put all this part in a function, it must also add an array of bools called conductive to make it easier to identify conductive elements
    vtkSmartPointer<vtkUnstructuredGrid> partitionedGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    // std::set<subdomain_id_type> elem_types;
    {
        vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkExtractCells> extractCells = vtkSmartPointer<vtkExtractCells>::New();
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            // elem_types.insert(static_cast<subdomain_id_type>(elem->type()));
            cellIds->InsertNextId(elem->id());
        }
        extractCells->SetInputData(vtkGrid);
        extractCells->SetCellList(cellIds);
        extractCells->Update();
        partitionedGrid = extractCells->GetOutput();
    }
    // init.comm().set_union(elem_types);
    // std::map<subdomain_id_type, subdomain_id_type> subdomain_map;
    // subdomain_id_type subdomain_id = 0;
    // for (const auto &elem_type : elem_types)
    // {
    //     subdomain_map[elem_type] = subdomain_id;
    //     subdomain_id++;
    // }
    vtkGrid = nullptr;

    libMesh::out << "Reading loads data..." << std::endl;
    vtkCellData *cell_data = partitionedGrid->GetCellData();
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

    // Use a pointer to the local variable
    MaterialsMap *materials_ptr = &materials;

    // Also here you can set the subdomain_id for the conductive elements, just choose
    libMesh::out << "Identifying conductive domain..." << std::endl;
    unsigned int k = 0;
    // std::set<subdomain_id_type> conductive_subdomain;
    // subdomain_id_type elem_type;
    // Do this before partitioning the mesh
    // conductive_subdomain.insert(1);

    for (auto &elem : mesh.active_local_element_ptr_range())
    {
        // elem_type = static_cast<int>(elem->type());
        matNum = matArray->GetTuple1(k);
        auto it = materials_ptr->find(matNum);
        if (it == materials_ptr->end())
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
        k++;
    }
    // init.comm().set_union(conductive_subdomain);

    // Print information about the mesh to the screen.
    std::cout << "Mesh Information: " << std::endl;
    mesh.print_info();

    libMesh::out << "Building equations systems..." << std::endl;
    EquationSystems equation_systems(mesh);
    // equation_systems.parameters.set<std::set<subdomain_id_type> *>("conductive_subdomain") = &conductive_subdomain;
    equation_systems.parameters.set<Session *>("session") = &session;
    // Declare the system and its variables
    SystemBuilder system_builder(equation_systems);

    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    // Get a reference to the system we want to solve
    TransientLinearImplicitSystem &system = equation_systems.get_system<TransientLinearImplicitSystem>("PlasmaArcLinearSystem");
    system.zero_out_matrix_and_rhs = false;
    // Initialize the data structures for the equation system.
    equation_systems.init();
    // Print information about the system to the screen.
    equation_systems.print_info();

    // Setup transient options
    libMesh::out << "Setting up transient analysis options..." << std::endl;
    // Initializing variables for transient analysis options
    auto timeStepper = std::make_shared<TimeStepper>(session.config, system);
    equation_systems.parameters.set<std::shared_ptr<TimeStepper>>("time_stepper") = timeStepper;
    auto element_data_manager = std::make_shared<ElementDataManager>(materials, cell_data, timeStepper);
    equation_systems.parameters.set<std::shared_ptr<ElementDataManager>>("element_data_manager") = element_data_manager;
    // Iterate through the time steps

    // Extract the time_steps array
    // Check if the user wants to compute the equilibrium first step
    // Extract the time_steps array
    // The list of time steps doesn't include the initial time, t=0 always is assumed
    // A static initial step may be required by the user however to get equilibrium

    // Get Point Data from the vtkGrid, it contains the boundary conditions to be applied to mesh nodes
    libMesh::out << "Reading boundary conditions data..." << std::endl;
    vtkPointData *pointData = partitionedGrid->GetPointData();
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

    // Assign boundary conditions
    libMesh::out << "Setting boundary conditions object..." << std::endl;
    NodalBCConstraint node_dof_constraints(system, pointData, timeStepper);
    system.attach_constraint_object(node_dof_constraints);

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
    // solver->set_preconditioner_type(Utility::string_to_enum<PreconditionerType>("ILU"));
    for (auto &step : *timeStepper)
    {
        libMesh::out << "Solving time step " << step.get_step() + 1 << std::endl;

        libMesh::out << "Time increment at the end of time step: " << step.get_delta_t() << std::endl;

        libMesh::out << "Number of time substeps: " << step.get_n_substeps() << std::endl;

        libMesh::out << "Setting boundary conditions..." << std::endl;

        libMesh::out << "Reinitializing constrained DOFs..." << std::endl;
        // This will reinitialize the constraints for the system
        element_data_manager->update_time_step();

        if (step.get_step() > 0)
        {
            node_dof_constraints.check_boundary_conditions();
        }
        if (node_dof_constraints.update_boundary_conditions)
        {
            system.reinit_constraints();
            equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->updated_bcs = true;
        }
        else
        {
            equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->updated_bcs = false;
        }
        // system.get_dof_map().print_dof_constraints();
        // Iterate through the substeps of the current time step
        for (auto &substep : step)
        {
            libMesh::out << "Substep: " << substep.get_index() + 1 << ", Time: " << substep.get_time()
                         << ", dt: " << substep.get_dt()
                         << ", Write Results: " << substep.get_write_results()
                         << ", Transient Effect: " << substep.get_transient_effect() << std::endl;
            // Assemble & solve the linear system
            libMesh::out << "Building matrix and solving system..." << std::endl;
            system.solve();
            // // How many iterations were required to solve the linear system?
            // const unsigned int n_linear_iterations = system.n_linear_iterations();

            // // What was the final residual of the linear system?
            // const Real final_linear_residual = system.final_linear_residual();

            // If the solver did no work (sometimes -ksp_converged_reason
            // says "Linear solve converged due to CONVERGED_RTOL
            // iterations 0") but the nonlinear residual norm is above
            // the tolerance, we need to pick an even lower linear
            // solver tolerance and try again.  Note that the tolerance
            // is relative to the norm of the RHS, which for this
            // particular problem does not go to zero, since we are
            // solving for the full solution rather than the update.
            //
            // Similarly, if the solver did no work and this is the 0th
            // nonlinear step, it means that the delta between solutions
            // is being inaccurately measured as "0" since the solution
            // did not change.  Decrease the tolerance and try again.
            // while (system.n_linear_iterations() == 0)
            // {
            //     Real old_linear_solver_tolerance = equation_systems.parameters.get<Real>("linear solver tolerance");
            //     equation_systems.parameters.set<Real>("linear solver tolerance") = 1.e-3 * old_linear_solver_tolerance;
            //     system.solve();
            // }
            solver->print_converged_reason();
            if (substep.get_write_results() || substep.get_index() == step.get_n_substeps() - 1)
            {
                libMesh::out << "Computing derived quantities..." << std::endl;
                // Call the function to compute the magnetic field

                equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->compute_derived_quantities();
                // init.comm().barrier();
                // equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystem").solution->print(std::cout);
                // libMesh::out << "Writing time substep results..." << std::endl;
                // output_filename = mesh_filename_without_extension + "_out_" + std::to_string(timeStepper->get_cumulative_substep());
                // output_filename = output_filename + ".pvtu";
                // final_output_path = output_path / output_filename;
                std::set<std::string> system_names = {"PlasmaArcLinearSystem", "DerivedQuantitiesSystem"};
                // // Write the results
                // VTKIO vtk_output(mesh);
                // exo_output.write_equation_systems(final_output_path, equation_systems, &system_names);
                libMesh::out << "Writing time substep results..." << std::endl;
                session.set_output_filepath(timeStepper->get_cumulative_substep());
                // // Write the results
                VTKIO_Updated vtk_output(mesh, &equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystemElem"));
                vtk_output.write_equation_systems(session.final_output_path, equation_systems, &system_names);
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

            libMesh::out << "Time substep " << substep.get_index() + 1 << " of " << step.get_n_substeps() << " completed." << std::endl;
        }
    }
    // All done.
    return 0;
}

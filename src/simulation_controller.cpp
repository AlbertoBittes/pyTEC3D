#include "simulation_controller.h"
#include "vtkio_updated.h"
#include "element_types.h"
#include "system_builder.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/numeric_vector.h"

SimulationController::SimulationController(libMesh::Parallel::Communicator &comm, bool batch_mode)
    : comm(comm), batch_mode(batch_mode)
{
    session = std::make_shared<Session>();
    mesh = std::make_shared<Mesh>(comm);
    equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
    equation_systems->init_element_builder();
    load_stepper = std::make_shared<LoadStepper>();
    // Initialize other services as needed
}

void SimulationController::run_cli(libMesh::Parallel::Communicator &comm)
{
    std::string input;
    while (true)
    {
        std::cout << "Enter command: ";
        std::getline(std::cin, input);

        Command cmd = parse_command(input);
        switch (cmd)
        {
        case Command::READ_MESH:
        {
            std::string filepath = extract_argument(input);
            session->set_input_mesh_filepath(filepath);
            read_mesh();
            break;
        }
        case Command::CLEAR:
            clear();
            break;
        case Command::PHYSICS:
            physics();
            break;
        case Command::SOLVE:
        {
            std::string step = extract_argument(input);
            load_stepper->read_new_load_step(step);
            solve();
            break;
        }
        case Command::FINISH:
            finish();
            break;
        case Command::EXIT:
            exit();
            return;
        default:
            libMesh::out << "Unknown command.\n";
        }
    }
}

void SimulationController::read_mesh()
{

    libMesh::out << "Reading mesh..." << std::endl;
    std::string mesh_filepath = session->get_input_mesh_filepath();
    // mesh.clear();
    VTKIO_Updated vtkio(*mesh);
    if (mesh->comm().rank() == 0)
    {
        vtkio.read(mesh_filepath);
    }

    auto element_types_map = vtkio.get_element_types_map_csv();
    mesh->comm().broadcast<ElementTypesMapCSV>(*element_types_map);
    std::shared_ptr<ElementTypeFactory> element_types = std::make_shared<ElementTypeFactory>(*element_types_map);
    element_types->create_element_types(*element_types_map);

    // Add element_types to mesh_data_input repository
    mesh_input_data->set_object(element_types);
    element_types = std::make_shared<ElementTypeFactory>(*element_types_map);
    element_types->create_element_types(*element_types_map);

    MeshCommunication().broadcast(*mesh);

    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh->prepare_for_use();
    MeshTools::Modification::orient_elements(*mesh);

    // Print information about the mesh to the screen.
    std::cout << "Mesh Information: " << std::endl;
    mesh->print_info();
}

void SimulationController::physics()
{
    libMesh::out << "Building system of equations..." << std::endl;

    System *systems[NUM_SYSTEMS] = {nullptr};

    libMesh::out << "Checking system options..." << std::endl;

    auto system_builder_view = std::make_shared<
        TypedContextView<ElementTypeFactory>>(
        TypedContextView<ElementTypeFactory>::from(
            mesh_input_data));
    for (int i = 0; i < NUM_SYSTEMS; ++i)
    {
        SystemType type = static_cast<SystemType>(i);
        const std::string sys_name = get_system_name(type);

        if (type == MAIN)
            systems[i] = &equation_systems->add_system<ApecMLinearSystem>(sys_name);
        else
            systems[i] = &equation_systems->add_system<ExplicitSystem>(sys_name);

        SystemBuilder builder(systems[i], system_builder_view);
        systems[i]->attach_init_object(builder);
    }

    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    equation_systems->init();
    equation_systems->print_info();
}

void SimulationController::solve()
{
    auto load_step = load_stepper->get_active_load_step();

    auto prev_load_step = load_stepper->get_old_load_step();

    auto time_step = load_step->get_object<TimeStep>();

    libMesh::out << "Solving time step " << time_step->get_step() + 1 << std::endl;

    libMesh::out << "Time increment at the end of time step: " << time_step->get_delta_t() << std::endl;

    libMesh::out << "Number of time substeps: " << time_step->get_n_substeps() << std::endl;

    libMesh::out << "Setting boundary conditions..." << std::endl;

    libMesh::out << "Reinitializing constrained DOFs..." << std::endl;
    // read material properties

    // Read loads
    // Here you have to prepare the various contexts for the systems
    // Get the system
    ApecMLinearSystem &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    // Read boundary conditions
    // 1. prepare the context for element builder
    std::shared_ptr<ElementBuilderContext> element_builder_context =
        std::make_shared<ElementBuilderContext>(element_types,
                                                load_step->get_object<MaterialsObject>(),
                                                load_step->get_object<RealConstantsObject>(),
                                                load_step->get_object<LoadsObject>(),
                                                session);
    auto element_builder = equation_systems->get_element_builder();
    element_builder->set_context(element_builder_context);

    // 2. prepare context for assembler
    std::shared_ptr<AssemblerContext> assembler_context = std::make_shared<AssemblerContext>(time_step, session);
    auto assembler = system.get_assembler();
    assembler->set_context(assembler_context);

    // 3. Prepare context for boundary conditions object and reinit constraints
    std::shared_ptr<BoundaryConditionsContext> boundary_conditions_context = std::make_shared<BoundaryConditionsContext>(load_step->get_object<NodeBCObject>(),
                                                                                                                         time_step);
    auto boundary_conditions_object = system.get_boundary_conditions_object();
    boundary_conditions_object->set_context(boundary_conditions_context);
    boundary_conditions_object->set_boundary_conditions(prev_load_step->get_object<NodeBCObject>());

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
        // Update the old solution of the system

        system.solve();

        // Update solution
        *system.old_local_solution = *system.current_local_solution;

        system.get_linear_solver()->print_converged_reason();
        if (substep.get_write_results() || substep.get_index() == time_step->get_n_substeps() - 1)
        {
            libMesh::out << "Computing derived quantities..." << std::endl;
            // Call the function to compute the magnetic field

            // equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->compute_derived_quantities();
            // init.comm().barrier();
            // equation_systems.get_system<ExplicitSystem>("NodalSolutionSystem").solution->print(std::cout);
            // libMesh::out << "Writing time substep results..." << std::endl;
            // output_filename = mesh_filename_without_extension + "_out_" + std::to_string(timeStepper->get_cumulative_substep());
            // output_filename = output_filename + ".pvtu";
            // final_output_path = output_path / output_filename;
            std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};
            // // Write the results
            // VTKIO vtk_output(mesh);
            // exo_output.write_equation_systems(final_output_path, equation_systems, &system_names);
            libMesh::out << "Writing time substep results..." << std::endl;
            // session->set_output_filepath(time_step.get_time_stepper().get_cumulative_substep());
            // // Write the results
            // VTKIO_Updated vtk_output(equation_systems->get_mesh(), &equation_systems->get_system<ExplicitSystem>("ElementSolutionSystem"));
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
            // ExplicitSystem &derived_system = equation_systems.get_system<ExplicitSystem>("NodalSolutionSystem");
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

void SimulationController::run_batch_mode()
{
    // libMesh::out << "Running in batch mode with mesh file: " << mesh_filename << std::endl;

    // Create a session with the mesh filename

    // equation_systems.parameters.set<std::set<subdomain_id_type> *>("conductive_subdomain") = &conductive_subdomain;
    // Read and prepare the mesh
    read_mesh();

    // Build and initialize equation systems
    physics();

    // Solve the system for each time step
    for (unsigned int n = 0; n < session->get_n_load_steps(); n++)
    {
        // Create a new load step
        const std::string load_step_filename = session->get_load_steps_folderpath() / ("load_step_" + std::to_string(n + 1) + ".bin");
        load_stepper->read_new_load_step(session->get_load_steps_folderpath() / load_step_filename);
        // Solve the system for the current load step
        solve();
    }
}
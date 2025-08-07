#include "simulation_controller.h"
#include "vtkio_updated.h"
#include "element_types.h"
#include "element_builder.h"
#include "system_builder.h"
#include <libmesh/mesh_communication.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/system_subset_by_subdomain.h>

SimulationController::SimulationController(libMesh::Parallel::Communicator &comm, bool batch_mode, const std::string &mesh_filename)
    : apec_m_model_data(std::make_shared<FEMDataRepository>()),
      batch_mode(batch_mode), comm(comm)
{
    // Initialize the session and mesh
    apec_m_model_data->set_object<Session>(std::make_shared<Session>());
    apec_m_model_data->set_object<ApecMMesh>(std::make_shared<ApecMMesh>(comm));
    auto mesh = apec_m_model_data->get_object<ApecMMesh>();
    equation_systems = std::make_shared<ApecMEquationSystems>(*mesh);
    load_stepper = std::make_shared<LoadStepper>(equation_systems.get());
    if (!mesh_filename.empty())
    {
        this->apec_m_model_data->get_object<Session>()->set_input_mesh_filepath(mesh_filename);
    }
}

void SimulationController::set_transient_effects(ApecMLinearSystem &system)
{
    auto time_step = load_stepper->get_active_load_step()->get_object<TimeStep>();
    auto element_type_factory = apec_m_model_data->get_object<ElementTypeFactory>();
    // Later, this will be a method of system, which will be a viewed object, I still have to think about it
    libMesh::out << "Setting transient effects..." << std::endl;
    if (!time_step->get_transient_effect())
    {
        if (system.has_variable("VOLT"))
        {
            system.release_subset_solve();
            std::set<subdomain_id_type> subdomains;
            for (const auto &pair : *element_type_factory->get_data_map())
            {
                subdomains.insert(pair.first); // pair.first is the subdomain_id_type (element ID)
            }
            SystemSubsetBySubdomain::SubdomainSelectionByList selection(subdomains);
            std::vector<unsigned int> all_var_nums;
            system.get_all_variable_numbers(all_var_nums);
            // Convert std::vector to set
            std::set<unsigned int> all_var_nums_set(all_var_nums.begin(), all_var_nums.end());

            all_var_nums_set.erase(system.variable_number("VOLT"));

            SystemSubsetBySubdomain *subset_ptr = nullptr;
            subset_ptr = new SystemSubsetBySubdomain(system, selection, &all_var_nums_set);
            system.restrict_solve_to(subset_ptr, SUBSET_ZERO);
        }
    }
    else
    {
        // If we have a transient effect, we do not restrict the solve to a subset
        system.release_subset_solve();
        system.restrict_solve_to(nullptr, SUBSET_ZERO);
    }
}

void SimulationController::reinitialize_constrained_dofs(ApecMLinearSystem &system)
{
    auto load_step = load_stepper->get_active_load_step();

    auto prev_load_step = load_stepper->get_old_load_step();

    libMesh::out << "Reinitializing constrained DOFs..." << std::endl;
    // 3. Prepare context view for boundary conditions object and reinit constraints
    auto bc_view = std::make_shared<
        TypedContextView<NodeBCObject, TimeStep>>(
        TypedContextView<NodeBCObject, TimeStep>::from(
            load_step));
    auto boundary_conditions_object = system.get_boundary_conditions_object();
    boundary_conditions_object->set_context_view(bc_view);
    if (prev_load_step != nullptr)
    {
        // Check if we need to update the boundary conditions
        boundary_conditions_object->set_boundary_conditions(prev_load_step->get_object<NodeBCObject>());
    }
    else
    {
        // If there is no previous load step, we assume we are at the beginning
        boundary_conditions_object->set_boundary_conditions(nullptr);
    }
}

void SimulationController::set_load_step_data()
{
    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto &nodal_solution_system = equation_systems->get_system<ApecMSolutionSystem>("NodalSolutionSystem");
    auto &element_solution_system = equation_systems->get_system<ApecMSolutionSystem>("ElementSolutionSystem");
    auto load_step = load_stepper->get_active_load_step();

    libMesh::out << "Setting load step data..." << std::endl;

    auto element_builder_view = std::make_shared<
        TypedContextView<Session, ElementTypeFactory, TimeStep, MaterialsObject, RealConstantsObject, LoadsObject>>(
        TypedContextView<Session, ElementTypeFactory, TimeStep, MaterialsObject, RealConstantsObject, LoadsObject>::from(
            apec_m_model_data,
            load_step));
    auto element_builder = equation_systems->get_element_builder();
    element_builder->set_context_view(element_builder_view);

    // 2. Prepare context view for assembler
    auto assembler_view = std::make_shared<
        TypedContextView<TimeStep, Session>>(
        TypedContextView<TimeStep, Session>::from(
            apec_m_model_data,
            load_step));
    auto assembler = system.get_assembler();

    assembler->set_context_view(assembler_view);
    auto solution_assembler_view = std::make_shared<
        TypedContextView<ElementTypeFactory, TimeStep, Session>>(
        TypedContextView<ElementTypeFactory, TimeStep, Session>::from(
            apec_m_model_data,
            load_step));
    auto assembler_solution = nodal_solution_system.get_assembler();
    assembler_solution->set_context_view(solution_assembler_view);
    auto assembler_element_solution = element_solution_system.get_assembler();
    assembler_element_solution->set_context_view(solution_assembler_view);
}

void SimulationController::set_solver_settings()
{
    auto load_step = load_stepper->get_active_load_step();

    libMesh::out << "Setting solver..." << std::endl;

    auto solver_conf = load_step->get_object<PetscSolverConfiguration>();
    solver_conf->configure_solver();
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
            apec_m_model_data->get_object<Session>()->set_input_mesh_filepath(filepath);
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
    const std::string &mesh_filepath = apec_m_model_data->get_object<Session>()->get_input_mesh_filepath();
    // mesh.clear();
    auto mesh = apec_m_model_data->get_object<ApecMMesh>();

    auto mesh_base = static_cast<MeshBase *>(mesh.get());

    VTKIO_Updated vtkio(*mesh, mesh_filepath);

    apec_m_model_data->set_object<ElementTypeFactory>(std::make_shared<ElementTypeFactory>(vtkio.get_element_types_map_csv()));

    MeshCommunication().broadcast(*mesh);

    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh->prepare_for_use(true);
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
            apec_m_model_data));
    for (int i = 0; i < NUM_SYSTEMS; ++i)
    {
        SystemType type = static_cast<SystemType>(i);
        const std::string sys_name = get_system_name(type);

        if (type == MAIN)
            systems[i] = &equation_systems->add_system<ApecMLinearSystem>(sys_name);
        else
            systems[i] = &equation_systems->add_system<ApecMSolutionSystem>(sys_name);

        SystemBuilder *builder = new SystemBuilder(systems[i], system_builder_view);
        builder->init_variables();
        systems[i]->attach_init_object(*builder);
    }

    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    equation_systems->init();
    equation_systems->print_info();
    equation_systems->get_element_builder()->init();
}

void SimulationController::solve()
{
    auto load_step = load_stepper->get_active_load_step();

    auto time_step = load_step->get_object<TimeStep>();

    auto session = apec_m_model_data->get_object<Session>();

    time_step->print_solve_msg();

    auto element_builder = equation_systems->get_element_builder();
    // read material properties

    // Read loads
    // Here you have to prepare the various contexts for the systems
    // Get the system
    ApecMLinearSystem &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    this->set_solver_settings();

    this->set_load_step_data();

    this->set_transient_effects(system);

    this->reinitialize_constrained_dofs(system);
    // system.get_dof_map().print_dof_constraints();
    // Iterate through the substeps of the current time step
    system.get_system_matrix().zero();

    for (auto &substep : *time_step)
    {
        equation_systems->get_element_builder()->set_elem_builder_misc_data();

        if (substep.get_index() == 0)
        {
            auto prev_load_step = load_stepper->get_old_load_step();

            if (prev_load_step == nullptr)
            {
                // If there is a previous load step, we update the old solution
                element_builder->build_write_element_matrices_rhs(true);
            }
            else
            {
                // If there is no previous load step, we initialize the old solution to zero
                auto prev_element_builder_view = std::make_shared<
                    TypedContextView<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>>(
                    TypedContextView<Session, ElementTypeFactory, MaterialsObject, RealConstantsObject, LoadsObject>::from(
                        apec_m_model_data,
                        prev_load_step));
                // element_builder->set_context_view(prev_element_builder_view);
                element_builder->build_write_element_matrices_rhs(element_builder->check_compute_matrix(prev_element_builder_view));
            }
        }
        substep.print_solve_msg();
        // update time in element data object
        // Assemble & solve the linear system
        libMesh::out << "Solve begins..." << std::endl;
        // Update the old solution of the system

        system.solve();

        system.compute_solution_time_derivative(substep);

        system.rhs->zero();

        // system.get_linear_solver()->print_converged_reason();
        if (substep.get_write_results() || substep.get_index() == time_step->get_n_substeps() - 1)
        {
            libMesh::out << "Computing derived quantities..." << std::endl;
            // Call the function to compute the magnetic field

            equation_systems->get_element_builder()->build_write_solution_vectors("NodalSolutionSystem");
            equation_systems->get_element_builder()->build_write_solution_vectors("ElementSolutionSystem");
            auto nodal_solution_system = &equation_systems->get_system<ApecMSolutionSystem>("NodalSolutionSystem");
            auto element_solution_system = &equation_systems->get_system<ApecMSolutionSystem>("ElementSolutionSystem");
            nodal_solution_system->get_assembler()->compute_derived_quantities();
            element_solution_system->get_assembler()->compute_derived_quantities();
            // equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->compute_derived_quantities();
            std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};
            // // Write the results
            // VTKIO vtk_output(mesh);
            // exo_output.write_equation_systems(final_output_path, equation_systems, &system_names);
            libMesh::out << "Writing time substep results..." << std::endl;
            // session->set_output_filepath(time_step.get_time_stepper().get_cumulative_substep());
            // // Write the results
            session->set_output_filepath(load_stepper->get_cumulative_substep_index());
            VTKIO_Updated vtk_output(equation_systems->get_mesh(), &equation_systems->get_system<ApecMSolutionSystem>("ElementSolutionSystem"));
            vtk_output.write_equation_systems(session->get_path("final_output_path"), *equation_systems, &system_names);
            load_stepper->update_cumulative_substep_index();
        }
        // Update solution
        *system.old_local_solution = *system.current_local_solution;
        // auto &old_time_derivative = system.get_vector("solution_derivative_old");
        // auto &current_time_derivative = system.get_vector("solution_derivative");
        // old_time_derivative.swap(current_time_derivative);

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
    auto session = apec_m_model_data->get_object<Session>();
    for (unsigned int n = 0; n < session->get_n_load_steps(); n++)
    {
        // Create a new load step
        const std::string load_step_filename = session->get_path("load_steps_folderpath") / ("load_step_" + std::to_string(n + 1) + ".bin");
        load_stepper->read_new_load_step(load_step_filename);
        // Solve the system for the current load step
        solve();
    }
}

void SimulationController::clear()
{
    libMesh::out << "Dummy clear() called." << std::endl;
}

void SimulationController::finish()
{
    libMesh::out << "Dummy finish() called." << std::endl;
}

void SimulationController::exit()
{
    libMesh::out << "Dummy exit() called." << std::endl;
}
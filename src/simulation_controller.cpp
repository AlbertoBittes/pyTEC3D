#include "simulation_controller.h"
#include "apec_m_mesh.h"
#include "globals.h"
#include "element_types.h"
#include "element_builder.h"
#include "load_stepper.h"
#include "element_types_factory.h"
#include "solver_configuration.h"
#include "global_constraints_object.h"
#include "apec_m_eq_system_object.h"
#include "apec_m_periodic_boundary_factory.h"
#include "vtkio_updated.h"

#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>

#include <filesystem>

SimulationController::SimulationController(libMesh::Parallel::Communicator &comm, bool batch_mode, const std::string &mesh_filename)
    : apec_m_model_data(std::make_shared<FEMDataRepository>()),
      batch_mode(batch_mode), comm(comm)
{
    // Initialize the session and mesh
    apec_m_model_data->set_object<Session>(std::make_shared<Session>());
    if (!mesh_filename.empty())
    {
        this->apec_m_model_data->get_object<Session>()->set_input_mesh_filepath(mesh_filename);
    }
}

void SimulationController::set_load_step_data()
{
    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto load_step = load_stepper->context<LoadStep>();

    libMesh::out << "Setting load step data..." << std::endl;

    auto element_builder_view = std::make_shared<
        TypedContextView<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>>(
        TypedContextView<Session, ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>::from(
            apec_m_model_data,
            load_step));
    auto element_builder = equation_systems->get_element_builder();
    element_builder->update_context_views(element_builder_view);

    auto bc_view = std::make_shared<
        TypedContextView<TimeStep>>(
        TypedContextView<TimeStep>::from(
            load_step));
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->update_context_views(bc_view);

    // 2. Prepare context view for assembler
    auto assembler_view = std::make_shared<
        TypedContextView<TimeStep, Session>>(
        TypedContextView<TimeStep, Session>::from(
            apec_m_model_data,
            load_step));
    auto assembler = system.get_assembler();
    assembler->set_context_view(assembler_view);
}

void SimulationController::set_solver_settings()
{
    auto load_step = load_stepper->context<LoadStep>();

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
            if (!load_stepper)
            {
                libMesh::out << "Load stepper not initialized. Please read a mesh first.\n";
                break;
            }
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
    auto session = apec_m_model_data->get_object<Session>();
    const std::string &mesh_filepath = session->get_input_mesh_filepath();

    auto mesh = ApecMMesh::read(comm, *apec_m_model_data, mesh_filepath);
    auto elem_types_factory = apec_m_model_data->get_object<ElementTypesFactory>();
    mesh->generate_infinite_elements(elem_types_factory);
    mesh->apply_boundary_info();

    equation_systems = ApecMEquationSystems::build(*mesh, apec_m_model_data);
    load_stepper = std::make_shared<LoadStepper>(equation_systems.get());

    auto periodic_factory = apec_m_model_data->get_object<ApecMPeriodicBoundaryFactory>();
    periodic_factory->initialize(equation_systems.get());

    mesh->prepare_for_use(true);

    libMesh::out << "Mesh Information: " << std::endl;
    mesh->print_info();
}

void SimulationController::physics()
{
    if (!equation_systems)
    {
        libMesh::out << "Equation systems not initialized. Please read a mesh first.\n";
        return;
    }
    equation_systems->initialize_physics();
}

void SimulationController::solve()
{
    if (!equation_systems || !load_stepper)
    {
        libMesh::out << "Solver prerequisites missing. Please read a mesh first.\n";
        return;
    }

    auto load_step = load_stepper->context<LoadStep>();

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

    // system.get_dof_map().print_dof_constraints();
    // Iterate through the substeps of the current time step
    system.get_system_matrix().zero();

    for (auto &substep : *time_step)
    {
        equation_systems->apply_ramped_effects(substep);

        system.reinit_constraints();

        if (substep.get_index() == 0)
        {
            system.delete_rows_cols_constrained_dofs();
            system.check_constrained_dofs();
        }

        element_builder->check_compute_matrix();

        element_builder->build_write_element_matrices_rhs();

        substep.print_solve_msg();
        // update time in element data object
        // Assemble & solve the linear system
        libMesh::out << "Solve begins..." << std::endl;
        // Update the old solution of the system

        system.solve();

        system.compute_solution_time_derivative(substep);

        system.get_dof_map().enforce_constraints_exactly(system);

        system.rhs->zero();

        // system.get_linear_solver()->print_converged_reason();
        if (substep.get_write_results() || substep.get_index() == time_step->get_n_substeps() - 1)
        {
            libMesh::out << "Computing derived quantities..." << std::endl;
            // Call the function to compute the magnetic field

            equation_systems->get_element_builder()->build_write_solution_vectors();

            // equation_systems.parameters.get<std::shared_ptr<AssemblerLinear>>("assembler")->compute_derived_quantities();
            std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};
            // // Write the results
            // VTKIO vtk_output(mesh);
            // exo_output.write_equation_systems(final_output_path, equation_systems, &system_names);
            libMesh::out << "Writing time substep results..." << std::endl;
            // session->set_output_filepath(time_step.get_time_stepper().get_cumulative_substep());
            // // Write the results
            session->set_output_filepath(load_stepper->get_cumulative_substep_index());
            VTKIO_Updated vtk_output(equation_systems->get_mesh(), &equation_systems->get_system<ApecMExplicitSystem>("ElementSolutionSystem"));
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
    unsigned int step_index = 1;
    while (true)
    {
        const auto load_step_path = session->get_path("load_steps_folderpath") /
                                    ("load_step_" + std::to_string(step_index) + ".bin");
        if (!std::filesystem::exists(load_step_path))
            break;

        load_stepper->read_new_load_step(load_step_path);
        solve();
        ++step_index;
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

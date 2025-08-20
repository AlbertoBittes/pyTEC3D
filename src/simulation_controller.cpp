#include "simulation_controller.h"
#include "vtkio_updated.h"
#include "globals.h"
#include "element_types.h"
#include "element_builder.h"
#include "system_builder.h"
#include "load_stepper.h"
#include "inf_elem_builder_updated.h"
#include "element_types_factory.h"
#include "solver_configuration.h"
#include "global_constraints_object.h"
#include "apec_m_eq_system_object.h"

#include <libmesh/boundary_info.h>
#include <libmesh/inf_elem_builder.h>
#include <libmesh/string_to_enum.h>
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

void SimulationController::set_transient_effects()
{
    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto time_step = load_stepper->context<LoadStep>()->get_object<TimeStep>();
    auto element_type_factory = apec_m_model_data->get_object<ElementTypesFactory>();
    // Later, this will be a method of system, which will be a viewed object, I still have to think about it
    libMesh::out << "Setting transient effects..." << std::endl;
    if (!time_step->get_transient_effect())
    {
        if (system.has_variable("VOLT"))
        {
            system.release_subset_solve();
            std::set<subdomain_id_type> subdomains;
            for (const auto &pair : *element_type_factory->get_element_types_map())
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

void SimulationController::check_constrained_dofs()
{
    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto global_constraints_object = system.get_global_constraints_object();
    libMesh::out << "Checking constrained DOFs..." << std::endl;
    // 3. Prepare context view for boundary conditions object and reinit constraints
    global_constraints_object->check_update_boundary_conditions();
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

void SimulationController::set_boundary_info(const std::set<std::string> &bc_arrays)
{
    auto mesh = apec_m_model_data->get_object<ApecMMesh>();

    BoundaryInfo &boundary_info = mesh->get_boundary_info();

    boundary_info.clear();

    for (auto &array : bc_arrays)
    {
        unsigned int node_array_index = mesh->get_node_integer_index(array);
        Globals::BoundaryConditionType bc_enum = Globals::get_key_by_predicate_one(
            Globals::BOUNDARY_CONDITION_INFO_MAP,
            [array](const Globals::BoundaryConditionInfo &info)
            {
                return info.label == array;
            });

        // for (const Node *node_ptr : mesh->local_node_ptr_range())
        for (const Node *node_ptr : mesh->node_ptr_range())
        {
            boundary_id_type boundary_id_mesh = node_ptr->get_extra_integer(node_array_index);
            if (boundary_id_mesh != 0)
            {
                switch (bc_enum)
                {
                case Globals::BoundaryConditionType::INFINITE_BOUNDARY:
                    boundary_info.add_node(node_ptr, boundary_id_mesh);
                    break;
                case Globals::BoundaryConditionType::FLUX_NORMAL:
                    boundary_info.add_node(node_ptr, Globals::FLUX_NORMAL_ID_BC);
                    break;
                case Globals::BoundaryConditionType::FLUX_PARALLEL:
                    boundary_info.add_node(node_ptr, Globals::FLUX_PARALLEL_ID_BC);
                    break;
                default:
                    // Handle unknown
                    break;
                }
            }
        }
        boundary_info.build_side_list_from_node_list();
        // boundary_info.parallel_sync_node_ids();
        // boundary_info.parallel_sync_side_ids();
        // boundary_info.clear_boundary_node_ids();
    }

    boundary_info.print_info();
    boundary_info.print_summary();

    return;
}

void SimulationController::generate_infinite_elements()
{
    auto mesh = apec_m_model_data->get_object<ApecMMesh>();
    auto &boundary_info = mesh->get_boundary_info();
    // In the future change this by looping over the elementtypes
    if (mesh->has_node_integer("infinite_boundary"))
    {

        const std::set<libMesh::boundary_id_type> &side_boundary_ids = boundary_info.get_side_boundary_ids();
        unsigned int etype_num_index = mesh->get_elem_integer_index("ansys_etype");
        unsigned int material_num_index = mesh->get_elem_integer_index("ansys_material_type");
        unsigned int real_constants_num_index = mesh->get_elem_integer_index("ansys_real_constant");
        for (const auto &id : side_boundary_ids)
        {
            // The infinite elements are attached to all elements that build the outer surface of the FEM-region.
            InfElemBuilderUpdated builder(*mesh);
            builder.set_infinite_elem_side_boundary_id(id);
            const Point origin(0.0, 0.0, 0.0);
            builder.build_inf_elem(origin);
            // find the neighbours; for correct linking the two areas
            // mesh->find_neighbors();
            // Reassign subdomain_id() of all infinite elements.
            // and their neighbours. This makes finding the surface
            // between these elemests much easier.
            for (auto &elem : mesh->element_ptr_range())
                if (elem->infinite())
                {
                    elem->subdomain_id() = id;
                    elem->set_extra_integer(etype_num_index, id);
                    elem->set_extra_integer(material_num_index, 0);
                    elem->set_extra_integer(real_constants_num_index, 0);
                }
        }
    }

    return;
}

void SimulationController::read_mesh()
{

    libMesh::out << "Reading mesh..." << std::endl;
    const std::string &mesh_filepath = apec_m_model_data->get_object<Session>()->get_input_mesh_filepath();
    // mesh.clear();
    auto mesh = apec_m_model_data->get_object<ApecMMesh>();

    VTKIO_Updated vtkio(*mesh, mesh_filepath);

    auto element_types_factory = vtkio.get_element_types_factory();

    element_types_factory->populate_element_types_map_from_dataframe();

    apec_m_model_data->set_object<ElementTypesFactory>(element_types_factory);

    MeshCommunication().broadcast(*mesh);

    // orient before generating infinite elems
    MeshTools::Modification::orient_elements(*mesh);

    this->set_boundary_info({"infinite_boundary"});

    this->generate_infinite_elements();

    this->set_boundary_info({"flux_parallel", "flux_normal"});
    // Partition after generating infinite elems
    libMesh::out << "Finding elements' neighbors and partitioning mesh..." << std::endl;
    mesh->prepare_for_use(true);

    // Print information about the mesh to the screen.
    std::cout << "Mesh Information: " << std::endl;
    mesh->print_info();

    return;
}

void SimulationController::physics()
{
    libMesh::out << "Building system of equations..." << std::endl;

    System *systems[Globals::NUM_SYSTEMS] = {nullptr};

    libMesh::out << "Checking system options..." << std::endl;

    auto system_builder_view = std::make_shared<
        TypedContextView<ElementTypesFactory>>(
        TypedContextView<ElementTypesFactory>::from(
            apec_m_model_data));
    for (int i = 0; i < Globals::NUM_SYSTEMS; ++i)
    {
        Globals::SystemType type = static_cast<Globals::SystemType>(i);
        const std::string sys_name = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, type).name;

        if (type == Globals::DOF_SOLUTION_SYSTEM)
            systems[i] = &equation_systems->add_system<ApecMLinearSystem>(sys_name);
        else
            systems[i] = &equation_systems->add_system<ApecMExplicitSystem>(sys_name);

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

    this->set_transient_effects();

    this->check_constrained_dofs();
    // system.get_dof_map().print_dof_constraints();
    // Iterate through the substeps of the current time step
    system.get_system_matrix().zero();

    for (auto &substep : *time_step)
    {
        equation_systems->apply_ramped_effects(substep);

        system.reinit_constraints();

        element_builder->check_compute_matrix();

        element_builder->build_write_element_matrices_rhs();
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
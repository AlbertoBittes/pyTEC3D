#include "apec_m.h"

#include "apec_m_periodic_boundary_factory.h"
#include "element_types_factory.h"
#include "binary_reader.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "solver_configuration.h"
#include "time_step.h"

#include <filesystem>
#include <exception>
#include <utility>
#include <system_error>
#include <typeindex>
#include <unordered_map>

APEC_M::APEC_M(libMesh::Parallel::Communicator &comm)
    : comm(comm),
      model_data_repository(std::make_shared<FEMDataRepository>())
{
    model_data_repository->set_object(std::make_shared<Session>());
}

FEMDataRepository &APEC_M::ensure_mesh_repository(const std::string &name)
{
    auto it = mesh_repository_lookup.find(name);
    if (it != mesh_repository_lookup.end())
        return *mesh_repositories[it->second].second;

    auto repository = std::make_shared<FEMDataRepository>();
    mesh_repositories.emplace_back(name, repository);
    mesh_repository_lookup[name] = mesh_repositories.size() - 1;
    return *repository;
}

std::shared_ptr<FEMDataRepository> APEC_M::find_mesh_repository(const std::string &name) const
{
    auto it = mesh_repository_lookup.find(name);
    if (it == mesh_repository_lookup.end())
        return nullptr;

    return mesh_repositories[it->second].second;
}

std::shared_ptr<ApecMMesh> APEC_M::get_mesh(const std::string &name) const
{
    auto repository = find_mesh_repository(name);
    if (!repository)
        return nullptr;
    return repository->get_object<ApecMMesh>();
}

std::shared_ptr<ApecMMesh> APEC_M::get_active_mesh() const
{
    auto session = model_data_repository->get_object<Session>();
    if (!session || !session->has_active_mesh())
        return nullptr;

    const Session &session_ref = *session;
    return get_mesh(session_ref.get_active_mesh());
}

std::shared_ptr<ApecMEquationSystems> APEC_M::get_equation_system(const std::string &name) const
{
    auto repository = find_mesh_repository(name);
    if (!repository)
        return nullptr;
    return repository->get_object<ApecMEquationSystems>();
}

std::shared_ptr<ApecMEquationSystems> APEC_M::get_active_equation_system() const
{
    auto session = model_data_repository->get_object<Session>();
    if (!session || !session->has_active_mesh())
        return nullptr;

    return get_equation_system(session->get_active_mesh());
}

std::shared_ptr<AnalysisBase> APEC_M::get_analysis(const std::string &name) const
{
    auto repository = find_mesh_repository(name);
    if (!repository)
        return nullptr;
    return repository->get_object<AnalysisBase>();
}

std::shared_ptr<AnalysisBase> APEC_M::get_active_analysis() const
{
    auto session = model_data_repository->get_object<Session>();
    if (!session || !session->has_active_mesh())
        return nullptr;

    return get_analysis(session->get_active_mesh());
}

bool APEC_M::has_mesh(const std::string &name) const
{
    auto repository = find_mesh_repository(name);
    return repository && static_cast<bool>(repository->get_object<ApecMMesh>());
}

bool APEC_M::set_active_mesh(const std::string &name)
{
    auto repository = find_mesh_repository(name);
    if (!repository || !repository->get_object<ApecMMesh>())
        return false;

    auto session = model_data_repository->get_object<Session>();
    if (!session)
        return false;

    session->set_active_mesh(name);

    return true;
}

bool APEC_M::load_model_data(const std::string &filepath)
{
    if (filepath.empty())
        return false;

    auto session = model_data_repository->get_object<Session>();
    if (!session)
        return false;

    auto &session_ref = *session;
    std::filesystem::path input_path(filepath);
    std::filesystem::path base_dir = session_ref.get_working_directory();
    std::error_code ec;
    if (!base_dir.is_absolute())
    {
        auto absolute_base = std::filesystem::absolute(base_dir, ec);
        if (!ec)
            base_dir = absolute_base;
    }

    std::filesystem::path candidate_path = input_path.is_absolute()
                                               ? input_path
                                               : base_dir / input_path;

    candidate_path = candidate_path.lexically_normal();

    if (!std::filesystem::exists(candidate_path, ec) || ec)
    {
        auto message_path = candidate_path;
        if (!message_path.is_absolute())
        {
            auto absolute_message = std::filesystem::absolute(message_path, ec);
            if (!ec)
                message_path = absolute_message;
        }

        libMesh::out << "Model data file must exist in the working directory: " << message_path.lexically_normal() << std::endl;
        return false;
    }

    std::filesystem::path model_path;
    auto absolute_path = std::filesystem::absolute(candidate_path, ec);
    if (!ec)
        model_path = absolute_path.lexically_normal();
    else
        model_path = candidate_path;

    try
    {
        BinaryReader reader;
        reader.read_binary_fem_data(model_path, model_data_repository);
    }
    catch (const std::exception &ex)
    {
        libMesh::out << "Failed to load model data '" << model_path << "': " << ex.what() << std::endl;
        return false;
    }
    catch (...)
    {
        libMesh::out << "Failed to load model data '" << model_path << "' due to an unknown error." << std::endl;
        return false;
    }

    session_ref.set_model_data_filepath(model_path);
    session_ref.check_load_steps_folder();

    return true;
}

void APEC_M::print_model_data() const
{
    if (!model_data_repository)
    {
        libMesh::out << "Model data repository is not initialized." << std::endl;
        return;
    }

    const auto &objects = model_data_repository->get_all_objects();
    if (objects.empty())
    {
        libMesh::out << "Model data repository is empty." << std::endl;
        return;
    }

    static const std::unordered_map<std::type_index, std::string> friendly_names = {
        {typeid(Session), "Session"},
        {typeid(MaterialsObject), "Materials"},
        {typeid(RealConstantsObject), "Real Constants"},
        {typeid(ElementTypesFactory), "Element Types"},
        {typeid(TimeStep), "Time Parameters"},
        {typeid(ApecMPeriodicBoundaryFactory), "Periodic Boundaries"},
        {typeid(PetscSolverConfiguration), "Solver Settings"},
        {typeid(ApecMSettingsObject), "APEC-M Settings"},
        {typeid(DofInterpolationSettingsObject), "DOF Interpolation Settings"}};

    libMesh::out << "Model data contents:" << std::endl;
    for (const auto &entry : objects)
    {
        const auto &type_idx = entry.first;
        const auto &object = entry.second;
        const auto it = friendly_names.find(type_idx);
        const std::string &name = (it != friendly_names.end()) ? it->second : type_idx.name();

        libMesh::out << "- " << name << std::endl;
        if (object)
            object->print_data_info();
        else
            libMesh::out << "  (object not initialized)" << std::endl;
    }
}

bool APEC_M::distribute_boundary_ids()
{
    auto session = model_data_repository->get_object<Session>();
    if (!session)
        return false;

    Session &session_ref = *session;
    auto mesh = get_active_mesh();
    if (!mesh)
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

    mesh->apply_boundary_info();
    libMesh::out << "Boundary IDs distributed for mesh '" << session_ref.get_active_mesh() << "'." << std::endl;
    return true;
}

bool APEC_M::generate_infinite_elements()
{
    auto session = model_data_repository->get_object<Session>();
    if (!session)
        return false;

    Session &session_ref = *session;
    auto mesh = get_active_mesh();
    if (!mesh)
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

    mesh->apply_boundary_info();

    auto element_factory = model_data_repository->get_object<ElementTypesFactory>();
    if (!element_factory)
    {
        libMesh::out << "No ElementTypesFactory available. Unable to generate infinite elements." << std::endl;
        return false;
    }

    ElementTypesMap *map = element_factory->get_element_types_map();
    if (!map || map->empty())
    {
        libMesh::out << "ElementTypesFactory has no element types. Skipping infinite element generation." << std::endl;
        return false;
    }

    auto &boundary_info = mesh->get_boundary_info();
    const auto &side_boundary_ids = boundary_info.get_side_boundary_ids();

    bool has_valid_boundary = false;
    for (const auto &pair : *map)
    {
        const auto boundary_id = pair.first;
        if (side_boundary_ids.find(boundary_id) == side_boundary_ids.end())
        {
            libMesh::out << "Warning: Boundary ID " << boundary_id
                         << " not found on mesh. Skipping infinite elements for this ID." << std::endl;
            continue;
        }
        has_valid_boundary = true;
    }

    if (!has_valid_boundary)
    {
        libMesh::out << "No matching boundary IDs found. Infinite elements were not generated." << std::endl;
        return false;
    }

    mesh->generate_infinite_elements(element_factory);
    // mesh->apply_boundary_info();
    libMesh::out << "Generated infinite elements for mesh '" << session_ref.get_active_mesh() << "'." << std::endl;
    return true;
}

bool APEC_M::generate_nodal_couplings_for_mesh(const std::string &name,
                                               const std::shared_ptr<ApecMMesh> &mesh)
{
    if (!mesh)
    {
        libMesh::out << "No mesh provided for nodal coupling generation." << std::endl;
        return false;
    }

    libMesh::out << "Generating nodal couplings (stub) for mesh '" << name << "'." << std::endl;
    return true;
}

bool APEC_M::generate_nodal_couplings()
{
    auto session = model_data_repository->get_object<Session>();
    if (!session || !session->has_active_mesh())
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

    auto &session_ref = *session;
    auto mesh = get_active_mesh();
    return generate_nodal_couplings_for_mesh(session_ref.get_active_mesh(), mesh);
}

bool APEC_M::generate_periodic_boundaries()
{
    auto equation_system = get_active_equation_system();

    auto session = model_data_repository->get_object<Session>();
    if (!session || !session->has_active_mesh())
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

    auto &session_ref = *session;
    if (!equation_system)
    {
        libMesh::out << "No equation system available for active mesh." << std::endl;
        return false;
    }

    auto periodic_factory = model_data_repository->get_object<ApecMPeriodicBoundaryFactory>();
    if (!periodic_factory)
    {
        libMesh::out << "No periodic boundary factory available." << std::endl;
        return false;
    }

    periodic_factory->initialize(equation_system.get());
    libMesh::out << "Periodic boundaries generated for mesh '" << session_ref.get_active_mesh() << "'." << std::endl;
    return true;
}

std::size_t APEC_M::load_step_file_count() const
{
    auto analysis = std::dynamic_pointer_cast<TransientAnalysis>(get_active_analysis());
    if (!analysis)
        return 0;
    return analysis->load_step_file_count();
}

std::vector<std::filesystem::path> APEC_M::list_load_step_files() const
{
    auto analysis = std::dynamic_pointer_cast<TransientAnalysis>(get_active_analysis());
    if (!analysis)
        return {};
    return analysis->list_load_step_files();
}

bool APEC_M::solve_load_step(const std::string &filename)
{
    auto analysis = std::dynamic_pointer_cast<TransientAnalysis>(get_active_analysis());
    if (!analysis)
    {
        libMesh::out << "No load step analysis available for the active mesh." << std::endl;
        return false;
    }

    if (filename.empty())
    {
        libMesh::out << "Please provide a load step filename with -solve." << std::endl;
        return false;
    }

    return analysis->solve(filename);
}

bool APEC_M::solve_transient_load_steps(const std::string &base_name)
{
    auto analysis = std::dynamic_pointer_cast<TransientAnalysis>(get_active_analysis());
    if (!analysis)
    {
        libMesh::out << "No load step analysis available for the active mesh." << std::endl;
        return false;
    }

    return analysis->solve_transient(base_name);
}

bool APEC_M::set_mode(Session::Mode mode)
{
    auto session = model_data_repository->get_object<Session>();
    if (!session)
        return false;

    Session &session_ref = *session;
    auto current = session_ref.get_mode();
    if (current == mode)
        return true;

    bool transition_ok = true;

    if (current == Session::Mode::PreProcessing && mode == Session::Mode::Analysis)
    {
        transition_ok = enter_analysis_mode();
    }
    else if (current == Session::Mode::Analysis && mode == Session::Mode::PreProcessing)
    {
        enter_preprocessing_mode();
    }

    if (!transition_ok)
        return false;

    session_ref.set_mode(mode);
    libMesh::out << "Program mode set to " << Session::mode_to_string(mode) << std::endl;
    return true;
}

void APEC_M::read_mesh(const std::string &name, const std::string &filepath)
{
    auto session = model_data_repository->get_object<Session>();
    if (session)
        session->set_input_mesh_filepath(filepath);

    auto &repository = ensure_mesh_repository(name);

    auto mesh = ApecMMesh::read_mesh(comm, repository, filepath);
    if (!mesh)
        return;

    mesh->set_parent_object(this);
    repository.set_object<ApecMMesh>(mesh);

    auto equation_system = ApecMEquationSystems::build(*mesh);
    repository.set_object<ApecMEquationSystems>(equation_system);

    set_active_mesh(name);
}

void APEC_M::clear_mesh(const std::string &name)
{
    auto it = mesh_repository_lookup.find(name);
    if (it == mesh_repository_lookup.end())
        return;

    const auto removed_index = it->second;
    mesh_repository_lookup.erase(it);
    mesh_repositories.erase(mesh_repositories.begin() + removed_index);

    for (auto &entry : mesh_repository_lookup)
        if (entry.second > removed_index)
            --entry.second;

    auto session = model_data_repository->get_object<Session>();
    if (session && session->has_active_mesh() && session->get_active_mesh() == name)
    {
        if (!mesh_repositories.empty())
        {
            const auto &next_mesh = mesh_repositories.front().first;
            if (!set_active_mesh(next_mesh))
                session->set_active_mesh("");
        }
        else
            session->set_active_mesh("");
    }
}

bool APEC_M::enter_analysis_mode()
{
    if (mesh_repositories.empty())
    {
        libMesh::out << "No meshes available. Cannot enter analysis mode." << std::endl;
        return false;
    }

    bool initialized_any = false;

    auto element_factory = model_data_repository->get_object<ElementTypesFactory>();
    auto periodic_factory = model_data_repository->get_object<ApecMPeriodicBoundaryFactory>();

    for (auto &entry : mesh_repositories)
    {
        const auto &name = entry.first;
        auto &repository = entry.second;

        if (!repository)
            continue;

        auto mesh = repository->get_object<ApecMMesh>();
        if (!mesh)
            continue;

        libMesh::out << "Preparing mesh '" << name << "' for analysis." << std::endl;

        if (element_factory)
        {
            mesh->generate_infinite_elements(element_factory);
            libMesh::out << "Generated infinite elements for mesh '" << name << "'." << std::endl;
        }
        else
        {
            libMesh::out << "No ElementTypesFactory available. Skipping infinite elements for mesh '" << name << "'." << std::endl;
        }

        mesh->apply_boundary_info();

        generate_nodal_couplings_for_mesh(name, mesh);

        auto equation_system = repository->get_object<ApecMEquationSystems>();
        if (!equation_system)
        {
            equation_system = ApecMEquationSystems::build(*mesh);
            repository->set_object<ApecMEquationSystems>(equation_system);
        }

        equation_system->initialize_systems(model_data_repository);

        if (periodic_factory)
        {
            periodic_factory->initialize(equation_system.get());
            libMesh::out << "Periodic boundaries generated for mesh '" << name << "'." << std::endl;
        }
        else
        {
            libMesh::out << "No periodic boundary factory available for mesh '" << name << "'." << std::endl;
        }

        mesh->prepare_for_use(true);

        mesh->print_info();

        equation_system->initialize_physics();

        auto analysis = repository->get_object<AnalysisBase>();
        if (!analysis)
        {
            auto legacy_load_step = repository->get_object<TransientAnalysis>();
            if (legacy_load_step)
            {
                analysis = legacy_load_step;
                repository->set_object<AnalysisBase>(analysis);
            }
        }
        if (!analysis)
        {
            analysis = AnalysisBase::build(equation_system.get(), model_data_repository);
            if (!analysis)
            {
                libMesh::out << "Failed to create analysis object for mesh '" << name << "'." << std::endl;
                continue;
            }
            repository->set_object<AnalysisBase>(analysis);
        }
        else
        {
            analysis->set_parent_object(equation_system.get());
        }

        if (!analysis->initialize_model_context(model_data_repository))
        {
            libMesh::out << "Failed to initialize analysis context for mesh '" << name << "'." << std::endl;
            continue;
        }

        equation_system->set_analysis(analysis);

        initialized_any = true;
    }

    if (!initialized_any)
    {
        libMesh::out << "No equation systems initialized. Cannot enter analysis mode." << std::endl;
        return false;
    }

    return true;
}

void APEC_M::enter_preprocessing_mode()
{
    mesh_repositories.clear();
    mesh_repository_lookup.clear();
    auto session = model_data_repository->get_object<Session>();
    if (session)
        session->set_active_mesh("");
}

void APEC_M::set_cli(std::unique_ptr<CommandLineInterface> cli_instance)
{
    cli = std::move(cli_instance);
}

void APEC_M::run_cli()
{
    if (cli)
        cli->run_interactive(*this);
}

void APEC_M::run_batch_file(const std::string &filepath)
{
    if (cli)
        cli->run_batch(*this, filepath);
}

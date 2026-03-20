#include "apec_m.h"

#include "augment_sparsity_on_nodes.h"
#include "element_types_factory.h"
#include "fe_data_object.h"
#include "binary_reader.h"
#include "solver_configuration.h"
#include "time_step.h"
#include "vtk_data_map_object.h"

#include <filesystem>
#include <optional>
#include <exception>
#include <utility>
#include <system_error>
#include <set>
#include <algorithm>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <variant>

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
    if (!session->has_active_mesh())
        return nullptr;

    return get_equation_system(session->get_active_mesh());
}

std::vector<std::shared_ptr<ApecMEquationSystems>> APEC_M::get_all_equation_systems() const
{
    std::vector<std::shared_ptr<ApecMEquationSystems>> systems;
    systems.reserve(mesh_repositories.size());
    for (const auto &entry : mesh_repositories)
    {
        const auto &repository = entry.second;
        if (!repository)
            continue;
        auto equation_system = repository->get_object<ApecMEquationSystems>();
        if (equation_system)
            systems.push_back(equation_system);
    }
    return systems;
}

std::vector<std::pair<std::string, std::shared_ptr<ApecMEquationSystems>>> APEC_M::get_all_equation_systems_with_names() const
{
    std::vector<std::pair<std::string, std::shared_ptr<ApecMEquationSystems>>> systems;
    systems.reserve(mesh_repositories.size());
    for (const auto &entry : mesh_repositories)
    {
        const auto &name = entry.first;
        const auto &repository = entry.second;
        if (!repository)
            continue;
        auto equation_system = repository->get_object<ApecMEquationSystems>();
        if (equation_system)
            systems.emplace_back(name, equation_system);
    }
    return systems;
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

bool APEC_M::load_data(const std::vector<std::string> &filepaths)
{
    if (filepaths.empty())
    {
        libMesh::out << "No data files provided." << std::endl;
        return false;
    }

    auto session = model_data_repository->get_object<Session>();

    auto base_dir = session->get_working_directory();
    std::error_code ec;
    if (!base_dir.is_absolute())
    {
        auto absolute_base = std::filesystem::absolute(base_dir, ec);
        if (!ec)
            base_dir = absolute_base;
    }

    auto analysis = get_analysis();

    bool any_loaded = false;
    for (const auto &raw_path : filepaths)
    {
        std::filesystem::path input_path(raw_path);
        std::filesystem::path resolved;
        if (input_path.is_absolute())
        {
            resolved = input_path;
        }
        else
        {
            std::filesystem::path effective_base = base_dir;
            if (session && input_path.parent_path().empty() && input_path.filename() == "model_data.bin")
            {
                const auto &model_data_folder = session->get_path("model_data_folderpath");
                if (!model_data_folder.empty())
                    effective_base = model_data_folder;
            }
            resolved = effective_base / input_path;
        }
        resolved = resolved.lexically_normal();

        auto equation_systems = get_active_equation_system();

        try
        {
            BinaryReader reader;
            reader.read_binary_fem_data(resolved.string(), model_data_repository, equation_systems);
        }
        catch (const std::exception &ex)
        {
            libMesh::out << "Failed to load data from '" << resolved << "': " << ex.what() << std::endl;
            continue;
        }
        catch (...)
        {
            libMesh::out << "Failed to load data from '" << resolved << "' due to an unknown error." << std::endl;
            continue;
        }

        any_loaded = true;
        libMesh::out << "Loaded data from '" << resolved << "'." << std::endl;
    }

    return any_loaded;
}

bool APEC_M::solve_load_step()
{
    auto analysis = get_analysis();
    auto eq_system = get_active_equation_system();
    if (!analysis || !eq_system)
    {
        libMesh::out << "No analysis available for the active mesh." << std::endl;
        return false;
    }

    analysis->set_parent_object(this);
    return analysis->solve();
}

bool APEC_M::finish_analysis()
{
    auto analysis = get_analysis();
    if (!analysis)
    {
        libMesh::out << "No analysis available for the active mesh." << std::endl;
        return false;
    }

    analysis->set_parent_object(this);
    return analysis->finish();
}

bool APEC_M::set_mode(Session::Mode mode, std::optional<AnalysisBase::Type> analysis_type)
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
        if (!analysis_type)
        {
            libMesh::out << "Analysis type must be specified when entering analysis mode." << std::endl;
            return false;
        }

        transition_ok = enter_analysis_mode(analysis_type);
    }
    else if (current == Session::Mode::Analysis && mode == Session::Mode::PreProcessing)
    {
        enter_preprocessing_mode();
    }
    else if (analysis_type)
    {
        libMesh::out << "Ignoring analysis type when not entering analysis mode." << std::endl;
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

    auto equation_system = ApecMEquationSystems::build(*mesh, this);
    repository.set_object<ApecMEquationSystems>(equation_system);

    set_active_mesh(name);
    if (session)
        session->initialize_folder_paths();
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

bool APEC_M::enter_analysis_mode(std::optional<AnalysisBase::Type> analysis_type)
{
    if (mesh_repositories.empty())
    {
        libMesh::out << "No meshes available. Cannot enter analysis mode." << std::endl;
        return false;
    }

    if (!analysis_type)
    {
        libMesh::out << "Analysis type not provided. Cannot enter analysis mode." << std::endl;
        return false;
    }

    shared_analysis.reset();

    if (!validate_model_for_analysis())
    {
        libMesh::out << "Model consistency check failed. Cannot enter analysis mode." << std::endl;
        return false;
    }

    bool initialized_any = false;

    auto analysis_context = std::make_shared<AnalysisBase::ContextView>(AnalysisBase::ContextView::from(model_data_repository));

    if (!analysis_context)
    {
        libMesh::out << "Failed to create analysis context for analysis mode." << std::endl;
        return false;
    }

    if (!shared_analysis)
    {
        shared_analysis = AnalysisBase::build(this, analysis_context, *analysis_type);
        if (!shared_analysis)
        {
            libMesh::out << "Failed to create analysis object for mesh." << std::endl;
        }
    }

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

        mesh->generate_nodal_couplings();

        auto equation_system = repository->get_object<ApecMEquationSystems>();
        if (!equation_system)
        {
            equation_system = ApecMEquationSystems::build(*mesh, this);
            repository->set_object<ApecMEquationSystems>(equation_system);
        }
        equation_system->initialize_systems(model_data_repository);

        auto &main_system = equation_system->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
        AugmentSparsityOnNodes *gf = new AugmentSparsityOnNodes(*mesh);
        main_system.get_dof_map().add_coupling_functor(*gf);
        libMesh::out << "Added DOF coupling sparsity functor for mesh '" << name << "'." << std::endl;

        mesh->prepare_for_use(true);

        mesh->print_info();

        equation_system->initialize_physics();

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
    shared_analysis.reset();
    mesh_repositories.clear();
    mesh_repository_lookup.clear();
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

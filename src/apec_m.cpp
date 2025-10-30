#include "apec_m.h"

#include "apec_m_periodic_boundary_factory.h"
#include "element_types_factory.h"

#include <filesystem>
#include <utility>

APEC_M::APEC_M(libMesh::Parallel::Communicator &comm)
    : comm(comm),
      model_data_repository(std::make_shared<FEMDataRepository>()),
      session_object(std::make_shared<Session>())
{
}

APEC_M::MeshContext &APEC_M::ensure_context(const std::string &name)
{
    return mesh_registry.try_emplace(name).first->second;
}

Session &APEC_M::session()
{
    return *session_object;
}

const Session &APEC_M::session() const
{
    return *session_object;
}

std::shared_ptr<ApecMMesh> APEC_M::get_mesh(const std::string &name) const
{
    auto it = mesh_registry.find(name);
    if (it == mesh_registry.end())
        return nullptr;
    return it->second.mesh;
}

std::shared_ptr<ApecMEquationSystems> APEC_M::get_equation_system(const std::string &name) const
{
    auto it = mesh_registry.find(name);
    if (it == mesh_registry.end())
        return nullptr;
    return it->second.equation_system;
}

bool APEC_M::has_mesh(const std::string &name) const
{
    auto it = mesh_registry.find(name);
    return it != mesh_registry.end() && static_cast<bool>(it->second.mesh);
}

bool APEC_M::set_active_mesh(const std::string &name)
{
    auto it = mesh_registry.find(name);
    if (it == mesh_registry.end() || !it->second.mesh)
        return false;

    session_object->set_active_mesh(name);

    if (it->second.equation_system)
    {
        const auto &active_analysis = it->second.equation_system->get_active_analysis();
        if (active_analysis)
            session_object->set_active_analysis(*active_analysis);
    }

    return true;
}

std::optional<ApecMEquationSystems::AnalysisActivationState> APEC_M::activate_analysis(const std::string &mesh_name, const std::string &analysis_name)
{
    auto it = mesh_registry.find(mesh_name);
    if (it == mesh_registry.end() || !it->second.mesh || !it->second.equation_system)
        return std::nullopt;

    const auto state = it->second.equation_system->activate_analysis(analysis_name);

    if (!set_active_mesh(mesh_name))
        return std::nullopt;

    return state;
}

bool APEC_M::load_model_data(const std::string &filepath)
{
    if (filepath.empty())
        return false;

    std::filesystem::path model_path = session_object->get_working_directory() / filepath;

    if (!std::filesystem::exists(model_path))
        return false;

    session_object->set_model_data_filepath(model_path);

    // TODO: populate model_data_repository with actual objects (materials, element types, etc.)

    return true;
}

void APEC_M::read_mesh(const std::string &name, const std::string &filepath)
{
    session_object->set_input_mesh_filepath(filepath);

    auto &context = ensure_context(name);

    context.mesh = ApecMMesh::read(comm, *model_data_repository, filepath);

    if (auto element_factory = model_data_repository->get_object<ElementTypesFactory>())
        context.mesh->generate_infinite_elements(element_factory);

    context.mesh->apply_boundary_info();

    context.equation_system = ApecMEquationSystems::build(*context.mesh, model_data_repository);
    context.equation_system->clear_active_analysis();

    if (auto periodic_factory = model_data_repository->get_object<ApecMPeriodicBoundaryFactory>())
        periodic_factory->initialize(context.equation_system.get());

    context.mesh->prepare_for_use(true);

    set_active_mesh(name);
}

void APEC_M::clear_mesh(const std::string &name)
{
    mesh_registry.erase(name);
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

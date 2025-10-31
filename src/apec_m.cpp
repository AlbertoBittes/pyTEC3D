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

std::shared_ptr<ApecMMesh> APEC_M::get_active_mesh() const
{
    if (!session_object->has_active_mesh())
        return nullptr;
    return get_mesh(session_object->get_active_mesh());
}

std::shared_ptr<ApecMEquationSystems> APEC_M::get_equation_system(const std::string &name) const
{
    auto it = mesh_registry.find(name);
    if (it == mesh_registry.end())
        return nullptr;
    return it->second.equation_system;
}

std::shared_ptr<ApecMEquationSystems> APEC_M::get_active_equation_system() const
{
    if (!session_object->has_active_mesh())
        return nullptr;
    return get_equation_system(session_object->get_active_mesh());
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

    return true;
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

bool APEC_M::distribute_boundary_ids()
{
    auto mesh = get_active_mesh();
    if (!mesh)
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

    mesh->apply_boundary_info();
    libMesh::out << "Boundary IDs distributed for mesh '" << session_object->get_active_mesh() << "'." << std::endl;
    return true;
}

bool APEC_M::generate_infinite_elements()
{
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
    mesh->apply_boundary_info();
    libMesh::out << "Generated infinite elements for mesh '" << session_object->get_active_mesh() << "'." << std::endl;
    return true;
}

bool APEC_M::generate_nodal_couplings()
{
    if (!session_object->has_active_mesh())
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

    libMesh::out << "Generating nodal couplings (stub) for mesh '" << session_object->get_active_mesh() << "'." << std::endl;
    return true;
}

bool APEC_M::generate_periodic_boundaries()
{
    auto equation_system = get_active_equation_system();

    if (!session_object->has_active_mesh())
    {
        libMesh::out << "No active mesh selected." << std::endl;
        return false;
    }

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
    libMesh::out << "Periodic boundaries generated for mesh '" << session_object->get_active_mesh() << "'." << std::endl;
    return true;
}

bool APEC_M::set_mode(Session::Mode mode)
{
    auto current = session_object->get_mode();
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

    session_object->set_mode(mode);
    libMesh::out << "Program mode set to " << Session::mode_to_string(mode) << std::endl;
    return true;
}

void APEC_M::read_mesh(const std::string &name, const std::string &filepath)
{
    session_object->set_input_mesh_filepath(filepath);

    auto &context = ensure_context(name);

    context.mesh = ApecMMesh::read(comm, *model_data_repository, filepath);
    if (context.mesh)
        context.mesh->set_parent_object(this);

    if (auto element_factory = model_data_repository->get_object<ElementTypesFactory>())
        context.mesh->generate_infinite_elements(element_factory);

    context.mesh->apply_boundary_info();

    context.equation_system = ApecMEquationSystems::build(*context.mesh, model_data_repository);

    if (auto periodic_factory = model_data_repository->get_object<ApecMPeriodicBoundaryFactory>())
        periodic_factory->initialize(context.equation_system.get());

    context.mesh->prepare_for_use(true);

    set_active_mesh(name);
}

void APEC_M::clear_mesh(const std::string &name)
{
    mesh_registry.erase(name);

    if (session_object->has_active_mesh() && session_object->get_active_mesh() == name)
    {
        if (!mesh_registry.empty())
        {
            auto next_mesh = mesh_registry.begin()->first;
            if (!set_active_mesh(next_mesh))
                session_object->set_active_mesh("");
        }
        else
            session_object->set_active_mesh("");
    }
}

bool APEC_M::enter_analysis_mode()
{
    if (mesh_registry.empty())
    {
        libMesh::out << "No meshes available. Cannot enter analysis mode." << std::endl;
        return false;
    }

    bool initialized_any = false;

    for (auto &entry : mesh_registry)
    {
        auto &name = entry.first;
        auto &context = entry.second;

        if (!context.mesh)
            continue;

        libMesh::out << "Preparing mesh '" << name << "' for analysis." << std::endl;
        context.mesh->prepare_for_use(true);

        if (!context.equation_system)
        {
            context.equation_system = ApecMEquationSystems::build(*context.mesh, model_data_repository);
        }

        context.equation_system->initialize_physics();
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
    mesh_registry.clear();
    session_object->set_active_mesh("");
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

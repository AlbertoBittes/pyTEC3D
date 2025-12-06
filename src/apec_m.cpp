#include "apec_m.h"

#include "apec_m_periodic_boundary_factory.h"
#include "apec_m_settings.h"
#include "element_types_factory.h"
#include "fe_data_object.h"
#include "binary_reader.h"
#include "materials_object.h"
#include "real_constants_object.h"
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

std::shared_ptr<AnalysisBase> APEC_M::get_analysis() const
{
    return shared_analysis;
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
    if (!session)
        return false;

    const bool analysis_mode = session->get_mode() == Session::Mode::Analysis;
    std::unordered_set<BinaryLabelId> analysis_disallowed_labels;
    if (analysis_mode)
    {
        analysis_disallowed_labels = {BinaryLabelId::ELEMENT_TYPES,
                                      BinaryLabelId::MATERIALS,
                                      BinaryLabelId::REAL_CONSTANTS,
                                      BinaryLabelId::SOLVER_SETTINGS,
                                      BinaryLabelId::PERIODIC_BOUNDARIES,
                                      BinaryLabelId::DOF_COUPLINGS};
    }

    auto base_dir = session->get_working_directory();
    std::error_code ec;
    if (!base_dir.is_absolute())
    {
        auto absolute_base = std::filesystem::absolute(base_dir, ec);
        if (!ec)
            base_dir = absolute_base;
    }

    auto analysis = get_analysis();
    if (analysis)
        analysis->set_parent_object(this);

    bool any_loaded = false;
    for (const auto &raw_path : filepaths)
    {
        std::filesystem::path input_path(raw_path);
        std::filesystem::path resolved = input_path.is_absolute() ? input_path : base_dir / input_path;
        resolved = resolved.lexically_normal();

        auto equation_systems = get_active_equation_system();

        try
        {
            BinaryReader reader;
            if (analysis_mode)
                reader.set_disallowed_labels(analysis_disallowed_labels);
            reader.read_binary_fem_data(resolved.string(), model_data_repository, equation_systems);
        }
        catch (const BinaryLabelDisallowedError &ex)
        {
            libMesh::out << "Failed to load data from '" << resolved << "': " << ex.what() << std::endl;
            throw;
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

bool APEC_M::solve_transient_load_steps()
{
    auto analysis = std::dynamic_pointer_cast<TransientAnalysis>(get_analysis());
    auto eq_system = get_active_equation_system();
    if (!analysis || !eq_system)
    {
        libMesh::out << "No load step analysis available for the active mesh." << std::endl;
        return false;
    }

    analysis->set_parent_object(this);
    return analysis->solve_transient();
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

bool APEC_M::validate_model_for_analysis()
{
    bool valid = true;

    auto element_factory = model_data_repository->get_object<ElementTypesFactory>();
    auto materials_object = model_data_repository->get_object<MaterialsObject>();
    auto real_constants_object = model_data_repository->get_object<RealConstantsObject>();
    auto periodic_factory = model_data_repository->get_object<ApecMPeriodicBoundaryFactory>();

    std::unordered_set<int> available_etypes;
    std::unordered_set<int> available_materials;
    std::unordered_set<int> available_real_constants;

    if (!element_factory || !element_factory->get_element_types_map() || element_factory->get_element_types_map()->empty())
    {
        libMesh::out << "No element types loaded. Cannot enter analysis mode." << std::endl;
        valid = false;
    }
    else
    {
        for (const auto &entry : *element_factory->get_element_types_map())
            available_etypes.insert(static_cast<int>(entry.first));
    }

    if (!materials_object || materials_object->get_materials_map().empty())
    {
        libMesh::out << "No materials loaded. Cannot enter analysis mode." << std::endl;
        valid = false;
    }
    else
    {
        for (const auto &entry : materials_object->get_materials_map())
            available_materials.insert(static_cast<int>(entry.first));
    }

    if (!real_constants_object || real_constants_object->get_constants_map().empty())
    {
        libMesh::out << "No real constants loaded. Cannot enter analysis mode." << std::endl;
        valid = false;
    }
    else
    {
        for (const auto &entry : real_constants_object->get_constants_map())
            available_real_constants.insert(static_cast<int>(entry.first));
    }

    if (periodic_factory)
        periodic_factory->populate_periodic_boundaries_from_map();

    for (const auto &entry : mesh_repositories)
    {
        const auto &mesh_name = entry.first;
        auto mesh = entry.second ? entry.second->get_object<ApecMMesh>() : nullptr;
        if (!mesh)
            continue;

        const auto etype_index = mesh->get_elem_integer_index("ansys_etype");
        const auto material_index = mesh->get_elem_integer_index("ansys_material_type");
        const auto real_constant_index = mesh->get_elem_integer_index("ansys_real_constant");

        auto collect_used_ids = [&](unsigned int index) {
            std::unordered_set<int> used;
            if (index == libMesh::invalid_uint)
                return used;
            for (const auto elem : mesh->element_ptr_range())
                used.insert(elem->get_extra_integer(index));
            return used;
        };

        if (etype_index == libMesh::invalid_uint)
        {
            libMesh::out << "Mesh '" << mesh_name << "' is missing element integer array 'ansys_etype'." << std::endl;
            valid = false;
        }
        if (material_index == libMesh::invalid_uint)
        {
            libMesh::out << "Mesh '" << mesh_name << "' is missing element integer array 'ansys_material_type'." << std::endl;
            valid = false;
        }
        if (real_constant_index == libMesh::invalid_uint)
        {
            libMesh::out << "Mesh '" << mesh_name << "' is missing element integer array 'ansys_real_constant'." << std::endl;
            valid = false;
        }

        const auto used_etypes = collect_used_ids(etype_index);
        const auto used_materials = collect_used_ids(material_index);
        const auto used_real_constants = collect_used_ids(real_constant_index);

        auto report_missing = [&](const std::string &label,
                                  const std::unordered_set<int> &used,
                                  const std::unordered_set<int> &available) {
            std::vector<int> missing;
            missing.reserve(used.size());
            for (int value : used)
                if (!available.count(value))
                    missing.push_back(value);
            if (!missing.empty())
            {
                std::sort(missing.begin(), missing.end());
                missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
                libMesh::out << "Mesh '" << mesh_name << "' references undefined " << label << " IDs: ";
                for (std::size_t i = 0; i < missing.size(); ++i)
                {
                    if (i)
                        libMesh::out << ", ";
                    libMesh::out << missing[i];
                }
                libMesh::out << std::endl;
                valid = false;
            }
        };

        if (!available_etypes.empty())
            report_missing("element type", used_etypes, available_etypes);
        if (!available_materials.empty())
            report_missing("material", used_materials, available_materials);
        if (!available_real_constants.empty())
            report_missing("real constant", used_real_constants, available_real_constants);

        if (periodic_factory)
        {
            const auto &periodic_boundaries = periodic_factory->get_periodic_boundaries();
            if (!periodic_boundaries.empty())
            {
                auto &boundary_info = mesh->get_boundary_info();
                if (boundary_info.get_side_boundary_ids().empty())
                    mesh->apply_boundary_info();

                const auto &side_ids = boundary_info.get_side_boundary_ids();
                if (side_ids.empty())
                {
                    libMesh::out << "Mesh '" << mesh_name << "' has no boundary IDs but periodic boundaries are defined." << std::endl;
                    valid = false;
                }

                for (const auto &pb : periodic_boundaries)
                {
                    if (side_ids.find(pb.myboundary) == side_ids.end())
                    {
                        libMesh::out << "Periodic boundary main ID " << pb.myboundary << " not found on mesh '" << mesh_name << "'." << std::endl;
                        valid = false;
                    }
                    if (side_ids.find(pb.pairedboundary) == side_ids.end())
                    {
                        libMesh::out << "Periodic boundary paired ID " << pb.pairedboundary << " not found on mesh '" << mesh_name << "'." << std::endl;
                        valid = false;
                    }
                }
            }
        }
    }

    return valid;
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

    auto element_factory = model_data_repository->get_object<ElementTypesFactory>();
    auto periodic_factory = model_data_repository->get_object<ApecMPeriodicBoundaryFactory>();

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
            equation_system = ApecMEquationSystems::build(*mesh, this);
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

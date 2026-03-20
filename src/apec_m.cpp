#include "apec_m.h"

#include "augment_sparsity_on_nodes.h"
#include "element_types_factory.h"
#include "fe_data_object.h"
#include "binary_reader.h"
#include "solver_configuration.h"
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

void APEC_M::read_mesh(const std::string &name)
{
    auto mesh = ApecMMesh::read_mesh(comm, repository, filepath);

    mesh->set_parent_object(this);
    repository.set_object<ApecMMesh>(mesh);

    auto equation_system = ApecMEquationSystems::build(*mesh, this);
    repository.set_object<ApecMEquationSystems>(equation_system);

    set_active_mesh(name);
    if (session)
        session->initialize_folder_paths();
}

bool APEC_M::enter_analysis_mode()
{
    shared_analysis.reset();

    auto analysis_context = std::make_shared<AnalysisBase::ContextView>(AnalysisBase::ContextView::from(model_data_repository));

    shared_analysis = AnalysisBase::build(this, analysis_context);

    libMesh::out << "Preparing mesh '" << name << "' for analysis." << std::endl;

    mesh->generate_nodal_couplings();

    equation_system->initialize_systems(model_data_repository);

    auto &main_system = equation_system->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    AugmentSparsityOnNodes *gf = new AugmentSparsityOnNodes(*mesh);
    main_system.get_dof_map().add_coupling_functor(*gf);
    libMesh::out << "Added DOF coupling sparsity functor for mesh '" << name << "'." << std::endl;

    mesh->prepare_for_use(true);

    mesh->print_info();

    equation_system->initialize_physics();

    return true;
}
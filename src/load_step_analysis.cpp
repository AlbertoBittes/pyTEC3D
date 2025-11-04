#include "load_step_analysis.h"

#include <algorithm>
#include <filesystem>
#include <set>

#include <libmesh/libmesh.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>

#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"
#include "binary_reader.h"
#include "element_builder.h"
#include "element_types_factory.h"
#include "globals.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "time_step.h"
#include "vtkio_updated.h"

namespace
{
    bool has_bin_extension(const std::filesystem::path &path)
    {
        return path.has_extension() && path.extension() == ".bin";
    }
}

LoadStepAnalysis::LoadStepAnalysis(ApecMEquationSystems *eq_systems,
                                   std::shared_ptr<FEMDataRepository> model_data_repo)
    : ApecMAnalysis(),
      FEParentObjectBase<ApecMEquationSystems>(eq_systems)
{
    if (model_data_repo)
    {
        model_data_handle = model_data_repo;
        session_handle = model_data_repo->get_object<Session>();
    }
}

void LoadStepAnalysis::set_model_repository(const std::shared_ptr<FEMDataRepository> &repo)
{
    model_data_handle = repo;
    if (repo)
        session_handle = repo->get_object<Session>();
}

void LoadStepAnalysis::set_session(const std::shared_ptr<Session> &session_ptr)
{
    session_handle = session_ptr;
}

std::shared_ptr<FEMDataRepository> LoadStepAnalysis::model_repository() const
{
    return model_data_handle.lock();
}

std::shared_ptr<Session> LoadStepAnalysis::session_ptr() const
{
    if (auto session = session_handle.lock())
        return session;
    if (auto repo = model_repository())
        return repo->get_object<Session>();
    return nullptr;
}

std::shared_ptr<ApecMEquationSystems> LoadStepAnalysis::equation_systems_ptr() const
{
    auto *raw_ptr = this->get_parent_object();
    if (!raw_ptr)
        return nullptr;
    return std::shared_ptr<ApecMEquationSystems>(raw_ptr, [](ApecMEquationSystems *) {});
}

std::filesystem::path LoadStepAnalysis::load_steps_folder() const
{
    auto session = session_ptr();
    if (!session)
        return {};

    const auto &folder = session->get_path("load_steps_folderpath");
    return folder;
}

std::vector<std::filesystem::path> LoadStepAnalysis::list_load_step_files() const
{
    std::vector<std::filesystem::path> files;
    const auto folder = load_steps_folder();
    if (folder.empty() || !std::filesystem::exists(folder))
        return files;

    for (const auto &entry : std::filesystem::directory_iterator(folder))
    {
        if (!entry.is_regular_file())
            continue;
        const auto &path = entry.path();
        if (!has_bin_extension(path))
            continue;
        files.push_back(path);
    }

    std::sort(files.begin(), files.end());
    return files;
}

std::size_t LoadStepAnalysis::load_step_file_count() const
{
    return list_load_step_files().size();
}

std::shared_ptr<TimeStep> LoadStepAnalysis::read_time_step(const std::filesystem::path &filepath)
{
    if (filepath.empty())
    {
        libMesh::out << "Load step path is empty." << std::endl;
        return nullptr;
    }

    if (!model_repository())
    {
        libMesh::out << "Model repository not available for load step analysis." << std::endl;
        return nullptr;
    }

    if (!std::filesystem::exists(filepath))
    {
        libMesh::out << "Load step file '" << filepath << "' not found." << std::endl;
        return nullptr;
    }

    auto equation_systems = equation_systems_ptr();
    if (!equation_systems)
    {
        libMesh::out << "Equation systems not available for load step analysis." << std::endl;
        return nullptr;
    }

    try
    {
        auto load_step_repo = std::make_shared<FEMDataRepository>();
        BinaryReader reader;
        reader.read_binary_fem_data(filepath, load_step_repo, equation_systems);
        auto time_step = load_step_repo->get_object<TimeStep>();
        if (!time_step)
        {
            libMesh::out << "Time step data missing in load step binary '" << filepath << "'." << std::endl;
            return nullptr;
        }

        std::shared_ptr<TimeStep> previous_time_step;
        if (active_load_step_repository)
            previous_time_step = active_load_step_repository->get_object<TimeStep>();

        time_step->populate_time_step(previous_time_step);
        set_active_load_step_repository(load_step_repo);
        return time_step;
    }
    catch (const std::exception &ex)
    {
        libMesh::out << "Failed to read load step '" << filepath << "': " << ex.what() << std::endl;
        return nullptr;
    }
    catch (...)
    {
        libMesh::out << "Failed to read load step '" << filepath << "' due to an unknown error." << std::endl;
        return nullptr;
    }

    return nullptr;
}

void LoadStepAnalysis::set_active_load_step_repository(const std::shared_ptr<FEMDataRepository> &repo)
{
    active_load_step_repository = repo;

    auto model_repo = model_repository();
    if (!model_repo || !repo)
    {
        this->ctx_view.reset();
        return;
    }

    auto new_view = std::make_shared<ContextView>(ContextView::from(model_repo, repo));
    this->update_context_views(new_view);
}

void LoadStepAnalysis::configure_services()
{
    auto equation_systems = equation_systems_ptr();
    auto model_repo = model_repository();
    if (!equation_systems || !model_repo || !active_load_step_repository)
        return;

    auto time_step = active_load_step_repository->get_object<TimeStep>();
    if (!time_step)
        return;

    auto element_builder_view = std::make_shared<
        TypedContextView<ElementTypesFactory, MaterialsObject, RealConstantsObject>>(
        TypedContextView<ElementTypesFactory, MaterialsObject, RealConstantsObject>::from(
            model_repo,
            active_load_step_repository));

    auto element_builder = equation_systems->get_element_builder();
    element_builder->update_context_views(element_builder_view);

    auto system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    auto bc_view = std::make_shared<TypedContextView<TimeStep>>(
        TypedContextView<TimeStep>::from(
            active_load_step_repository));
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->update_context_views(bc_view);

    auto assembler_view = std::make_shared<
        TypedContextView<TimeStep, Session>>(
        TypedContextView<TimeStep, Session>::from(
            active_load_step_repository,
            model_repo));

    auto assembler = system.get_assembler();
    assembler->update_context_views(assembler_view);
}

bool LoadStepAnalysis::execute_solve()
{
    auto equation_systems = equation_systems_ptr();
    if (!equation_systems)
        return false;

    auto load_repo = active_load_step_repository;
    if (!load_repo)
    {
        libMesh::out << "No active load step repository available for analysis." << std::endl;
        return false;
    }

    auto time_step = load_repo->get_object<TimeStep>();
    if (!time_step)
    {
        libMesh::out << "Load step does not contain time stepping information." << std::endl;
        return false;
    }

    auto session = session_ptr();
    if (!session)
    {
        libMesh::out << "Session data not available for load step analysis." << std::endl;
        return false;
    }

    configure_services();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();
    if (!assembler)
    {
        libMesh::out << "Assembler not available for linear system." << std::endl;
        return false;
    }

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

        assembler->check_update_system_matrix(&substep);
        assembler->write_element_matrices_and_rhs();

        substep.print_solve_msg();
        libMesh::out << "Solve begins..." << std::endl;

        system.solve();

        system.compute_solution_time_derivative(substep);
        system.get_dof_map().enforce_constraints_exactly(system);
        system.rhs->zero();

        const bool should_write = substep.get_write_results() ||
                                  substep.get_index() == time_step->get_n_substeps() - 1;

        if (should_write)
        {
            libMesh::out << "Computing derived quantities..." << std::endl;
            assembler->write_solution_vectors();

            std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};

            session->set_output_filepath(cumulative_substep_index);
            VTKIO_Updated vtk_output(equation_systems->get_mesh(),
                                     &equation_systems->get_system<ApecMExplicitSystem>("ElementSolutionSystem"));
            vtk_output.write_equation_systems(session->get_path("final_output_path"),
                                              *equation_systems,
                                              &system_names);
            ++cumulative_substep_index;
        }

        *system.old_local_solution = *system.current_local_solution;

        libMesh::out << "Time substep " << substep.get_index() + 1
                     << " of " << time_step->get_n_substeps()
                     << " completed." << std::endl;
    }

    return true;
}

std::filesystem::path LoadStepAnalysis::resolve_load_step_path(const std::string &filename) const
{
    auto folder = load_steps_folder();
    if (folder.empty())
        return {};

    if (filename.empty())
        return {};

    auto candidate = std::filesystem::path(filename);
    if (!candidate.is_absolute())
        candidate = folder / candidate;

    if (!candidate.has_extension())
        candidate.replace_extension(".bin");

    return candidate;
}

std::filesystem::path LoadStepAnalysis::resolve_default_base_name() const
{
    auto session = session_ptr();
    if (!session)
        return {};

    const auto folder = load_steps_folder();
    if (folder.empty())
        return {};

    auto token = session->get_mesh_filename_without_extension();
    if (token.empty())
        token = "load_step";

    std::filesystem::path base = folder / (token + "_loadstep");
    return base;
}

bool LoadStepAnalysis::solve(const std::string &filename)
{
    const auto path = resolve_load_step_path(filename);
    if (path.empty())
    {
        libMesh::out << "Provide a load step filename to solve." << std::endl;
        return false;
    }

    if (!read_time_step(path))
        return false;

    return execute_solve();
}

bool LoadStepAnalysis::solve_transient(const std::string &base_name)
{
    std::filesystem::path base;
    if (base_name.empty())
        base = resolve_default_base_name();
    else
    {
        base = resolve_load_step_path(base_name);
        base = base.parent_path() / base.stem();
    }

    if (base.empty())
    {
        libMesh::out << "Unable to determine base name for transient solve." << std::endl;
        return false;
    }

    unsigned int step_index = 1;
    bool solved_any = false;
    while (true)
    {
        auto next_path = base;
        next_path.replace_filename(base.stem().string() + "_" + std::to_string(step_index));
        next_path.replace_extension(".bin");

        if (!std::filesystem::exists(next_path))
            break;

        libMesh::out << "Solving transient load step '" << next_path << "'." << std::endl;
        if (!read_time_step(next_path))
            return false;

        if (!execute_solve())
            return false;

        solved_any = true;
        ++step_index;
    }

    if (!solved_any)
    {
        libMesh::out << "No transient load steps found for base '" << base << "'." << std::endl;
        return false;
    }

    return true;
}

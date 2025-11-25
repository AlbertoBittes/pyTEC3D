#include "static_analysis.h"
#include "binary_reader.h"
#include "apec_m_system_objects.h"
#include "vtkio_updated.h"

#include <set>

#include <libmesh/libmesh.h>

StaticAnalysis::StaticAnalysis(APEC_M *program,
                               std::shared_ptr<AnalysisBase::ContextView> context_view,
                               Type type,
                               std::string analysis_name)
    : AnalysisBase(type, program, std::move(analysis_name), std::move(context_view))
{
}

StaticAnalysis::StaticAnalysis(Type type,
                               APEC_M *program,
                               std::shared_ptr<AnalysisBase::ContextView> context_view,
                               std::string analysis_name)
    : AnalysisBase(type, program, std::move(analysis_name), std::move(context_view))
{
}

bool StaticAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "StaticAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

bool DOFInterpolationAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "DOFInterpolationAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

bool StaticAnalysis::read_grid_data(const std::string &filepath)
{
    auto equation_systems = active_equation_system_ptr();
    if (!equation_systems)
    {
        libMesh::out << "No active equation systems available for grid data import." << std::endl;
        return false;
    }

    if (filepath.empty())
    {
        libMesh::out << "Grid data filepath is empty." << std::endl;
        return false;
    }

    try
    {
        BinaryReader reader;
        reader.read_binary_fem_data(filepath, nullptr, equation_systems);
    }
    catch (const std::exception &ex)
    {
        libMesh::out << "Failed to read grid data '" << filepath << "': " << ex.what() << std::endl;
        return false;
    }
    catch (...)
    {
        libMesh::out << "Failed to read grid data '" << filepath << "' due to an unknown error." << std::endl;
        return false;
    }

    return true;
}

std::shared_ptr<TimeStep> StaticAnalysis::read_time_step(const std::filesystem::path &filepath)
{
    if (filepath.empty())
    {
        libMesh::out << "Load step path is empty." << std::endl;
        return nullptr;
    }

    if (!std::filesystem::exists(filepath))
    {
        libMesh::out << "Load step file '" << filepath << "' not found." << std::endl;
        return nullptr;
    }

    auto equation_systems = this->active_equation_system_ptr();
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

        auto previous_time_step = this->context<TimeStep>();
        time_step->populate_time_step(previous_time_step);

        auto settings = load_step_repo->get_object<ApecMSettingsObject>();
        update_time_data_context(time_step, settings);
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

std::filesystem::path StaticAnalysis::resolve_load_step_path(const std::string &filename) const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("load_steps_folderpath");
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

std::filesystem::path StaticAnalysis::resolve_default_base_name() const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("load_steps_folderpath");
    if (folder.empty())
        return {};

    auto token = session->get_mesh_filename_without_extension();
    if (token.empty())
        token = "load_step";

    std::filesystem::path base = folder / (token + "_loadstep");
    return base;
}

void StaticAnalysis::configure_services() const
{
    auto *equation_systems = this->active_equation_system();

    auto analysis_ctx = this->get_context_view();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto bc_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep>() : nullptr;
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->update_context_views(bc_view);

    auto assembler_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep, Session>() : nullptr;

    auto assembler = system.get_assembler();
    assembler->update_context_views(assembler_view);
}

bool StaticAnalysis::execute_solve()
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return false;

    auto time_step = this->context<TimeStep>();
    if (!time_step)
    {
        libMesh::out << "Load step does not contain time stepping information." << std::endl;
        return false;
    }

    configure_services();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    for (auto &substep : *time_step)
    {
        prepare_time_substep_solve(substep);

        if (!run_picard_iterations(substep))
            return false;

        print_output(substep);

        *system.old_local_solution = *system.current_local_solution;
    }

    return true;
}

void StaticAnalysis::prepare_time_substep_solve(Substep &substep)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();
    if (!assembler)
    {
        libMesh::out << "Assembler not available for linear system." << std::endl;
        return;
    }

    equation_systems->apply_ramped_effects(substep);

    system.get_dof_map().swap_dof_constraints();

    system.reinit_constraints();

    assembler->write_element_matrices_and_rhs();
}

void StaticAnalysis::solve_equation_systems(Substep &substep)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();
    if (!assembler)
    {
        libMesh::out << "Assembler not available for linear system." << std::endl;
        return;
    }

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.compute_solution_time_derivative(substep);
    assembler->write_solution_vectors(true);
}

void StaticAnalysis::print_output(Substep &substep)
{
    auto *equation_systems = this->active_equation_system();

    auto time_step = this->context<TimeStep>();

    auto session = this->context<Session>();

    auto assembler = equation_systems ? equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem").get_assembler() : nullptr;

    const bool should_write = substep.get_write_results() ||
                              substep.get_index() == substep.get_parent_object()->get_n_substeps() - 1;

    if (!should_write || !session || !equation_systems || !assembler)
        return;

    std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};

    session->set_output_filepath(cumulative_substeps_completed);
    VTKIO_Updated vtk_output(equation_systems->get_mesh(),
                             &equation_systems->get_system<ApecMExplicitSystem>("ElementSolutionSystem"));
    vtk_output.write_equation_systems(session->get_path("final_output_path"),
                                      *equation_systems,
                                      &system_names);
    ++cumulative_substeps_completed;
}

bool StaticAnalysis::run_picard_iterations(Substep &substep)
{
    // Normal transient analysis does not use Picard iterations, always converges in one step
    solve_equation_systems(substep);
    return true;
}

bool StaticAnalysis::solve(const std::string &filename)
{
    std::filesystem::path path;
    if (filename.empty())
    {
        auto base = resolve_default_base_name();
        if (base.empty())
        {
            libMesh::out << "Unable to determine default load step base name." << std::endl;
            return false;
        }

        const auto next_index = load_steps_completed_count + 1;
        path = base;
        path.replace_filename(base.stem().string() + "_" + std::to_string(next_index));
        path.replace_extension(".bin");
    }
    else
        path = resolve_load_step_path(filename);

    if (path.empty())
    {
        libMesh::out << "Unable to resolve load step filename for solve request." << std::endl;
        return false;
    }

    if (!std::filesystem::exists(path))
    {
        libMesh::out << "Load step file '" << path << "' not found." << std::endl;
        return false;
    }

    if (!read_time_step(path))
        return false;

    if (!execute_solve())
        return false;

    load_steps_completed_count++;
    return true;
}

bool StaticAnalysis::solve_transient(const std::string &base_name)
{
    libMesh::out << "Analysis '" << name
                 << "' does not support transient solves (requested base name '"
                 << base_name << "')." << std::endl;
    return false;
}

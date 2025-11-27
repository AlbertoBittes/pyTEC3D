#include "static_analysis.h"
#include "binary_reader.h"
#include "apec_m_system_objects.h"
#include "apec_m.h"
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

bool DOFInterpolationAnalysis::solve()
{
    libMesh::out << "DOFInterpolationAnalysis is not implemented yet." << std::endl;
    return false;
}

std::filesystem::path StaticAnalysis::resolve_grid_step_path(const std::string &filename) const
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

std::filesystem::path StaticAnalysis::resolve_time_step_path(const std::string &filename) const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("mesh_folderpath");
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

std::filesystem::path StaticAnalysis::resolve_default_grid_base_name() const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("load_steps_folderpath");
    if (folder.empty())
        return {};

    auto token = session->get_mesh_filename_without_extension();
    if (token.empty())
        token = "grid";

    std::filesystem::path base = folder / (token + "_grid");
    return base;
}

std::filesystem::path StaticAnalysis::resolve_default_time_base_name() const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("mesh_folderpath");
    if (folder.empty())
        return {};

    auto token = session->get_mesh_filename_without_extension();
    if (token.empty())
        token = "timestep";

    std::filesystem::path base = folder / (token + "_timestep");
    return base;
}

void StaticAnalysis::configure_services() const
{
    auto *equation_systems = this->active_equation_system();

    auto analysis_ctx = this->get_context_view();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto bc_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep>() : nullptr;
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->set_context_view(bc_view);

    auto assembler_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep, Session>() : nullptr;

    auto assembler = system.get_assembler();
    assembler->set_context_view(assembler_view);
}

bool StaticAnalysis::execute_solve()
{
    auto apec_m = this->get_parent_object();

    auto model_data_repo = apec_m ? apec_m->model_data() : nullptr;

    this->set_context_view(std::make_shared<ContextView>(ContextView::from(model_data_repo)));

    auto equation_systems = this->active_equation_system_ptr();

    auto time_step = this->context<TimeStep>();

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

bool StaticAnalysis::solve()
{
    if (!execute_solve())
        return false;

    load_steps_completed_count++;
    return true;
}

bool StaticAnalysis::solve_transient()
{
    libMesh::out << "Analysis '" << name << "' does not support transient solves." << std::endl;
    return false;
}

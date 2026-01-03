#include "static_analysis.h"
#include "binary_reader.h"
#include "apec_m_system.h"
#include "apec_m.h"
#include "vtkio_updated.h"
#include "toroidal_projection.h"

#include <optional>
#include <set>

#include <libmesh/libmesh.h>

bool DOFInterpolationAnalysis::solve()
{
    libMesh::out << "DOFInterpolationAnalysis is not implemented yet." << std::endl;
    return false;
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
        {
            print_output(substep);
            return false;
        }
        print_output(substep);
        *system.old_local_solution = *system.current_local_solution;
    }

    return true;
}

void StaticAnalysis::prepare_time_substep_solve(Substep &substep)
{
    auto *equation_systems = this->active_equation_system();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();

    equation_systems->apply_ramped_effects(substep);

    system.get_dof_map().swap_dof_constraints();

    system.reinit_constraints();

    // assembler->write_element_matrices_and_rhs();
}

void StaticAnalysis::solve_equation_systems(Substep &substep, bool nodal_solution)
{
    auto *equation_systems = this->active_equation_system();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    auto assembler = system.get_assembler();

    assembler->write_element_matrices_and_rhs();

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.compute_solution_time_derivative(substep);
    assembler->write_solution_vectors(nodal_solution);
}

void StaticAnalysis::print_output(Substep &substep, bool picard_iteration_print)
{
    libMesh::out << "Writing output for substep " << substep.get_index() << std::endl;

    auto *equation_systems = this->active_equation_system();

    auto time_step = this->context<TimeStep>();

    auto session = this->context<Session>();

    auto assembler = equation_systems ? equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem").get_assembler() : nullptr;

    const bool should_write = substep.get_write_results() ||
                              substep.get_index() == substep.get_parent_object()->get_n_substeps() - 1;

    if (!should_write)
        return;

    const auto base_substep_index = time_step ? time_step->count_substeps_before_current() : 0;
    std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};

    std::optional<unsigned int> picard_iteration;
    if (picard_iteration_print)
        picard_iteration = substep.get_number_picard_iterations();
    else
        picard_iteration = std::nullopt;

    session->set_output_filepath(base_substep_index + substep.get_index(), picard_iteration);
    VTKIO_Updated vtk_output(equation_systems->get_mesh(),
                             equation_systems);
    vtk_output.write_equation_systems(session->get_path("final_output_path"),
                                      *equation_systems,
                                      &system_names);
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

    return true;
}

bool StaticAnalysis::solve_transient()
{
    libMesh::out << "This analysis does not support transient solves." << std::endl;
    return false;
}

#include "static_analysis.h"
#include "apec_m.h"

#include <libmesh/libmesh.h>

bool DOFInterpolationAnalysis::solve()
{
    libMesh::out << "DOFInterpolationAnalysis is not implemented yet." << std::endl;
    return false;
}

void StaticAnalysis::configure_services() const
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->configure_services(this->get_context_view());
}

bool StaticAnalysis::execute_solve()
{
    auto apec_m = this->get_parent_object();
    if (!apec_m)
        return false;

    auto model_data_repo = apec_m ? apec_m->model_data() : nullptr;

    this->set_context_view(std::make_shared<ContextView>(ContextView::from(model_data_repo)));

    const auto analysis_ctx = this->get_context_view();

    auto equation_systems = this->active_equation_system_ptr();
    if (!equation_systems)
        return false;

    auto time_step = this->context<TimeStep>();

    apply_to_all_equation_systems([&](ApecMEquationSystems &eq_sys)
                                  { eq_sys.configure_services(analysis_ctx); });

    for (auto &substep : *time_step)
    {
        apply_to_all_equation_systems([&](ApecMEquationSystems &eq_sys)
                                      { eq_sys.prepare_time_substep_solve(substep); });

        if (!run_picard_iterations(substep))
        {
            apply_to_all_equation_systems([&](ApecMEquationSystems &eq_sys)
                                          { eq_sys.print_output(substep); });
            return false;
        }
        apply_to_all_equation_systems([&](ApecMEquationSystems &eq_sys)
                                      { eq_sys.print_output(substep); });
        apply_to_all_equation_systems([&](ApecMEquationSystems &eq_sys)
                                      { eq_sys.update_old_local_solution(); });
    }

    return true;
}

void StaticAnalysis::prepare_time_substep_solve(Substep &substep)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->prepare_time_substep_solve(substep);
}

void StaticAnalysis::solve_equation_systems(Substep &substep, bool nodal_solution)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->solve_equation_systems(substep, nodal_solution);
}

void StaticAnalysis::print_output(Substep &substep, bool picard_iteration_print)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->print_output(substep, picard_iteration_print);
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

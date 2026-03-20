#include "analysis_base.h"

#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "apec_m.h"

#include <libmesh/libmesh.h>

AnalysisBase::AnalysisBase(
    APEC_M *program,
    std::shared_ptr<ContextView> context_view)
    : FEParentObjectBase<APEC_M>(program),
      FEMService<ElementTypesFactory, TimeStep, Materials, Constants>(std::move(context_view))
{
}

bool AnalysisBase::finish()
{
    auto program = this->get_parent_object();
    auto equation_systems = program->get_active_equation_system();

    equation_systems->zero_all_systems();
    program->model_data_repository->set_object(std::shared_ptr<TimeStep>()); // drop current time step and its history

    return true;
}

void AnalysisBase::configure_services() const
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->configure_services(this->get_context_view());
}

bool AnalysisBase::execute_solve()
{
    auto apec_m = this->get_parent_object();

    this->set_context_view(std::make_shared<ContextView>(ContextView::from(apec_m->model_data_repository)));

    const auto analysis_ctx = this->get_context_view();

    auto equation_systems = this->get_parent_object()->get_active_equation_system();

    auto time_step = this->context<TimeStep>();

    equation_systems->configure_services(analysis_ctx);

    for (auto &substep : *time_step)
    {
        equation_systems->prepare_time_substep_solve(substep);
        equation_systems->print_output(substep);
        equation_systems->update_old_local_solution();
    }

    return true;
}

void AnalysisBase::prepare_time_substep_solve(Substep &substep)
{
}

void AnalysisBase::solve_equation_systems(Substep &substep, bool nodal_solution)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->solve_equation_systems(substep, nodal_solution);
}

void AnalysisBase::print_output(Substep &substep, bool picard_iteration_print)
{
    auto *equation_systems = this->active_equation_system();
    if (!equation_systems)
        return;
    equation_systems->print_output(substep, picard_iteration_print);
}

bool AnalysisBase::solve()
{
    if (!execute_solve())
        return false;

    return true;
}

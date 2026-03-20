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

    equation_systems->zero_all_systems();
    program->model_data_repository->set_object(std::shared_ptr<TimeStep>()); // drop current time step and its history

    return true;
}

bool AnalysisBase::execute_solve()
{
    auto apec_m = this->get_parent_object();

    this->set_context_view(std::make_shared<ContextView>(ContextView::from(apec_m->model_data_repository)));

    const auto analysis_ctx = this->get_context_view();

    auto time_step = this->context<TimeStep>();

    equation_systems->configure_services(analysis_ctx);

    equation_systems->prepare_time_step_solve(*time_step);
    equation_systems->print_output(*time_step);
    equation_systems->update_old_local_solution();

    return true;
}

bool AnalysisBase::solve()
{
    if (!execute_solve())
        return false;

    return true;
}

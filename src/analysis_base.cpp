#include "analysis_base.h"
#include "apec_m_eq_system_object.h"

#include <libmesh/libmesh.h>

AnalysisBase::AnalysisBase(Type type,
                           ApecMEquationSystems *eq_systems,
                           std::string analysis_name)
    : FEParentObjectBase<ApecMEquationSystems>(eq_systems),
      FEMServiceWithOldContext<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session>(),
      analysis_type(type),
      name(std::move(analysis_name))
{
}

bool AnalysisBase::initialize_model_context(const std::shared_ptr<FEMDataRepository> &model_repo)
{
    if (!model_repo)
    {
        libMesh::out << "AnalysisBase: model data repository is not available for context initialization." << std::endl;
        return false;
    }

    auto dummy_time_repo = std::make_shared<FEMDataRepository>();
    dummy_time_repo->set_object(std::make_shared<TimeStep>());

    auto context = std::make_shared<ContextView>(
        ContextView::from(model_repo, dummy_time_repo));

    std::get<std::shared_ptr<TimeStep>>(context->ptrs).reset();

    this->set_context_view(context);
    this->set_old_context_view(nullptr);

    auto element_factory = context->template get<ElementTypesFactory>();
    auto materials = context->template get<MaterialsObject>();
    auto real_constants = context->template get<RealConstantsObject>();
    auto session = context->template get<Session>();
    auto settings = context->template get<ApecMSettingsObject>();

    const bool needs_settings = (analysis_type == Type::APECM_EQUILIBRIUM || analysis_type == Type::APECM_EVOLUTION);
    const bool ready = static_cast<bool>(element_factory && materials && real_constants && session &&
                                         (!needs_settings || (settings && !settings->empty())));
    if (!ready)
        libMesh::out << "AnalysisBase: incomplete model context detected." << std::endl;
    return ready;
}

bool AnalysisBase::solve_transient(const std::string &base_name)
{
    libMesh::out << "Analysis '" << name
                 << "' does not support transient solves (requested base name '"
                 << base_name << "')." << std::endl;
    return false;
}

void AnalysisBase::update_time_step_context(const std::shared_ptr<TimeStep> &time_step)
{
    auto current = this->get_context_view();
    if (!current)
        return;

    auto new_view = std::make_shared<ContextView>(*current);
    std::get<std::shared_ptr<TimeStep>>(new_view->ptrs) = time_step;
    this->update_context_views(new_view);
}

bool AnalysisBase::finish()
{
    reset_load_step_counter();
    if (auto *equation_systems = this->get_parent_object())
        equation_systems->zero_all_systems();
    return true;
}

void AnalysisBase::mark_load_step_completed()
{
    ++load_steps_completed_count;
}

void AnalysisBase::reset_load_step_counter()
{
    load_steps_completed_count = 0;
}

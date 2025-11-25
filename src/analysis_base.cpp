#include "analysis_base.h"
#include "transient_analysis.h"

#include "apec_m_equilibrium_analysis.h"
#include "apec_m_evolution_analysis.h"

#include "apec_m_eq_system_object.h"
#include "apec_m.h"

#include <libmesh/libmesh.h>

AnalysisBase::AnalysisBase(Type type,
                           APEC_M *program,
                           std::string analysis_name,
                           std::shared_ptr<ContextView> context_view)
    : FEParentObjectBase<APEC_M>(program),
      FEMServiceWithOldContext<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>(),
      analysis_type(type),
      name(std::move(analysis_name))
{
    if (context_view)
        initialize_context(std::move(context_view));
}

std::shared_ptr<ApecMEquationSystems> AnalysisBase::active_equation_system_ptr() const
{
    auto program = this->get_parent_object();
    if (!program)
        return nullptr;
    return program->get_active_equation_system();
}

ApecMEquationSystems *AnalysisBase::active_equation_system() const
{
    auto ptr = active_equation_system_ptr();
    return ptr.get();
}

std::shared_ptr<AnalysisBase::ContextView> AnalysisBase::create_context_view(const std::shared_ptr<FEMDataRepository> &model_repo)
{
    if (!model_repo)
    {
        libMesh::out << "AnalysisBase: model data repository is not available for context initialization." << std::endl;
        return nullptr;
    }

    auto dummy_time_repo = std::make_shared<FEMDataRepository>();
    dummy_time_repo->set_object(std::make_shared<TimeStep>());
    auto dummy_settings_repo = std::make_shared<FEMDataRepository>();
    dummy_settings_repo->set_object(std::make_shared<ApecMSettingsObject>());

    auto context = std::make_shared<ContextView>(
        ContextView::from(model_repo, dummy_time_repo, dummy_settings_repo));

    std::get<std::shared_ptr<TimeStep>>(context->ptrs).reset();
    std::get<std::shared_ptr<ApecMSettingsObject>>(context->ptrs).reset();

    return context;
}

bool AnalysisBase::initialize_context(std::shared_ptr<ContextView> context_view)
{
    if (!context_view)
    {
        libMesh::out << "AnalysisBase: context view is not available for initialization." << std::endl;
        context_initialized = false;
        return false;
    }

    this->set_context_view(std::move(context_view));
    this->set_old_context_view(nullptr);

    auto context = this->get_context_view();
    auto element_factory = context->template get<ElementTypesFactory>();
    auto materials = context->template get<MaterialsObject>();
    auto real_constants = context->template get<RealConstantsObject>();
    auto session = context->template get<Session>();
    const bool ready = static_cast<bool>(element_factory && materials && real_constants && session);
    if (!ready)
        libMesh::out << "AnalysisBase: incomplete model context detected." << std::endl;
    context_initialized = ready;
    return ready;
}

void AnalysisBase::update_time_data_context(const std::shared_ptr<TimeStep> &time_step,
                                            const std::shared_ptr<ApecMSettingsObject> &settings)
{
    auto current = this->get_context_view();
    if (!current)
        return;

    auto new_view = std::make_shared<ContextView>(*current);
    std::get<std::shared_ptr<TimeStep>>(new_view->ptrs) = time_step;
    std::get<std::shared_ptr<ApecMSettingsObject>>(new_view->ptrs) = settings;
    this->update_context_views(new_view);
}

bool AnalysisBase::finish()
{
    load_steps_completed_count = 0;
    if (auto equation_systems = this->active_equation_system_ptr())
        equation_systems->zero_all_systems();
    return true;
}

std::shared_ptr<AnalysisBase> AnalysisBase::build(APEC_M *program,
                                                  const std::shared_ptr<FEMDataRepository> &model_repo,
                                                  std::shared_ptr<ContextView> context_view,
                                                  Type requested_type)
{
    const auto type = requested_type;

    if (!context_view)
        context_view = create_context_view(model_repo);

    std::shared_ptr<AnalysisBase> analysis;

    switch (type)
    {
    case AnalysisBase::Type::STATIC:
        analysis = std::make_shared<StaticAnalysis>(program, context_view);
        break;
    case AnalysisBase::Type::APECM_EQUILIBRIUM:
        analysis = std::make_shared<APECMEquilibriumAnalysis>(program, context_view);
        break;
    case AnalysisBase::Type::APECM_EVOLUTION:
        analysis = std::make_shared<APECMEvolutionAnalysis>(program, context_view);
        break;
    case AnalysisBase::Type::TRANSIENT:
        analysis = std::make_shared<TransientAnalysis>(program, context_view);
        break;
    case AnalysisBase::Type::DOF_INTERPOLATION:
        analysis = std::make_shared<DOFInterpolationAnalysis>(program, context_view);
        break;
    }

    if (analysis && context_view && !analysis->has_valid_context())
    {
        libMesh::out << "Failed to initialize analysis context for '" << analysis->get_name() << "'." << std::endl;
        return nullptr;
    }

    return analysis;
}

#include "analysis_base.h"
#include "transient_analysis.h"

#include "apec_m_equilibrium_analysis.h"
#include "apec_m_evolution_analysis.h"

#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"
#include "apec_m.h"

#include <libmesh/libmesh.h>

AnalysisBase::AnalysisBase(Type type,
                           APEC_M *program,
                           std::string analysis_name,
                           std::shared_ptr<ContextView> context_view)
    : FEParentObjectBase<APEC_M>(program),
      FEMService<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>(),
      analysis_type(type),
      name(std::move(analysis_name))
{
    if (context_view)
        this->set_context_view(std::move(context_view));
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

bool AnalysisBase::finish()
{
    load_steps_completed_count = 0;
    if (auto equation_systems = this->active_equation_system_ptr())
        equation_systems->zero_all_systems();

    if (auto program = this->get_parent_object())
    {
        auto model_repo = program->model_data();
        if (model_repo)
            model_repo->set_object(std::shared_ptr<TimeStep>()); // drop current time step and its history
    }
    return true;
}

std::shared_ptr<AnalysisBase> AnalysisBase::build(APEC_M *program,
                                                  std::shared_ptr<ContextView> context_view,
                                                  Type requested_type)
{
    const auto type = requested_type;

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

    return analysis;
}

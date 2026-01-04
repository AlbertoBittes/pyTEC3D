#include "analysis_base.h"
#include "transient_analysis.h"

#include "apec_m_equilibrium_analysis.h"
#include "apec_m_evolution_analysis.h"

#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "apec_m.h"
#include "session.h"

#include <libmesh/libmesh.h>

AnalysisBase::AnalysisBase(
    APEC_M *program,
    std::shared_ptr<ContextView> context_view)
    : FEParentObjectBase<APEC_M>(program),
      FEMService<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>(std::move(context_view))
{
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

void AnalysisBase::apply_to_all_equation_systems(const std::function<void(ApecMEquationSystems &)> &func) const
{
    auto program = this->get_parent_object();
    if (!program)
        return;

    const auto all_equation_systems = program->get_all_equation_systems_with_names();
    const auto model_data = program->model_data();
    const auto session = model_data ? model_data->get_object<Session>() : nullptr;
    const auto initial_active_mesh = session ? session->get_active_mesh() : std::string{};

    for (const auto &entry : all_equation_systems)
    {
        if (!entry.second)
            continue;
        if (!entry.first.empty())
            program->set_active_mesh(entry.first);
        func(*entry.second);
    }

    if (!initial_active_mesh.empty())
        program->set_active_mesh(initial_active_mesh);
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

#include "apec_m_analysis.h"

#include <libmesh/libmesh.h>

StaticAnalysis::StaticAnalysis(ApecMEquationSystems *eq_systems,
                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : AnalysisBase(Type::STATIC, eq_systems, "StaticAnalysis", std::move(context_view))
{
}

bool StaticAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "StaticAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

DOFInterpolationAnalysis::DOFInterpolationAnalysis(ApecMEquationSystems *eq_systems,
                                                   std::shared_ptr<AnalysisBase::ContextView> context_view)
    : AnalysisBase(Type::DOF_INTERPOLATION, eq_systems, "DOFInterpolationAnalysis", std::move(context_view))
{
}

bool DOFInterpolationAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "DOFInterpolationAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

AnalysisBase::Type AnalysisBase::determine_analysis_type(const std::shared_ptr<FEMDataRepository> &model_repo)
{
    if (model_repo && model_repo->get_object<DofInterpolationSettingsObject>())
        return AnalysisBase::Type::DOF_INTERPOLATION;

    if (model_repo)
    {
        auto settings = model_repo->get_object<ApecMSettingsObject>();
        if (settings && !settings->empty())
        {
            if (settings->analysis_mode() == ApecMSettingsObject::AnalysisMode::EVOLUTION)
                return AnalysisBase::Type::APECM_EVOLUTION;
            else
                return AnalysisBase::Type::APECM_EQUILIBRIUM;
        }
    }

    return AnalysisBase::Type::TRANSIENT;
}

std::shared_ptr<AnalysisBase> AnalysisBase::build(ApecMEquationSystems *eq_systems,
                                                  const std::shared_ptr<FEMDataRepository> &model_repo,
                                                  std::shared_ptr<ContextView> context_view)
{
    const auto type = determine_analysis_type(model_repo);

    if (!context_view)
        context_view = create_context_view(model_repo);

    std::shared_ptr<AnalysisBase> analysis;

    switch (type)
    {
    case AnalysisBase::Type::STATIC:
        analysis = std::make_shared<StaticAnalysis>(eq_systems, context_view);
        break;
    case AnalysisBase::Type::APECM_EQUILIBRIUM:
        analysis = std::make_shared<APECMEquilibriumAnalysis>(eq_systems, context_view);
        break;
    case AnalysisBase::Type::APECM_EVOLUTION:
        analysis = std::make_shared<APECMEvolutionAnalysis>(eq_systems, context_view);
        break;
    case AnalysisBase::Type::TRANSIENT:
        analysis = std::make_shared<TransientAnalysis>(eq_systems, context_view);
        break;
    case AnalysisBase::Type::DOF_INTERPOLATION:
        analysis = std::make_shared<DOFInterpolationAnalysis>(eq_systems, context_view);
        break;
    }

    if (analysis && context_view && !analysis->has_valid_context())
    {
        libMesh::out << "Failed to initialize analysis context for '" << analysis->get_name() << "'." << std::endl;
        return nullptr;
    }

    return analysis;
}

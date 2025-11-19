#include "apec_m_analysis.h"

#include <libmesh/libmesh.h>

StaticAnalysis::StaticAnalysis(ApecMEquationSystems *eq_systems)
    : AnalysisBase(Type::STATIC, eq_systems, "StaticAnalysis")
{
}

bool StaticAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "StaticAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

DOFInterpolationAnalysis::DOFInterpolationAnalysis(ApecMEquationSystems *eq_systems)
    : AnalysisBase(Type::DOF_INTERPOLATION, eq_systems, "DOFInterpolationAnalysis")
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
                                                  const std::shared_ptr<FEMDataRepository> &model_repo)
{
    const auto type = determine_analysis_type(model_repo);

    switch (type)
    {
    case AnalysisBase::Type::STATIC:
        return std::make_shared<StaticAnalysis>(eq_systems);
    case AnalysisBase::Type::APECM_EQUILIBRIUM:
        return std::make_shared<APECMEquilibriumAnalysis>(eq_systems);
    case AnalysisBase::Type::APECM_EVOLUTION:
        return std::make_shared<APECMEvolutionAnalysis>(eq_systems);
    case AnalysisBase::Type::TRANSIENT:
        return std::make_shared<TransientAnalysis>(eq_systems);
    case AnalysisBase::Type::DOF_INTERPOLATION:
        return std::make_shared<DOFInterpolationAnalysis>(eq_systems);
    }

    return nullptr;
}

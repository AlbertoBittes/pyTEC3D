#include "apec_m_analysis.h"

#include <libmesh/libmesh.h>

#include "transient_analysis.h"

StaticAnalysis::StaticAnalysis(ApecMEquationSystems *eq_systems)
    : AnalysisBase(Type::Static, eq_systems, "StaticAnalysis")
{
}

bool StaticAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "StaticAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

APECMAnalysis::APECMAnalysis(ApecMEquationSystems *eq_systems)
    : AnalysisBase(Type::APECM, eq_systems, "APECMAnalysis")
{
}

bool APECMAnalysis::solve(const std::string &identifier)
{
    libMesh::out << "APECMAnalysis is not implemented yet (requested identifier '"
                 << identifier << "')." << std::endl;
    return false;
}

DOFInterpolationAnalysis::DOFInterpolationAnalysis(ApecMEquationSystems *eq_systems)
    : AnalysisBase(Type::DOFInterpolation, eq_systems, "DOFInterpolationAnalysis")
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
        return AnalysisBase::Type::DOFInterpolation;

    if (model_repo && model_repo->get_object<ApecMSettingsObject>())
        return AnalysisBase::Type::APECM;

    return AnalysisBase::Type::Transient;
}

std::shared_ptr<AnalysisBase> AnalysisBase::build(ApecMEquationSystems *eq_systems,
                                                  const std::shared_ptr<FEMDataRepository> &model_repo)
{
    const auto type = determine_analysis_type(model_repo);

    switch (type)
    {
    case AnalysisBase::Type::Static:
        return std::make_shared<StaticAnalysis>(eq_systems);
    case AnalysisBase::Type::APECM:
        return std::make_shared<APECMAnalysis>(eq_systems);
    case AnalysisBase::Type::Transient:
        return std::make_shared<TransientAnalysis>(eq_systems);
    case AnalysisBase::Type::DOFInterpolation:
        return std::make_shared<DOFInterpolationAnalysis>(eq_systems);
    }

    return nullptr;
}

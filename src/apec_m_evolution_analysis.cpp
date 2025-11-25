#include "apec_m_evolution_analysis.h"

#include <filesystem>

#include <libmesh/libmesh.h>

APECMEvolutionAnalysis::APECMEvolutionAnalysis(APEC_M *program,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : StaticAnalysis(AnalysisBase::Type::APECM_EVOLUTION, program, context_view, "APECMEvolutionAnalysis"),
      TransientAnalysis(AnalysisBase::Type::APECM_EVOLUTION, program, "APECMEvolutionAnalysis", context_view),
      APECMEquilibriumAnalysis(program, context_view)
{
    set_name("APECMEvolutionAnalysis");
}

bool APECMEvolutionAnalysis::solve_transient(const std::string &base_name)
{
    return this->solve(base_name);
}

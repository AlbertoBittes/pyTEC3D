#include "apec_m_evolution_analysis.h"

#include <filesystem>

#include <libmesh/libmesh.h>

APECMEvolutionAnalysis::APECMEvolutionAnalysis(APEC_M *program,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : StaticAnalysis(AnalysisBase::Type::APECM_EVOLUTION, program, context_view, "APECMEvolutionAnalysis"),
      APECMEquilibriumAnalysis(program,
                               context_view,
                               AnalysisBase::Type::APECM_EVOLUTION,
                               "APECMEvolutionAnalysis")
{
    set_name("APECMEvolutionAnalysis");
}

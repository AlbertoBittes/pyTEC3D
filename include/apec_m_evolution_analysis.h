#ifndef APEC_M_EVOLUTION_ANALYSIS_H
#define APEC_M_EVOLUTION_ANALYSIS_H

#include "apec_m_equilibrium_analysis.h"

class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis, public TransientAnalysis
{
public:
    explicit APECMEvolutionAnalysis(APEC_M *program = nullptr,
                                    std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
};

#endif // APEC_M_EVOLUTION_ANALYSIS_H

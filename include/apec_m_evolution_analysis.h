#ifndef APEC_M_EVOLUTION_ANALYSIS_H
#define APEC_M_EVOLUTION_ANALYSIS_H

#include "apec_m_equilibrium_analysis.h"

class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis
{
public:
    explicit APECMEvolutionAnalysis(ApecMEquationSystems *eq_systems = nullptr,
                                    std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);

    bool solve(const std::string &filename) override;
    bool solve_transient(const std::string &base_name) override;

private:
    bool solve_evolution_step(TimeStep &time_step, std::size_t profile_index);
};

#endif // APEC_M_EVOLUTION_ANALYSIS_H

#include "apec_m_equilibrium_analysis.h"

APECMEquilibriumAnalysis::APECMEquilibriumAnalysis(ApecMEquationSystems *eq_systems)
    : TransientAnalysis(AnalysisBase::Type::APECM, eq_systems, "APECMEquilibriumAnalysis")
{
}

APECMEvolutionAnalysis::APECMEvolutionAnalysis(ApecMEquationSystems *eq_systems)
    : APECMEquilibriumAnalysis(eq_systems)
{
    set_name("APECMEvolutionAnalysis");
}

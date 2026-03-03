#ifndef APEC_M_EVOLUTION_ANALYSIS_H
#define APEC_M_EVOLUTION_ANALYSIS_H

#include "apec_m_equilibrium_analysis.h"
#include "coordinate_system.h"

class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis
{
public:
    explicit APECMEvolutionAnalysis(APEC_M *program = nullptr,
                                    std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
    virtual bool picard_loop(Substep &substep) override;
    void set_toroidal_slices();
    void toroidal_average();
    void toroidal_projection_current();
    void set_poloidal_currents();

private:
    const CylindricalCoordinateSystem cylindrical_coordinate_system;
    std::vector<std::vector<dof_id_type>> toroidal_slices;
};

#endif // APEC_M_EVOLUTION_ANALYSIS_H

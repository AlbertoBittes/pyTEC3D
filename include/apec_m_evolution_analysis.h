#ifndef APEC_M_EVOLUTION_ANALYSIS_H
#define APEC_M_EVOLUTION_ANALYSIS_H

#include "apec_m_equilibrium_analysis.h"
#include "toroidal_projection.h"
#include "coordinate_system.h"

class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis
{
public:
    explicit APECMEvolutionAnalysis(APEC_M *program = nullptr,
                                    std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
    virtual bool picard_loop(Substep &substep) override;
    void set_toroidal_projection_map();
    void toroidal_average();
    void toroidal_projection_current();
    void set_poloidal_currents();

private:
    std::unique_ptr<ToroidalProjection> toroidal_projection;
    const CylindricalCoordinateSystem cylindrical_coordinate_system;
};

#endif // APEC_M_EVOLUTION_ANALYSIS_H

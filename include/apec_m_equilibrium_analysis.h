#ifndef APECM_EQUILIBRIUM_ANALYSIS_H
#define APECM_EQUILIBRIUM_ANALYSIS_H

#include "transient_analysis.h"
#include "apec_m_system_objects.h"
#include "element_builder.h"

#include <libmesh/node.h>
#include <libmesh/elem.h>
#include <limits>
#include <utility>
#include <vector>

class APECMEquilibriumAnalysis : public TransientAnalysis
{
public:
    explicit APECMEquilibriumAnalysis(ApecMEquationSystems *eq_systems = nullptr,
                                      std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);

protected:
    void apply_jplasma_gradshafranov();
    void constant_plasma_current_density();
    void find_separatrix();
    bool check_poloidal_flux_convergence();
    void apply_control_currents(std::pair<Real, Real> &correction_currents);
    std::pair<Real, Real> compute_correction_currents();

    Real check_limiter();
    Real check_X_point();

    Real get_psi_magnetic_axis() const { return psi_magnetic_axis; }
    void set_psi_magnetic_axis();

    DenseVector<Real> &get_b_target_axis() { return b_target_axis; };
    void set_b_target_axis();

    void identify_plasma_region();
    void reset_plasma_components();
    void remove_plasma_out_of_X_point();

    bool picard_loop(Substep &substep);

    inline bool check_distance_magnetic_axis(const Node *node)
    {
        return (node->operator-(magnetic_axis)).norm() < min_distance_to_axis();
    };

    std::shared_ptr<ApecMSettingsObject> settings() const;
    const std::vector<Real> &plasma_current_profile() const;
    const std::vector<Real> &beta_p_profile() const;
    Real min_distance_to_axis() const;
    virtual bool run_picard_iterations(Substep &substep) override;

private:
    Real psi_magnetic_axis = 0.0;
    Real psi_min = 0.0;

    Point magnetic_axis = Point();

    DenseVector<Real> b_target_axis;

    bool first_absolute_iteration = true;

    std::vector<dof_id_type> vacuum_chamber_psi_dofs;
    std::vector<Number> prev_it_vacuum_chamber_psi_vector;
    std::vector<Number> curr_it_vacuum_chamber_psi_vector;
};

#endif // APECM_EQUILIBRIUM_ANALYSIS_H

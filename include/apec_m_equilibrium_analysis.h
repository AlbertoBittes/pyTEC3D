#ifndef APECM_EQUILIBRIUM_ANALYSIS_H
#define APECM_EQUILIBRIUM_ANALYSIS_H

#include "transient_analysis.h"
#include "apec_m_system.h"
#include "element_builder.h"

#include <libmesh/node.h>
#include <libmesh/elem.h>
#include <limits>
#include <string>
#include <utility>
#include <vector>

class APECMEquilibriumAnalysis : public TransientAnalysis
{
public:
    using TransientAnalysis::TransientAnalysis;

protected:
    void apply_jplasma_gradshafranov();
    // void constant_plasma_current_density();
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

    virtual bool picard_loop(Substep &substep);

    // inline bool check_distance_magnetic_axis(const Point &point)
    // {
    //     return (point - magnetic_axis).norm() > min_distance_to_axis;
    // };

    virtual bool run_picard_iterations(Substep &substep) override;
    void append_iteration_diagnostics(Substep &substep, bool converged);

protected:
    Real psi_magnetic_axis = 0.0;
    Real psi_min = 0.0;

    Point magnetic_axis = Point();

    DenseVector<Real> b_target_axis = DenseVector<Real>(2);

    bool first_absolute_iteration = true;

    Real min_distance_to_axis = 0.3;

    std::vector<dof_id_type> vacuum_chamber_psi_dofs;
    std::vector<Number> prev_it_vacuum_chamber_psi_vector;
    std::vector<Number> curr_it_vacuum_chamber_psi_vector;

    Real last_psi_limiter = 0.0;
    Real last_psi_x_point = 0.0;
    std::string last_separatrix_configuration = "uninitialized";

    Point last_limiter_point = Point();
    Point last_x_point = Point();
    bool has_last_limiter_point = false;
    bool has_last_x_point = false;

    Real last_poloidal_flux_convergence_metric = std::numeric_limits<Real>::quiet_NaN();
    Real last_poloidal_flux_convergence_threshold = std::numeric_limits<Real>::quiet_NaN();

    Real last_radial_correction_current = 0.0;
    Real last_vertical_correction_current = 0.0;
    Real last_total_plasma_current_before_scaling = 0.0;
    Real last_plasma_current_scaling_lambda = 0.0;
    unsigned int last_plasma_element_count = 0;

    std::unique_ptr<NumericVector<Number>> old_coils_current_vector = nullptr;
    // Real prev_norm = std::numeric_limits<Real>::max();
};

#endif // APECM_EQUILIBRIUM_ANALYSIS_H

#ifndef APEC_M_EVOLUTION_ANALYSIS_H
#define APEC_M_EVOLUTION_ANALYSIS_H

#include "apec_m_equilibrium_analysis.h"
#include "coordinate_system.h"

#include <cstddef>
#include <unordered_map>
#include <vector>

class APECMEvolutionAnalysis : public APECMEquilibriumAnalysis
{
public:
    explicit APECMEvolutionAnalysis(APEC_M *program = nullptr,
                                    std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
    virtual bool picard_loop(Substep &substep) override;
    void set_toroidal_slices();
    void toroidal_average();
    void toroidal_projection_plasma_current();
    void set_poloidal_currents();
    void impressed_eddy_currents_coils_currents();
    void set_plasma_current();

private:
    enum class DofCacheKey
    {
        TAVG_EQ_PSI_DOFS,
        TAVG_EQ_B_RAD_DOFS,
        TAVG_EQ_B_VER_DOFS,
        TAVG_EQ_PSIE_DOFS,
        TAVG_EQ_JS_MAIN_DOFS,
        TAVG_EQ_LOCAL_ELEM_IDS,
        TAVG_EQ_COILS_LOCAL_ELEM_IDS,
        TAVG_EVOL_PSI_DOF_GROUPS,
        TAVG_EVOL_B_RAD_DOF_GROUPS,
        TAVG_EVOL_B_VER_DOF_GROUPS,
        TAVG_EVOL_PSIE_DOF_GROUPS,
        TAVG_EVOL_JT_MAIN_DOF_GROUPS,
        TAVG_EVOL_JT_LOCAL_ELEM_IDS,

        TPROJ_EQ_JS_MAIN_DOFS,
        TPROJ_STATIC_JS_MAIN_DOF_GROUPS,
        TPROJ_STATIC_JS_C1_DOF_GROUPS,
        TPROJ_STATIC_JS_C2_DOF_GROUPS,
        TPROJ_STATIC_LOCAL_ELEM_ID_GROUPS,

        IMPRESSED_EVOL_JT_Z,
        IMPRESSED_STATIC_JS_Z,
        IMPRESSED_EVOL_JT_X,
        IMPRESSED_STATIC_JS_X,
        IMPRESSED_EVOL_JT_Y,
        IMPRESSED_STATIC_JS_Y,

        SET_PLASMA_STATIC_JS_Z,
        SET_PLASMA_EVOL_JS_Z,
        SET_PLASMA_STATIC_JS_X,
        SET_PLASMA_EVOL_JS_X,
        SET_PLASMA_STATIC_JS_Y,
        SET_PLASMA_EVOL_JS_Y
    };

    struct DofCacheKeyHash
    {
        std::size_t operator()(DofCacheKey key) const noexcept
        {
            return static_cast<std::size_t>(key);
        }
    };

    struct DofsCache
    {
        bool is_initialized = false;
        std::unordered_map<DofCacheKey, std::vector<dof_id_type>, DofCacheKeyHash> vectors;
        std::unordered_map<DofCacheKey, std::vector<std::vector<dof_id_type>>, DofCacheKeyHash> vector_groups;
    };

    void cache_slice_dofs();
    const std::vector<dof_id_type> &cached_dofs(DofCacheKey key) const;
    const std::vector<std::vector<dof_id_type>> &cached_dof_groups(DofCacheKey key) const;

    const CylindricalCoordinateSystem cylindrical_coordinate_system;
    std::vector<std::vector<dof_id_type>> toroidal_slices_nodes_;
    std::vector<std::vector<dof_id_type>> toroidal_slices_elems_;
    DofsCache dofs_cache_;
};

#endif // APEC_M_EVOLUTION_ANALYSIS_H

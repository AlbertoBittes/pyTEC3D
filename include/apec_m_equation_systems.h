#ifndef APEC_M_EQUATION_SYSTEMS_H
#define APEC_M_EQUATION_SYSTEMS_H

#include <memory>
#include <tuple>
#include <vector>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>

#include "apec_m_periodic_boundary_factory.h"
#include "fem_service.h"
#include "fem_data_repository.h"
#include "time_step.h"
#include "solver_configuration.h"

using namespace libMesh;

class ApecMMesh;
class ElementTypesFactory;
class APEC_M;
class ApecMSystemBase;
class MaterialsObject;
class RealConstantsObject;
class Session;
class ApecMSettingsObject;
class ApecMEquationSystems
    : public EquationSystems,
      public FEDataObjectAbstract,
      public FEMService<ElementTypesFactory, PetscSolverConfiguration>,
      public FEParentObjectBase<APEC_M>

{
public:
    using AnalysisContextView = TypedContextView<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject, Session, ApecMSettingsObject>;
    using ContextView = TypedContextView<ElementTypesFactory, PetscSolverConfiguration>;

    void apply_ramped_effects(const Substep &time_substep);
    void configure_services(const std::shared_ptr<AnalysisContextView> &analysis_ctx);
    void prepare_time_substep_solve(Substep &substep);
    void solve_equation_systems(Substep &substep, bool nodal_solution = true);
    void print_output(Substep &substep, bool picard_iteration_print = false);
    void update_old_local_solution();
    void configure_systems();
    void initialize_system_contexts();
    void post_initialize_systems();
    void initialize_systems(const std::shared_ptr<FEMDataRepository> &model_data);
    void initialize_physics();
    void clear_preserving_solver_settings();
    void zero_all_systems();
    static std::shared_ptr<ApecMEquationSystems> build(ApecMMesh &mesh, APEC_M *program);

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : EquationSystems(std::forward<Args>(args)...)
    {
    }

private:
    std::shared_ptr<AnalysisContextView> analysis_context_view_;
    using RampVectorTuple = std::tuple<ApecMSystemBase *, NumericVector<Number> *, NumericVector<Number> *, NumericVector<Number> *>;
    std::vector<RampVectorTuple> collect_ramp_vectors();
};

#endif // APEC_M_EQUATION_SYSTEMS_H

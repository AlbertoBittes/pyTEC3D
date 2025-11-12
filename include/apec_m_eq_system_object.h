#ifndef APEC_M_EQ_SYSTEM_OBJECT_H
#define APEC_M_EQ_SYSTEM_OBJECT_H

#include <memory>
#include <vector>
#include <libmesh/equation_systems.h>
#include "globals.h"
#include "element_types_factory.h"
#include "fem_service_abstract.h"
#include "fem_data_repository.h"
#include "time_step.h"
#include "apec_m_analysis.h"
#include "solver_configuration.h"

using namespace libMesh;

class ApecMMesh;

class ApecMEquationSystems
    : public EquationSystems,
      public FEDataObjectAbstract,
      public FEMServiceAbstract<ElementTypesFactory, PetscSolverConfiguration>

{
public:
    using ContextView = TypedContextView<ElementTypesFactory, PetscSolverConfiguration>;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : EquationSystems(std::forward<Args>(args)...)
    {
    }

    void apply_ramped_effects(const Substep &time_substep);
    void configure_systems();
    void initialize_system_contexts();
    void post_initialize_systems();
    void initialize_physics(const std::shared_ptr<FEMDataRepository> &model_data);
    void clear_preserving_solver_settings();
    static std::shared_ptr<ApecMEquationSystems> build(ApecMMesh &mesh);
    void set_analysis(std::shared_ptr<AnalysisBase> analysis_object);
    std::shared_ptr<AnalysisBase> get_analysis() const { return analysis; }

private:
    using RampVectorTuple = std::tuple<NumericVector<Number> *, NumericVector<Number> *, NumericVector<Number> *>;
    std::vector<RampVectorTuple> collect_ramp_vectors();

    std::shared_ptr<AnalysisBase> analysis;
};

#endif // APEC_M_EQ_SYSTEM_OBJECT_H

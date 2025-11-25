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
#include "analysis_base.h"
#include "solver_configuration.h"

using namespace libMesh;

class ApecMMesh;

class ApecMEquationSystems
    : public EquationSystems,
      public FEDataObjectAbstract,
      public FEMServiceAbstract<ElementTypesFactory, PetscSolverConfiguration>,
      public FEParentObjectBase<APEC_M>

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
    void initialize_systems(const std::shared_ptr<FEMDataRepository> &model_data);
    void initialize_physics();
    void clear_preserving_solver_settings();
    void zero_all_systems();
    static std::shared_ptr<ApecMEquationSystems> build(ApecMMesh &mesh, APEC_M *program = nullptr);
private:
    using RampVectorTuple = std::tuple<NumericVector<Number> *, NumericVector<Number> *, NumericVector<Number> *>;
    std::vector<RampVectorTuple> collect_ramp_vectors();
};

#endif // APEC_M_EQ_SYSTEM_OBJECT_H

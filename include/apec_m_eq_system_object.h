#ifndef APEC_M_EQ_SYSTEM_OBJECT_H
#define APEC_M_EQ_SYSTEM_OBJECT_H

#include <memory>
#include <libmesh/equation_systems.h>
#include "globals.h"
#include "element_types_factory.h"
#include "fem_service_abstract.h"
#include "fem_data_repository.h"
#include "time_step.h"

using namespace libMesh;

// Forward declarations for any custom classes
class ElementBuilder;
class ApecMMesh;

class ApecMEquationSystems
    : public EquationSystems,
      public FEMServiceAbstract<ElementTypesFactory>

{
public:
    using ContextView = TypedContextView<ElementTypesFactory>;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : EquationSystems(std::forward<Args>(args)...), FEMServiceAbstract<ElementTypesFactory>(nullptr)
    {
        element_builder = std::make_shared<ElementBuilder>(this, nullptr, nullptr);
    }

    std::shared_ptr<ElementBuilder> get_element_builder() const
    {
        return element_builder;
    }

    void apply_ramped_effects(const Substep &time_substep);
    void configure_systems(const std::shared_ptr<FEMDataRepository> &model_data);
    void initialize_physics();
    void clear_preserving_solver_settings();
    static std::shared_ptr<ApecMEquationSystems> build(ApecMMesh &mesh,
                                                       const std::shared_ptr<FEMDataRepository> &model_data);

private:
    using RampVectorTuple = std::tuple<NumericVector<Number> *, NumericVector<Number> *, NumericVector<Number> *>;
    std::vector<RampVectorTuple> collect_ramp_vectors();

    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};

#endif // APEC_M_EQ_SYSTEM_OBJECT_H

#ifndef APEC_M_EQ_SYSTEM_OBJECT_H
#define APEC_M_EQ_SYSTEM_OBJECT_H

#include <memory>
#include <libmesh/equation_systems.h>
#include "globals.h"
#include "element_types_factory.h"
#include "fem_service_abstract.h"
#include "fem_data_repository.h"

using namespace libMesh;

// Forward declarations for any custom classes
class ElementBuilder;

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

private:
    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};

#endif // APEC_M_EQ_SYSTEM_OBJECT_H

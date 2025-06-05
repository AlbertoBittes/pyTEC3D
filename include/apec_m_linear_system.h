#ifndef APEC_M_LINEAR_SYSTEM_H
#define APEC_M_LINEAR_SYSTEM_H

#include <memory>
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "fem_service_abstract.h"
#include "fem_data_context.h"
#include "element_builder.h"
// Add other necessary includes here

class ApecMLinearSystem : public FEMServiceAbstract<ElementPhysicsContext>, public TransientSystem<LinearImplicitSystem>
{
public:
    using TransientSystem<LinearImplicitSystem>::TransientSystem;
};

class ApecMEquationSystems : public FEMServiceAbstract<ElementPhysicsContext>, public EquationSystems
{
public:
    using EquationSystems::EquationSystems;

    void reinit_element_builder(ElementBuilderContext *ctx)
    {
        element_builder = std::make_shared<ElementBuilder>(ctx);
    }

    ElementBuilder *get_element_builder() const
    {
        if (element_builder)
            return element_builder.get();
        return nullptr;
    }

private:
    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};
#endif // APEC_M_LINEAR_SYSTEM_H
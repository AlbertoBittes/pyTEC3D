#ifndef APEC_M_MODEL_OBJECTS_H
#define APEC_M_MODEL_OBJECTS_H

#include <memory>
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "fem_service_abstract.h"
#include "fem_data_repository.h"
#include "element_builder.h"
#include "assembler.h"
#include "boundary_conditions_object.h"

class ApecMLinearSystem
    : public FEMServiceAbstract<ElementTypeFactory>,
      public TransientSystem<LinearImplicitSystem>
{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args, std::shared_ptr<ContextView> ctx_view = nullptr)
        : TransientSystem<LinearImplicitSystem>(std::forward<Args>(args)...),
          FEMServiceAbstract<ElementTypeFactory>(ctx_view)
    {
        assembler = std::make_shared<AssemblerLinear>();
        boundary_conditions_object = std::make_shared<BoundaryConditionsObject>();
    }

    std::shared_ptr<BoundaryConditionsObject> get_boundary_conditions_object() const
    {
        return boundary_conditions_object;
    }

    std::shared_ptr<AssemblerLinear> get_assembler() const
    {
        return assembler;
    }

private:
    std::shared_ptr<BoundaryConditionsObject> boundary_conditions_object = nullptr;
    std::shared_ptr<AssemblerLinear> assembler = nullptr;
};

class ApecMEquationSystems
    : public FEMServiceAbstract<ElementTypeFactory>,
      public EquationSystems
{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args, std::shared_ptr<ContextView> ctx_view = nullptr)
        : EquationSystems(std::forward<Args>(args)...), FEMServiceAbstract<ElementTypeFactory>(ctx_view)
    {
        element_builder = std::make_shared<ElementBuilder>();
    }

    std::shared_ptr<ElementBuilder> get_element_builder() const
    {
        return element_builder;
    }

private:
    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};

#endif // APEC_M_MODEL_OBJECTS_H
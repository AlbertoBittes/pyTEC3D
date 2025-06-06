#ifndef APEC_M_LINEAR_SYSTEM_H
#define APEC_M_LINEAR_SYSTEM_H

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
    using TransientSystem<LinearImplicitSystem>::TransientSystem;

    // Constructor accepting a shared_ptr to a context view
    // ApecMLinearSystem(std::shared_ptr<ContextView> ctx_view = nullptr)
    //     : FEMServiceAbstract<ElementTypeFactory>(ctx_view)
    // {
    // }

    void set_boundary_conditions_object(std::shared_ptr<BoundaryConditionsObject> obj)
    {
        boundary_conditions_object = obj;
    }
    std::shared_ptr<BoundaryConditionsObject> get_boundary_conditions_object() const
    {
        return boundary_conditions_object;
    }

    void set_assembler(std::shared_ptr<AssemblerLinear> assembler_obj)
    {
        assembler = assembler_obj;
    }
    std::shared_ptr<AssemblerLinear> get_assembler() const
    {
        return assembler;
    }

    // Example usage of context:
    std::shared_ptr<ElementTypeFactory> element_type_factory()
    {
        return context<ElementTypeFactory>();
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
    using EquationSystems::EquationSystems;

    // Constructor accepting a shared_ptr to a context view
    // ApecMEquationSystems(std::shared_ptr<ContextView> ctx_view = nullptr)
    //     : FEMServiceAbstract<ElementTypeFactory>(ctx_view)
    // {
    // }

    void init_element_builder()
    {
        element_builder = std::make_shared<ElementBuilder>();
    }
    std::shared_ptr<ElementBuilder> get_element_builder() const
    {
        return element_builder;
    }

    // Example usage of context:
    std::shared_ptr<ElementTypeFactory> element_type_factory()
    {
        return context<ElementTypeFactory>();
    }

private:
    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};

#endif // APEC_M_LINEAR_SYSTEM_H
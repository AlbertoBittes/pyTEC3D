#ifndef APEC_M_MODEL_OBJECTS_H
#define APEC_M_MODEL_OBJECTS_H

#include <memory>
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "fem_service_abstract.h"
#include "element_builder.h"
#include "assembler.h"
#include "boundary_conditions_object.h"

using namespace libMesh;

class ApecMLinearSystem
    : public FEMServiceAbstract<ElementTypeFactory>,
      public TransientSystem<LinearImplicitSystem>
{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Variadic template constructor: forwards all args to TransientSystem<LinearImplicitSystem>
    template <typename... Args>
    ApecMLinearSystem(Args &&...args)
        : TransientSystem<LinearImplicitSystem>(std::forward<Args>(args)...),
          FEMServiceAbstract<ElementTypeFactory>(nullptr) // or whatever default you want
    {
        assembler = std::make_shared<AssemblerLinear>(this);
        this->attach_assemble_object(*assembler);
        boundary_conditions_object = std::make_shared<BoundaryConditionsObject>(this);
        this->attach_constraint_object(*boundary_conditions_object);
    }

    std::shared_ptr<BoundaryConditionsObject> get_boundary_conditions_object() const
    {
        return boundary_conditions_object;
    }

    std::shared_ptr<AssemblerLinear> get_assembler() const
    {
        return assembler;
    }

    void release_subset_solve()
    {
        if (this->_subset != nullptr)
        {
            delete this->_subset;
            this->_subset = nullptr;
        }
    }

private:
    std::shared_ptr<BoundaryConditionsObject> boundary_conditions_object;
    std::shared_ptr<AssemblerLinear> assembler;
};

class ApecMEquationSystems
    : public FEMServiceAbstract<ElementTypeFactory>,
      public EquationSystems
{
public:
    using ContextView = TypedContextView<ElementTypeFactory>;

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : EquationSystems(std::forward<Args>(args)...), FEMServiceAbstract<ElementTypeFactory>(nullptr)
    {
        element_builder = std::make_shared<ElementBuilder>(nullptr, this);
    }

    std::shared_ptr<ElementBuilder> get_element_builder() const
    {
        return element_builder;
    }

private:
    std::shared_ptr<ElementBuilder> element_builder = nullptr;
};

#endif // APEC_M_MODEL_OBJECTS_H
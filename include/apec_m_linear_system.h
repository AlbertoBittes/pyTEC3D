#ifndef APEC_M_LINEAR_SYSTEM_H
#define APEC_M_LINEAR_SYSTEM_H

#include <memory>
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "fem_service_abstract.h"
#include "fem_data_context.h"
#include "element_builder.h"
#include "assembler.h"
// Add other necessary includes here

class ApecMLinearSystem : public FEMServiceAbstract<ElementPhysicsContext>, public TransientSystem<LinearImplicitSystem>
{
public:
    using TransientSystem<LinearImplicitSystem>::TransientSystem;

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

private:
    std::shared_ptr<BoundaryConditionsObject> boundary_conditions_object = nullptr; // Pointer to boundary conditions object
    std::shared_ptr<AssemblerLinear> assembler = nullptr;                           // Pointer to assembler object
};

class ApecMEquationSystems : public FEMServiceAbstract<ElementPhysicsContext>, public EquationSystems
{
public:
    using EquationSystems::EquationSystems;

    void init_element_builder()
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
#endif // APEC_M_LINEAR_SYSTEM_H
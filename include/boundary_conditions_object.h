#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>
#include "fem_service_abstract.h"
#include "node_bc_object.h"
#include "time_step.h"

using namespace libMesh;

// Now expects a context view with NodeBCObject and TimeStep (shared_ptr-based)
class BoundaryConditionsObject
    : public System::Constraint,
      public FEMServiceAbstract<NodeBCObject, TimeStep>,
      public FEDataObjectBase<System>
{
public:
    using ContextView = TypedContextView<NodeBCObject, TimeStep>;
    using FEDataObjectBase<System>::FEDataObjectBase;

    BoundaryConditionsObject(std::shared_ptr<ContextView> ctx_view = nullptr)
        : FEMServiceAbstract<NodeBCObject, TimeStep>(ctx_view)
    {
    }

    virtual ~BoundaryConditionsObject() = default;

    void constrain() override;

    bool get_update_boundary_conditions() const
    {
        return update_boundary_conditions;
    }

    void check_update_boundary_conditions(std::shared_ptr<NodeBCObject> prev_node_bc_object);

    void set_boundary_conditions(std::shared_ptr<NodeBCObject> prev_node_bc_object);

    // Example: access context objects as shared_ptrs
    std::shared_ptr<NodeBCObject> node_bc_object() { return this->context<NodeBCObject>(); }
    std::shared_ptr<TimeStep> time_step() { return this->context<TimeStep>(); }

private:
    bool update_boundary_conditions = false;
    // Other private members related to boundary conditions can be added here
};

#endif // BOUNDARY_CONDITIONS_OBJECT_H
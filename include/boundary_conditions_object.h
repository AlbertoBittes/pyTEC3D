#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>
#include "fem_service_abstract.h"
#include "node_bc_object.h"
#include "time_step.h"

using namespace libMesh;

// The boundary conditions object:
// Checks if it is necessary to update the nodal degrees of freedom imposed by the user
// updates the constraints if necessary. The assembler will then constrain element matrices before putting
// them in the global system matrix
class BoundaryConditionsObject
    : public System::Constraint,
      public FEMServiceWithOldContext<NodeBCObject, TimeStep>,
      public FEParentObjectBase<System>
{
public:
    using ContextView = TypedContextView<NodeBCObject, TimeStep>;
    using FEParentObjectBase<System>::FEParentObjectBase;

    BoundaryConditionsObject(System *sys, std::shared_ptr<ContextView> ctx_view = nullptr, std::shared_ptr<ContextView> old_ctx_view = nullptr)
        : FEParentObjectBase<System>(sys), // Set the parent System
          FEMServiceWithOldContext<NodeBCObject, TimeStep>(ctx_view, old_ctx_view)
    {
    }

    virtual ~BoundaryConditionsObject() = default;

    void constrain() override;

    bool get_update_boundary_conditions() const
    {
        return update_boundary_conditions;
    }

    // Now uses old context view for comparison
    void check_update_boundary_conditions();

    // Uses old context view to decide if constraints need to be reinitialized
    void set_boundary_conditions();

    // Example: access context objects as shared_ptrs
    std::shared_ptr<NodeBCObject> node_bc_object() { return this->context<NodeBCObject>(); }
    std::shared_ptr<TimeStep> time_step() { return this->context<TimeStep>(); }

private:
    bool update_boundary_conditions = false;
    // Other private members related to boundary conditions can be added here
};

#endif // BOUNDARY_CONDITIONS_OBJECT_H
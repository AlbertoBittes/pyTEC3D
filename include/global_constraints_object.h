#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>

#include "fem_service_abstract.h"
#include "time_step.h"

using namespace libMesh;
// The boundary conditions object:
// Checks if it is necessary to update the nodal degrees of freedom imposed by the user
// updates the constraints if necessary. The assembler will then constrain element matrices before putting
// them in the global system matrix
class GlobalConstraintsObject
    : public System::Constraint,
      public FEMServiceWithOldContext<TimeStep>,
      public FEParentObjectBase<System>
{
public:
    using ContextView = TypedContextView<TimeStep>;
    using FEParentObjectBase<System>::FEParentObjectBase;

    GlobalConstraintsObject(System *sys, std::shared_ptr<ContextView> ctx_view = nullptr, std::shared_ptr<ContextView> old_ctx_view = nullptr)
        : FEParentObjectBase<System>(sys), // Set the parent System
          FEMServiceWithOldContext<TimeStep>(ctx_view, old_ctx_view)
    {
    }

    virtual ~GlobalConstraintsObject() = default;

    void constrain() override;

    // Now uses old context view for comparison
    bool check_update_boundary_conditions();
};

#endif // BOUNDARY_CONDITIONS_OBJECT_H
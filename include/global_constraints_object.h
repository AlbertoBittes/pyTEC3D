#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>

#include "fem_service.h"
#include "time_step.h"

#include <map>
#include <string>
#include <vector>

using namespace libMesh;

class DofCouplingsObject;
// The boundary conditions object:
// Checks if it is necessary to update the nodal degrees of freedom imposed by the user
// updates the constraints if necessary. The assembler will then constrain element matrices before putting
// them in the global system matrix
class GlobalConstraintsObject
    : public System::Constraint,
      public FEMService<TimeStep>,
      public FEParentObjectBase<System>
{
public:
  using ContextView = TypedContextView<TimeStep>;
  using FEParentObjectBase<System>::FEParentObjectBase;

  GlobalConstraintsObject(System *sys, std::shared_ptr<ContextView> ctx_view = nullptr)
      : FEParentObjectBase<System>(sys), // Set the parent System
        FEMService<TimeStep>(ctx_view)
  {
  }

  virtual ~GlobalConstraintsObject() = default;

  void constrain() override;

  // Now uses old context view for comparison
  bool check_update_boundary_conditions();

private:
  struct ConstraintEntry
  {
    DofConstraintRow row;
    Real value = 0.0;
  };

  std::map<dof_id_type, ConstraintEntry> build_time_integration_constraints();
  std::map<dof_id_type, ConstraintEntry> build_point_constraints();
  std::shared_ptr<DofCouplingsObject> resolve_dof_couplings();
  std::map<dof_id_type, ConstraintEntry> build_coupling_constraints();
};

#endif // BOUNDARY_CONDITIONS_OBJECT_H

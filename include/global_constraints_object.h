/**
 * @file global_constraints_object.h
 * @brief libMesh constraint handler for Dirichlet BCs, DOF couplings, and time-integration pins.
 *
 * GlobalConstraintsObject derives from libMesh::System::Constraint and is invoked
 * automatically each time libMesh needs to rebuild the constraint map.  It merges
 * three constraint sources:
 *   -# Time-integration homogeneous constraints (static formulation, OPTION4==1).
 *   -# Point/nodal Dirichlet constraints read from VTK point-data arrays.
 *   -# DOF-coupling (periodic boundary) constraints derived from CP_* arrays.
 */
#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>

#include "fe_parent_object.h"

#include <map>
#include <string>
#include <vector>

class DofCouplingsObject;
class ApecMLinearSystem;
/**
 * @brief Applies all boundary and coupling constraints to the DOF system.
 *
 * Registered with libMesh via System::attach_constraint_object().  On each
 * call to constrain(), the three constraint categories are collected and
 * inserted into the DofMap, then the assembler constrains element matrices
 * before scatter-adding them into the global system.
 */
class GlobalConstraintsObject
    : public libMesh::System::Constraint,
      public FEParentObjectBase<ApecMLinearSystem>
{
public:
  using FEParentObjectBase<ApecMLinearSystem>::FEParentObjectBase;

  GlobalConstraintsObject(ApecMLinearSystem *sys)
      : FEParentObjectBase<ApecMLinearSystem>(sys) // Set the parent System
  {
  }

  virtual ~GlobalConstraintsObject() = default;

  /** @brief Called by libMesh to rebuild the active constraint map. */
  void constrain() override;

  /**
   * @brief Returns true if the constraint set has changed since the last time step.
   *
   * When true the system matrix must be refactored; when false it can be reused.
   */
  bool check_update_boundary_conditions();

private:
  struct ConstraintEntry
  {
    libMesh::DofConstraintRow row;
    libMesh::Real value = 0.0;
  };

  /** @brief Homogeneous pin constraints imposed when time-integrated formulation is inactive. */
  std::map<libMesh::dof_id_type, ConstraintEntry> build_time_integration_constraints();

  /** @brief Dirichlet constraints from VTK point-data arrays (nodal BCs). */
  std::map<libMesh::dof_id_type, ConstraintEntry> build_point_constraints();

  /** @brief DOF-coupling constraints derived from CP_* periodic-boundary node pairs. */
  std::map<libMesh::dof_id_type, ConstraintEntry> build_coupling_constraints();
};

#endif // BOUNDARY_CONDITIONS_OBJECT_H

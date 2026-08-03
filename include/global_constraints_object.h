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
#ifndef GLOBAL_CONSTRAINTS_OBJECT_H
#define GLOBAL_CONSTRAINTS_OBJECT_H

#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>

#include "fe_parent_object.h"

#include <map>
#include <string>
#include <unordered_map>
#include <vector>

class DofCouplingsObject;
template <typename ImplicitSys>
class ApecMDOFSystem;
/**
 * @brief Applies all boundary and coupling constraints to the DOF system.
 *
 * Registered with libMesh via System::attach_constraint_object().  On each
 * call to constrain(), the three constraint categories are collected and
 * inserted into the DofMap, then the assembler constrains element matrices
 * before scatter-adding them into the global system.
 */
template <typename ImplSystem>
class GlobalConstraintsObject
    : public libMesh::System::Constraint,
      public FEParentObjectBase<ApecMDOFSystem<ImplSystem>>
{
public:
  using FEParentObjectBase<ApecMDOFSystem<ImplSystem>>::FEParentObjectBase;

  GlobalConstraintsObject(ApecMDOFSystem<ImplSystem> *sys)
      : FEParentObjectBase<ApecMDOFSystem<ImplSystem>>(sys) // Set the parent System
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

  struct DofOwnerInfo
  {
    libMesh::dof_id_type node_id = libMesh::invalid_uint;
    std::string variable;
  };

  /** @brief Homogeneous pin constraints imposed when time-integrated formulation is inactive. */
  std::map<libMesh::dof_id_type, ConstraintEntry> build_time_integration_constraints();

  /** @brief Dirichlet constraints from VTK point-data arrays (nodal BCs). */
  std::map<libMesh::dof_id_type, ConstraintEntry> build_point_constraints();

  /** @brief DOF-coupling constraints derived from CP_* periodic-boundary node pairs. */
  std::map<libMesh::dof_id_type, ConstraintEntry> build_coupling_constraints();

  /// Constraint cache: populated on first build and whenever the time step advances.
  std::map<libMesh::dof_id_type, ConstraintEntry> _constraint_cache;
  std::unordered_map<libMesh::dof_id_type, DofOwnerInfo> _point_constraint_owner_cache;
  bool _has_constraint_cache = false;
};

#endif // GLOBAL_CONSTRAINTS_OBJECT_H

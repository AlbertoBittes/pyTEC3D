/**
 * @file constraints.h
 * @brief libMesh constraint handler for Dirichlet BCs, DOF couplings, and time-integration pins.
 *
 * GlobalConstraintsObject derives from libMesh::System::Constraint and is invoked
 * automatically each time libMesh needs to rebuild the constraint map.  It merges
 * three constraint sources:
 *   -# Time-integration homogeneous constraints: on a static step, zero the DOFs whose fields have no row in K.
 *   -# Point/nodal Dirichlet constraints read from VTK point-data arrays.
 *   -# DOF-coupling (periodic boundary) constraints derived from CP_* arrays.
 */
#pragma once

#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>


#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

class DofCouplingsObject;
template <typename ImplicitSys>
class pyTEC3DDOFSystem;

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
    : public libMesh::System::Constraint
{
public:
  explicit GlobalConstraintsObject(pyTEC3DDOFSystem<ImplSystem> &system);

  virtual ~GlobalConstraintsObject() = default;

  /** @brief Called by libMesh to rebuild the active constraint map. */
  void constrain() override;

  /**
   * @brief Returns true if the constraint set has changed since the last time step.
   *
   * When true the system matrix must be refactored; when false it can be reused. Compared: the constrained DOFs
   * and their rows (masters and coefficients), not the prescribed values. Other DOFs restrict the solve to the new
   * unconstrained DOFs; other rows with the same DOFs force a full factorization (the pattern may change).
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

  /** @brief Caches CP_ID/VARIABLE coupling entries once at construction time. */
  void initialize_coupling_entries_cache();

  /** @brief Builds slave-vs-master midside edge-orientation parity on first constrain() call. */
  void build_edge_orientations_map();

  /// Constraint cache: populated on first build and whenever the time step advances.
  std::map<libMesh::dof_id_type, ConstraintEntry> _constraint_cache;
  std::unordered_map<libMesh::dof_id_type, DofOwnerInfo> _point_constraint_owner_cache;
  std::unique_ptr<std::map<std::pair<libMesh::dof_id_type, libMesh::dof_id_type>, bool>> _edge_orientations_map = nullptr;
  std::vector<int> _coupling_cp_id_values;
  std::vector<std::string> _coupling_variable_values;
  std::unordered_set<int> _edge_coupling_cp_ids;
  bool _has_constraint_cache = false;
  pyTEC3DDOFSystem<ImplSystem> &system_; ///< The DOF system constrained.
};

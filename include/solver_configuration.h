/**
 * @file solver_configuration.h
 * @brief PETSc linear-solver configuration for the APEC-M DOF system.
 *
 * Derives from libMesh::SolverConfiguration and applies solver settings
 * (KSP type, preconditioner, backend) to the PETSc KSP object before each
 * solve.  Currently hardwires a direct LU factorisation via MUMPS.
 */
#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>

#include <pybind11/pybind11.h>

#include "apec_m_visibility.h"

/**
 * @brief Applies PETSc KSP/PC options to the linear solver before each solve.
 *
 * Configures a direct solver (KSPPREONLY + PCLU) backed by MUMPS and
 * prints a one-time configuration summary table to the log.
 */
class APEC_M_HIDDEN PetscSolverConfiguration : public libMesh::SolverConfiguration
{
public:
  PetscSolverConfiguration(
      const pybind11::object &solver_settings_attribute = pybind11::none(),
      libMesh::PetscLinearSolver<libMesh::Number> *petsc_linear_solver = nullptr)
      : _petsc_linear_solver(petsc_linear_solver),
        _solver_settings_attribute(solver_settings_attribute)
  {
    if (_petsc_linear_solver)
      _petsc_linear_solver->set_solver_configuration(*this);
  }

  /** @brief Applies KSP and PC settings to the attached linear solver. */
  void configure_solver() override;

  void set_petsc_linear_solver(libMesh::PetscLinearSolver<libMesh::Number> *petsc_linear_solver)
  {
    _petsc_linear_solver = petsc_linear_solver;
    if (_petsc_linear_solver)
      _petsc_linear_solver->set_solver_configuration(*this);
  }

private:
  libMesh::PetscLinearSolver<libMesh::Number> *_petsc_linear_solver = nullptr; ///< Non-owning pointer to the PETSc solver.
  pybind11::object _solver_settings_attribute;
};

#endif // SOLVER_CONFIGURATION_H

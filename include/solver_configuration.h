/**
 * @file solver_configuration.h
 * @brief PETSc linear-solver configuration for the pyTEC3D DOF system.
 *
 * Derives from libMesh::SolverConfiguration and applies solver settings
 * (KSP type, preconditioner, backend) to the PETSc KSP object before each
 * solve.  Currently hardwires a direct LU factorisation via MUMPS.
 */
#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include <libmesh/linear_solver.h>
#include <libmesh/solver_configuration.h>

#include <pybind11/pybind11.h>

/**
 * @brief Applies PETSc KSP/PC options to the linear solver before each solve.
 *
 * Configures a direct solver (KSPPREONLY + PCLU) backed by MUMPS and
 * prints a one-time configuration summary table to the log.
 */
class PYBIND11_EXPORT PetscSolverConfiguration : public libMesh::SolverConfiguration
{
public:
  ~PetscSolverConfiguration() noexcept override {}

  PetscSolverConfiguration(
      const pybind11::object &solver_settings_attribute = pybind11::none(),
      libMesh::LinearSolver<libMesh::Number> *linear_solver = nullptr)
      : _linear_solver(linear_solver),
        _solver_settings_attribute(solver_settings_attribute)
  {
    if (_linear_solver)
      _linear_solver->set_solver_configuration(*this);
  }

  /** @brief Applies KSP and PC settings to the attached linear solver. */
  void configure_solver() override;

  void set_linear_solver(libMesh::LinearSolver<libMesh::Number> *linear_solver)
  {
    _linear_solver = linear_solver;
    if (_linear_solver)
      _linear_solver->set_solver_configuration(*this);
  }

private:
  libMesh::LinearSolver<libMesh::Number> *_linear_solver = nullptr; ///< Non-owning pointer to the active libMesh linear solver.
  pybind11::object _solver_settings_attribute;
};

#endif // SOLVER_CONFIGURATION_H

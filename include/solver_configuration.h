#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>

#include <pybind11/pybind11.h>

using namespace libMesh;

#if defined(__GNUC__) || defined(__clang__)
#define APEC_M_SOLVERCFG_CLASS_HIDDEN __attribute__((visibility("hidden")))
#else
#define APEC_M_SOLVERCFG_CLASS_HIDDEN
#endif

// This repository stores all data relate to mathematical solver settings, non-linear options (not iplemented yet),
// The user can choose the solver type (direct, iterative), the solver tolerances, and other similar parameters.
class APEC_M_SOLVERCFG_CLASS_HIDDEN PetscSolverConfiguration : public SolverConfiguration
{
public:
  PetscSolverConfiguration(
      const pybind11::object &solver_settings_attribute = pybind11::none(),
      PetscLinearSolver<Number> *petsc_linear_solver = nullptr)
      : _petsc_linear_solver(petsc_linear_solver),
        _solver_settings_attribute(solver_settings_attribute)
  {
    if (_petsc_linear_solver)
      _petsc_linear_solver->set_solver_configuration(*this);
  }

  // Getters and setters for solver settings
  void configure_solver() override;

  void set_petsc_linear_solver(PetscLinearSolver<Number> *petsc_linear_solver)
  {
    _petsc_linear_solver = petsc_linear_solver;
    if (_petsc_linear_solver)
      _petsc_linear_solver->set_solver_configuration(*this);
  }

private:
  PetscLinearSolver<Number> *_petsc_linear_solver = nullptr; // Pointer to the PetscLinearSolver
  pybind11::object _solver_settings_attribute;
};

#undef APEC_M_SOLVERCFG_CLASS_HIDDEN

#endif // SOLVER_CONFIGURATION_H

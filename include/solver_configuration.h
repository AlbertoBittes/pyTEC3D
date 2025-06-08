#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "fe_data_object.h"

using namespace libMesh;

class LoadStep;

class PetscSolverConfiguration : public SolverConfiguration, public SettingsFEDataObject
{
public:
  // Default constructor
  PetscSolverConfiguration(std::ifstream &file, PetscLinearSolver<Number> *petsc_linear_solver = nullptr)
      : SettingsFEDataObject(file),
        _petsc_linear_solver(petsc_linear_solver) {}

  // Getters and setters for solver settings
  void configure_solver() override;

private:
  PetscLinearSolver<Number> *_petsc_linear_solver; // Pointer to the PetscLinearSolver
};

#endif // SOLVER_CONFIGURATION_H
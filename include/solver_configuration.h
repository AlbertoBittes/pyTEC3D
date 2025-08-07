#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "fe_data_object.h"

using namespace libMesh;

class LoadStep;

// This repository stores all data relate to mathematical solver settings, non-linear options (not iplemented yet),
// The user can choose the solver type (direct, iterative), the solver tolerances, and other similar parameters.
class PetscSolverConfiguration : public SolverConfiguration, public SettingsFEDataObject
{
public:
  // Default constructor
  PetscSolverConfiguration(std::ifstream &file, PetscLinearSolver<Number> *petsc_linear_solver = nullptr)
      : SettingsFEDataObject(file),
        _petsc_linear_solver(petsc_linear_solver)
  {
    // PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());
    // SolverType solver_type = Utility::string_to_enum<SolverType>("GMRES");
    // solver->set_solver_type(GMRES);
    // PetscSolverConfiguration solver_configuration(*solver);
    // solver_configuration.configure_solver();
    petsc_linear_solver->set_solver_configuration(*this);
  }

  // Getters and setters for solver settings
  void configure_solver() override;

private:
  PetscLinearSolver<Number> *_petsc_linear_solver; // Pointer to the PetscLinearSolver
};

#endif // SOLVER_CONFIGURATION_H
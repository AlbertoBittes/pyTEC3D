#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include "libmesh/petsc_linear_solver.h"
#include "libmesh/solver_configuration.h"

#include <string>

using namespace libMesh;

class PetscSolverConfiguration : public SolverConfiguration
{
public:
    PetscSolverConfiguration(PetscLinearSolver<Number> &petsc_linear_solver) : _petsc_linear_solver(petsc_linear_solver)
    {
    }

    virtual void configure_solver();

    // The linear solver object that we are configuring
    PetscLinearSolver<Number> &_petsc_linear_solver;
};

#endif // SOLVER_CONFIGURATION_H
#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include "libmesh/petsc_linear_solver.h"
#include "libmesh/solver_configuration.h"
#include "../include/load_step.h"

#include <string>

using namespace libMesh;

class PetscSolverConfiguration : public SolverConfiguration
{
public:
    // Default constructor
    PetscSolverConfiguration()
        : solver_type(0), preconditioner_type(0), verbose(true), _petsc_linear_solver(nullptr), parent_load_step(nullptr) {}

    // Constructor with PetscLinearSolver
    explicit PetscSolverConfiguration(PetscLinearSolver<Number> &petsc_linear_solver)
        : solver_type(0), preconditioner_type(0), verbose(true), _petsc_linear_solver(&petsc_linear_solver), parent_load_step(nullptr) {}

    // Constructor with all attributes
    PetscSolverConfiguration(PetscLinearSolver<Number> &petsc_linear_solver, std::string solver_type,
                             std::string preconditioner_type, bool verbose, std::shared_ptr<LoadStep> parent_load_step = nullptr)
        : solver_type(solver_type), preconditioner_type(preconditioner_type), verbose(verbose),
          _petsc_linear_solver(&petsc_linear_solver), parent_load_step(parent_load_step) {}

    std::shared_ptr<LoadStep> get_parent_load_step() const
    {
        return parent_load_step;
    }

    // Method to configure the solver
    void configure_solver() override;

private:
    std::string solver_type;         // Solver type (e.g., GMRES, CG, etc.)
    std::string preconditioner_type; // Preconditioner type (e.g., ILU, Jacobi, etc.)
    unsigned int n_linear_iterations;
    Real theta;
    bool verbose; // Verbosity flag
    bool reuse_preconditioner;
    PetscLinearSolver<Number> *_petsc_linear_solver; // Pointer to the PetscLinearSolver
    std::shared_ptr<LoadStep> parent_load_step;      // Pointer to the parent LoadStep
};

#endif // SOLVER_CONFIGURATION_H
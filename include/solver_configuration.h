#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include "libmesh/petsc_linear_solver.h"
#include "libmesh/solver_configuration.h"

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
    PetscSolverConfiguration(PetscLinearSolver<Number> &petsc_linear_solver, unsigned int solver_type,
                             unsigned int preconditioner_type, bool verbose, LoadStep *parent_load_step = nullptr)
        : solver_type(solver_type), preconditioner_type(preconditioner_type), verbose(verbose),
          _petsc_linear_solver(&petsc_linear_solver), parent_load_step(parent_load_step) {}

    // Setters
    void set_petsc_linear_solver(PetscLinearSolver<Number> &petsc_linear_solver)
    {
        _petsc_linear_solver = &petsc_linear_solver;
    }

    void set_solver_type(unsigned int type)
    {
        solver_type = type;
    }

    void set_preconditioner_type(unsigned int type)
    {
        preconditioner_type = type;
    }

    void set_verbose(bool verbosity)
    {
        verbose = verbosity;
    }

    void set_parent_load_step(LoadStep *load_step)
    {
        parent_load_step = load_step;
    }

    // Getters
    PetscLinearSolver<Number> *get_petsc_linear_solver() const
    {
        return _petsc_linear_solver;
    }

    unsigned int get_solver_type() const
    {
        return solver_type;
    }

    unsigned int get_preconditioner_type() const
    {
        return preconditioner_type;
    }

    bool is_verbose() const
    {
        return verbose;
    }

    LoadStep *get_parent_load_step() const
    {
        return parent_load_step;
    }

    // Method to configure the solver
    void configure_solver() override;

private:
    unsigned int solver_type;                        // Solver type (e.g., GMRES, CG, etc.)
    unsigned int preconditioner_type;                // Preconditioner type (e.g., ILU, Jacobi, etc.)
    bool verbose;                                    // Verbosity flag
    PetscLinearSolver<Number> *_petsc_linear_solver; // Pointer to the PetscLinearSolver
    LoadStep *parent_load_step;                      // Pointer to the parent LoadStep
};

#endif // SOLVER_CONFIGURATION_H
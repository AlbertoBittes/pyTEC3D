#include "solver_configuration.h"

using namespace libMesh;

void PetscSolverConfiguration::configure_solver()
{
    PetscErrorCode ierr = 0;
    // ierr = KSPSetType(_petsc_linear_solver->ksp(), const_cast<KSPType>(KSPGMRES));
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PCSetType(_petsc_linear_solver->pc(), const_cast<PCType>(PCILU));
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // Initialize PETSc logging
    ierr = PetscLogDefaultBegin();
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
    // Set the solver type to GMRES
    ierr = KSPSetType(_petsc_linear_solver->ksp(), const_cast<KSPType>(KSPPREONLY));
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
    ierr = KSPSetInitialGuessNonzero(_petsc_linear_solver->ksp(), PETSC_FALSE);
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
    // Set the preconditioner type to MUMPS
    ierr = PCSetType(_petsc_linear_solver->pc(), const_cast<PCType>(PCLU));
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
    // ierr = PCSetType(_petsc_linear_solver->pc(), const_cast<PCType>(PCCHOLESKY));
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    ierr = PCFactorSetMatSolverType(_petsc_linear_solver->pc(), MATSOLVERMUMPS);
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
    // Set the verbosity level
    ierr = KSPSetFromOptions(_petsc_linear_solver->ksp());
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
    // Set additional options for verbosity
    ierr = PetscOptionsSetValue(NULL, "-ksp_monitor", "");
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    ierr = PetscOptionsSetValue(NULL, "-ksp_view", "");
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    ierr = PetscOptionsSetValue(NULL, "-ksp_converged_reason", "");
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    ierr = PetscOptionsSetValue(NULL, "-log_view", "");
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_true_residual", "");
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_singular_value", "");
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_solution", "");
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_short", "");
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_range", "");
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_lg_residualnorm", "");
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_lg_true_residualnorm", "");
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    // ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_lg_range", "");
    // CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);
}

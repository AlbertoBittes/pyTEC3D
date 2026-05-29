#include "solver_configuration.h"
#include "apec_m_logging.h"

#include <sstream>
#include <string>

using namespace libMesh;

namespace
{
    bool g_petsc_config_table_printed = false;
}

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

    const KSP ksp = _petsc_linear_solver->ksp();
    const PC pc = _petsc_linear_solver->pc();

    KSPType ksp_type = nullptr;
    ierr = KSPGetType(ksp, &ksp_type);
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    PCType pc_type = nullptr;
    ierr = PCGetType(pc, &pc_type);
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    MatSolverType mat_solver_type = nullptr;
    ierr = PCFactorGetMatSolverType(pc, &mat_solver_type);
    CHKERRABORT(_petsc_linear_solver->comm().get(), ierr);

    const std::string ksp_name = ksp_type ? std::string(ksp_type) : std::string("unknown");
    const std::string pc_name = pc_type ? std::string(pc_type) : std::string("unknown");
    const std::string backend_name = mat_solver_type ? std::string(mat_solver_type) : std::string("petsc-default");
    const std::string solver_class =
        (ksp_name == KSPPREONLY) ? "Direct" : "Iterative";

    if (g_petsc_config_table_printed)
        return;

    std::ostringstream table;
    table << "+------------------------------------------------------------+\n";
    table << "| PETSc Solver Configuration                                 |\n";
    table << "+----------------+-------------------------------------------+\n";
    table << fmt::format("| {:<14} | {:<41} |\n", "Backend", backend_name);
    table << fmt::format("| {:<14} | {:<41} |\n", "Type", solver_class);
    table << fmt::format("| {:<14} | {:<41} |\n", "KSP", ksp_name);
    table << fmt::format("| {:<14} | {:<41} |\n", "Preconditioner", pc_name);
    table << "+------------------------------------------------------------+";
    apec::logging::table_root(table.str());
    g_petsc_config_table_printed = true;
}

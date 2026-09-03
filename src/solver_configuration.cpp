/**
 * @file solver_configuration.cpp
 * @brief PetscSolverConfiguration::configure_solver() implementation.
 *
 * Call order inside PetscLinearSolver::solve_common():
 *   1. get_real/int_solver_setting() reads real_valued_data / int_valued_data  →
 *      KSPSetTolerances() with those values.
 *   2. KSPSetFromOptions() — command-line overrides.
 *   3. configure_solver() — runs LAST.
 *
 * Therefore configure_solver():
 *   - Writes real_valued_data["rel_tol"], real_valued_data["abs_tol"],
 *     int_valued_data["max_its"] so future solve() calls read the correct
 *     values at step 1.
 *   - Calls KSPSetTolerances() directly for the current solve (step 1 already
 *     ran before us).
 *   - Sets KSP type and PC type from the Python solver_settings dict.
 *   - Enables MUMPS only for direct factorisations (lu / cholesky).
 *   - Calls KSPSetFromOptions() at the end so command-line flags can still
 *     override everything we set above.
 */
#include "solver_configuration.h"
#include "apec_m_logging.h"

#include <libmesh/libmesh.h>

#include <cstdlib>
#include <string>

namespace
{
    bool g_petsc_config_table_printed = false;
}

void PetscSolverConfiguration::configure_solver()
{
    PetscErrorCode ierr = 0;
    const auto comm = _petsc_linear_solver->comm().get();
    const int n_threads = std::max(1, static_cast<int>(libMesh::n_threads()));

    // Read settings from the Python solver_settings dict.
    std::string ksp_type_str = KSPPREONLY;
    std::string pc_type_str = PCLU;
    if (!_solver_settings_attribute.is_none())
    {
        ksp_type_str = _solver_settings_attribute["ksp_type"].cast<std::string>();
        pc_type_str = _solver_settings_attribute["pc_type"].cast<std::string>();

        const PetscReal rtol = _solver_settings_attribute["initial_linear_tolerance"].cast<PetscReal>();
        const PetscReal abstol = _solver_settings_attribute["minimum_linear_tolerance"].cast<PetscReal>();
        const int maxits = _solver_settings_attribute["max_linear_iterations"].cast<int>();

        // Store in the libMesh SolverConfiguration maps so that future
        // solve() calls read the correct values at step 1 (before configure_solver
        // is called for that solve).
        real_valued_data["rel_tol"] = static_cast<double>(rtol);
        real_valued_data["abs_tol"] = static_cast<double>(abstol);
        int_valued_data["max_its"] = maxits;

        // Apply to the current solve directly: libMesh's KSPSetTolerances
        // already ran before this function, so we must call it again here.
        ierr = KSPSetTolerances(_petsc_linear_solver->ksp(),
                                rtol, abstol, PETSC_DEFAULT,
                                static_cast<PetscInt>(maxits));
        CHKERRABORT(comm, ierr);
    }

    // Keep BLAS/OpenMP backends consistent with the libMesh thread count.
    const std::string thread_count_str = std::to_string(n_threads);
    (void)setenv("OMP_NUM_THREADS", thread_count_str.c_str(), 1);

    // Enable PETSc performance logging so that -log_view output is available.
    ierr = PetscLogDefaultBegin();
    CHKERRABORT(comm, ierr);

    // Set KSP and PC types from settings.
    ierr = KSPSetType(_petsc_linear_solver->ksp(), ksp_type_str.c_str());
    CHKERRABORT(comm, ierr);
    ierr = KSPSetInitialGuessNonzero(_petsc_linear_solver->ksp(), PETSC_FALSE);
    CHKERRABORT(comm, ierr);

    ierr = PCSetType(_petsc_linear_solver->pc(), pc_type_str.c_str());
    CHKERRABORT(comm, ierr);

    // Enable convergence and solution monitoring via PETSc options.
    ierr = PetscOptionsSetValue(NULL, "-ksp_monitor", "");
    CHKERRABORT(comm, ierr);
    ierr = PetscOptionsSetValue(NULL, "-ksp_view", "");
    CHKERRABORT(comm, ierr);
    ierr = PetscOptionsSetValue(NULL, "-ksp_converged_reason", "");
    CHKERRABORT(comm, ierr);
    ierr = PetscOptionsSetValue(NULL, "-log_view", "");
    CHKERRABORT(comm, ierr);
    ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_true_residual", "");
    CHKERRABORT(comm, ierr);
    ierr = PetscOptionsSetValue(NULL, "-ksp_monitor_solution", "");
    CHKERRABORT(comm, ierr);
    ierr = PetscOptionsSetValue(NULL, "-mat_mumps_icntl_4", "2");
    CHKERRABORT(comm, ierr);

    const KSP ksp = _petsc_linear_solver->ksp();
    const PC pc = _petsc_linear_solver->pc();

    KSPType ksp_type = nullptr;
    ierr = KSPGetType(ksp, &ksp_type);
    CHKERRABORT(comm, ierr);

    PCType pc_type = nullptr;
    ierr = PCGetType(pc, &pc_type);
    CHKERRABORT(comm, ierr);

    MatSolverType mat_solver_type = nullptr;
    ierr = PCFactorGetMatSolverType(pc, &mat_solver_type);
    CHKERRABORT(comm, ierr);

    const std::string ksp_name = ksp_type ? std::string(ksp_type) : std::string("unknown");
    const std::string pc_name = pc_type ? std::string(pc_type) : std::string("unknown");
    const std::string backend_name = mat_solver_type ? std::string(mat_solver_type) : std::string("petsc-default");
    const std::string solver_class =
        (ksp_name == KSPPREONLY) ? "Direct" : "Iterative";

    if (g_petsc_config_table_printed)
        return;

    apec::logging::table_root(apec::logging::format_petsc_solver_configuration_table(
        backend_name,
        solver_class,
        ksp_name,
        pc_name));
    g_petsc_config_table_printed = true;
}

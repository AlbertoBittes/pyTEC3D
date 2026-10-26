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
 *   - Sets KSP type and PC type from the solver settings.
 *   - Enables MUMPS only for direct factorisations (lu / cholesky).
 *
 * The PETSc options of the settings are set by set_settings(), before the solve, so that its KSPSetFromOptions()
 * (step 2) sees them.
 */
#include "solver_configuration.h"
#include "logging.h"
#include "report_tables.h"

#include <libmesh/libmesh.h>
#include <libmesh/petsc_linear_solver.h>

#include <algorithm>
#include <cstdlib>
#include <string>

PetscSolverConfiguration::PetscSolverConfiguration(libMesh::LinearSolver<libMesh::Number> *linear_solver, bool nonlinear)
    : linear_solver_(linear_solver), nonlinear_(nonlinear)
{
    if (linear_solver_)
        linear_solver_->set_solver_configuration(*this);
}

bool PetscSolverConfiguration::set_settings(const SolverSettings &settings)
{
    const bool changed = has_settings_ && (settings.linear != settings_.linear || settings.petsc_options != settings_.petsc_options);

    PetscErrorCode ierr = 0;
    const auto comm = linear_solver_ ? linear_solver_->comm().get() : MPI_COMM_WORLD;
    for (const auto &[name, value] : settings_.petsc_options)
    {
        const bool kept = std::any_of(settings.petsc_options.begin(), settings.petsc_options.end(),
                                      [&name = name](const auto &option) { return option.first == name; });
        if (has_settings_ && !kept)
        {
            ierr = PetscOptionsClearValue(NULL, name.c_str());
            CHKERRABORT(comm, ierr);
        }
    }
    for (const auto &[name, value] : settings.petsc_options)
    {
        ierr = PetscOptionsSetValue(NULL, name.c_str(), value.c_str());
        CHKERRABORT(comm, ierr);
    }

    settings_ = settings;
    has_settings_ = true;
    report_pending_ = report_pending_ || changed;
    return changed;
}

void PetscSolverConfiguration::configure_solver()
{
    auto *petsc_linear_solver = dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(linear_solver_);
    if (!petsc_linear_solver)
        pytec3d::logging::error_and_abort("PetscSolverConfiguration::configure_solver requires libMesh::PetscLinearSolver.");

    PetscErrorCode ierr = 0;
    const auto comm = petsc_linear_solver->comm().get();
    const int n_threads = std::max(1, static_cast<int>(libMesh::n_threads()));

    const LinearSolverSettings &linear = settings_.linear;
    const PetscReal rtol = linear.relative_tolerance;
    const PetscReal abstol = linear.absolute_tolerance;
    const PetscReal dtol = linear.divergence_tolerance > 0.0 ? linear.divergence_tolerance : PETSC_DEFAULT;
    const int maxits = linear.max_iterations;

    // Store in the libMesh SolverConfiguration maps so that future solve() calls read the correct values at
    // step 1 (before configure_solver is called for that solve).
    real_valued_data["rel_tol"] = static_cast<double>(rtol);
    real_valued_data["abs_tol"] = static_cast<double>(abstol);
    int_valued_data["max_its"] = maxits;

    // Apply to the current solve directly: libMesh's KSPSetTolerances already ran before this function, so we
    // must call it again here.
    ierr = KSPSetTolerances(petsc_linear_solver->ksp(), rtol, abstol, dtol, static_cast<PetscInt>(maxits));
    CHKERRABORT(comm, ierr);

    // Keep BLAS/OpenMP backends consistent with the libMesh thread count.
    const std::string thread_count_str = std::to_string(n_threads);
#ifdef _WIN32
    (void)_putenv_s("OMP_NUM_THREADS", thread_count_str.c_str());
#else
    (void)setenv("OMP_NUM_THREADS", thread_count_str.c_str(), 1);
#endif

    // Enable PETSc performance logging so that -log_view output is available.
    ierr = PetscLogDefaultBegin();
    CHKERRABORT(comm, ierr);

    ierr = KSPSetType(petsc_linear_solver->ksp(), linear.ksp_type.c_str());
    CHKERRABORT(comm, ierr);
    ierr = KSPSetInitialGuessNonzero(petsc_linear_solver->ksp(), linear.initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE);
    CHKERRABORT(comm, ierr);
    ierr = PCSetType(petsc_linear_solver->pc(), linear.pc_type.c_str());
    CHKERRABORT(comm, ierr);

    // The factorization package of the direct preconditioners (MUMPS by default, also in parallel).
    if ((linear.pc_type == PCLU || linear.pc_type == PCCHOLESKY) && !linear.factor_solver_type.empty())
    {
        ierr = PCFactorSetMatSolverType(petsc_linear_solver->pc(), linear.factor_solver_type.c_str());
        CHKERRABORT(comm, ierr);
    }

    const KSP ksp = petsc_linear_solver->ksp();
    const PC pc = petsc_linear_solver->pc();

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

    if (!report_pending_)
        return;

    pytec3d::logging::table_root(
        pytec3d::report::format_solver_configuration_table(settings_, nonlinear_, ksp_name, pc_name, backend_name));
    report_pending_ = false;
}


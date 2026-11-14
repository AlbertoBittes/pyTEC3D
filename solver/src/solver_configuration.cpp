/**
 * @file solver_configuration.cpp
 * @brief PetscSolverConfiguration: the solver settings on the PETSc KSP.
 *
 * Order on each solve: tolerances, then KSPSetFromOptions() (the options database: the PETSc options of the settings
 * and the command line), then the KSP and PC types of the settings, which therefore win.
 */
#include "solver_configuration.h"

#include <algorithm>
#include <cstdlib>
#include <string>

#include <omp.h>

#include "logging.h"
#include "report_tables.h"

namespace
{
    void check(PetscErrorCode ierr, const char *what)
    {
        if (ierr)
            pytec3d::logging::error_and_abort(std::string("PETSc error in ") + what + ".");
    }
} // namespace

bool PetscSolverConfiguration::set_settings(const SolverSettings &settings, bool nonlinear)
{
    const bool changed = has_settings_ && (settings.linear != settings_.linear || settings.petsc_options != settings_.petsc_options);
    for (const auto &[name, value] : settings_.petsc_options)
    {
        const bool kept = std::any_of(settings.petsc_options.begin(), settings.petsc_options.end(),
                                      [&name = name](const auto &option) { return option.first == name; });
        if (has_settings_ && !kept)
            check(PetscOptionsClearValue(NULL, name.c_str()), "PetscOptionsClearValue");
    }
    for (const auto &[name, value] : settings.petsc_options)
        check(PetscOptionsSetValue(NULL, name.c_str(), value.c_str()), "PetscOptionsSetValue");
    settings_ = settings;
    has_settings_ = true;
    nonlinear_ = nonlinear;
    report_pending_ = report_pending_ || changed;
    return changed;
}

void PetscSolverConfiguration::configure(KSP ksp, double rtol)
{
    const LinearSolverSettings &linear = settings_.linear;
    const PetscReal dtol = linear.divergence_tolerance > 0.0 ? linear.divergence_tolerance : PETSC_DEFAULT;
    check(KSPSetTolerances(ksp, rtol, linear.absolute_tolerance, dtol, static_cast<PetscInt>(linear.max_iterations)),
          "KSPSetTolerances");

    // Keep the BLAS / OpenMP backends of the factorization consistent with the threads of the solver.
    const std::string threads = std::to_string(std::max(1, omp_get_max_threads()));
#ifdef _WIN32
    (void)_putenv_s("OMP_NUM_THREADS", threads.c_str());
#else
    (void)setenv("OMP_NUM_THREADS", threads.c_str(), 1);
#endif
    // Performance logging, so that -log_view has data.
    check(PetscLogDefaultBegin(), "PetscLogDefaultBegin");

    check(KSPSetFromOptions(ksp), "KSPSetFromOptions");
    check(KSPSetType(ksp, linear.ksp_type.c_str()), "KSPSetType");
    check(KSPSetInitialGuessNonzero(ksp, linear.initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE),
          "KSPSetInitialGuessNonzero");
    PC pc = nullptr;
    check(KSPGetPC(ksp, &pc), "KSPGetPC");
    check(PCSetType(pc, linear.pc_type.c_str()), "PCSetType");
    // The factorization package of the direct preconditioners (MUMPS by default, also in parallel).
    if ((linear.pc_type == PCLU || linear.pc_type == PCCHOLESKY) && !linear.factor_solver_type.empty())
        check(PCFactorSetMatSolverType(pc, linear.factor_solver_type.c_str()), "PCFactorSetMatSolverType");

    if (!report_pending_)
        return;
    KSPType ksp_type = nullptr;
    PCType pc_type = nullptr;
    MatSolverType mat_solver_type = nullptr;
    check(KSPGetType(ksp, &ksp_type), "KSPGetType");
    check(PCGetType(pc, &pc_type), "PCGetType");
    check(PCFactorGetMatSolverType(pc, &mat_solver_type), "PCFactorGetMatSolverType");
    pytec3d::logging::table_root(pytec3d::report::format_solver_configuration_table(
        settings_, nonlinear_, ksp_type ? ksp_type : "unknown", pc_type ? pc_type : "unknown",
        mat_solver_type ? mat_solver_type : "petsc-default"));
    report_pending_ = false;
}

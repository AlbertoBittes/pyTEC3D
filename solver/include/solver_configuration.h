/**
 * @file solver_configuration.h
 * @brief The solver settings of the DOF system and their application to the PETSc linear solver (KSP).
 *
 * The KSP type, preconditioner and tolerances of the settings are set on the KSP before each solve, after
 * KSPSetFromOptions(), so that they win over the options database; the PETSc options of the settings are put in the
 * options database. Direct factorisations (LU, Cholesky) use MUMPS unless the settings say otherwise.
 */
#pragma once

#include <string>
#include <utility>
#include <vector>

#include <petscksp.h>

/** @brief Settings of the PETSc linear solver (KSP), used alone or inside each Newton iteration. */
struct LinearSolverSettings
{
    std::string ksp_type = "preonly";
    std::string pc_type = "lu";
    double relative_tolerance = 1e-6;
    double absolute_tolerance = 1e-10;
    double divergence_tolerance = 0.0; ///< <= 0: the PETSc default.
    int max_iterations = 100;
    bool initial_guess_nonzero = false;
    std::string factor_solver_type = "mumps"; ///< Of the lu and cholesky preconditioners; empty: the PETSc default.
    /// "petsc": the PETSc KSP (the reduced matrix converted); "mfem": MFEM's solvers on the hypre matrix itself (MUMPSSolver
    /// for lu / cholesky; CG, GMRES, FGMRES, MINRES, BiCGSTAB with hypre's BoomerAMG, AMS or Jacobi).
    std::string backend = "petsc";

    bool operator==(const LinearSolverSettings &other) const
    {
        return ksp_type == other.ksp_type && pc_type == other.pc_type && relative_tolerance == other.relative_tolerance &&
               absolute_tolerance == other.absolute_tolerance && divergence_tolerance == other.divergence_tolerance &&
               max_iterations == other.max_iterations && initial_guess_nonzero == other.initial_guess_nonzero &&
               factor_solver_type == other.factor_solver_type && backend == other.backend;
    }
    bool operator!=(const LinearSolverSettings &other) const { return !(*this == other); }
};

/** @brief Settings of the Newton solver of the nonlinear (law) systems (newton_solver.h). */
struct NewtonSettings
{
    unsigned int max_iterations = 50;
    double absolute_residual_tolerance = 1e-35;
    double relative_residual_tolerance = 1e-8;
    double absolute_step_tolerance = 1e-8;
    double relative_step_tolerance = 1e-8;
    bool require_residual_reduction = false;
    bool require_finite_residual = true;
    bool brent_line_search = true;
    /// Picard iteration: the Jacobian of the laws is their secant (the tangent part is dropped), u+ = K_s(u)^-1 f.
    bool picard = false;
    double minimum_step_length = 1e-5;
    double linear_tolerance_multiplier = 1e-3;
    bool continue_after_max_iterations = true;
    bool continue_after_backtrack_failure = false;
    bool verbose = true;
};

/** @brief The solver settings of a solve (the Python solver_settings dict); they may change from one solve to the next. */
struct SolverSettings
{
    LinearSolverSettings linear;
    NewtonSettings nonlinear;
    /// PETSc options set before each linear solve, name ("-ksp_monitor") and value ("" for a flag).
    std::vector<std::pair<std::string, std::string>> petsc_options = {
        {"-ksp_monitor", ""}, {"-ksp_view", ""}, {"-ksp_converged_reason", ""}, {"-log_view", ""},
        {"-ksp_monitor_true_residual", ""}, {"-ksp_monitor_solution", ""}, {"-mat_mumps_icntl_4", "2"}};
};

/**
 * @brief Applies the settings to a KSP before each solve.
 *
 * Prints the configuration table to the log at the first solve and whenever the settings change.
 */
class PetscSolverConfiguration
{
public:
    /**
     * @brief The settings of the next solves.
     *
     * Sets their PETSc options now, and clears those of the previous settings they leave out, so that the
     * KSPSetFromOptions() of the next solve sees them.
     * @return True when the linear solver settings or the PETSc options differ from those of the previous call (not on
     * the first call): a factorization of the previous settings must not be reused.
     */
    bool set_settings(const SolverSettings &settings, bool nonlinear);
    /** @brief Applies the tolerances, the options database and the KSP and PC types to @p ksp. */
    void configure(KSP ksp, double rtol);
    const SolverSettings &settings() const { return settings_; }
    bool nonlinear() const { return nonlinear_; }
    /** @brief True once after new settings: the caller prints the configuration table (the MFEM backend). */
    bool take_report()
    {
        const bool pending = report_pending_;
        report_pending_ = false;
        return pending;
    }

private:
    SolverSettings settings_;
    bool has_settings_ = false;
    bool report_pending_ = true; ///< Print the configuration table at the next configure().
    bool nonlinear_ = false;
};

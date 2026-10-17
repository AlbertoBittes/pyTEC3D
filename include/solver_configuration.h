/**
 * @file solver_configuration.h
 * @brief PETSc linear-solver configuration for the pyTEC3D DOF system.
 *
 * Derives from libMesh::SolverConfiguration and applies the solver settings (KSP type, preconditioner,
 * tolerances) to the PETSc KSP object before each solve. Direct factorisations (LU, Cholesky) use MUMPS.
 */
#pragma once

#include <string>
#include <utility>
#include <vector>

#include <libmesh/linear_solver.h>
#include <libmesh/solver_configuration.h>

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
};

/** @brief Settings of the libMesh Newton solver of the nonlinear (law) systems. */
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
    double minimum_step_length = 1e-5;
    double linear_tolerance_multiplier = 1e-3;
    bool continue_after_max_iterations = true;
    bool continue_after_backtrack_failure = false;
    bool verbose = true;
};

/** @brief The solver settings of a model: the Python solver_settings dict, frozen when the model is loaded. */
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
 * @brief Applies PETSc KSP/PC options to the linear solver before each solve.
 *
 * Prints a one-time configuration summary table to the log.
 */
class PetscSolverConfiguration : public libMesh::SolverConfiguration
{
public:
    /** @param nonlinear The linear solver of the Newton iterations (for the report of the settings). */
    PetscSolverConfiguration(const SolverSettings &settings, libMesh::LinearSolver<libMesh::Number> *linear_solver,
                             bool nonlinear = false);

    /** @brief Applies KSP and PC settings to the attached linear solver. */
    void configure_solver() override;

private:
    libMesh::LinearSolver<libMesh::Number> *linear_solver_ = nullptr; ///< Non-owning.
    SolverSettings settings_;
    bool nonlinear_ = false;
};

/**
 * @file time_step.h
 * @brief Time-step state of the equation systems and the time scheme it gives.
 *
 * The DOF system solves M a + C v + K u = F at the end of the step, u, v, a the solution and its first and second
 * time derivatives. The scheme writes v and a as linear functions of the unknown u and of the state at the
 * start of the step (u_n, v_n, a_n):
 *
 *   v = c_v u + (hu_v u_n + hv_v v_n + ha_v a_n),   a = c_a u + (hu_a u_n + hv_a v_n + ha_a a_n)
 *
 * so the step is (K + c_v C + c_a M) u = F - C h_v - M h_a, and v, a follow from u after the solve.
 *  - First-order formulations (time_order 1), theta method: v = (u - u_n) / (theta dt) - (1 - theta) / theta v_n.
 *  - Second-order formulations (time_order 2), Newmark (beta, gamma):
 *      a = (u - u_n - dt v_n) / (beta dt^2) - (1 / (2 beta) - 1) a_n,   v = v_n + dt ((1 - gamma) a_n + gamma a).
 *    gamma = 1/2, beta = 1/4 is the unconditionally stable average acceleration scheme.
 *  - A static step (no transient effects, or time_order 0): v = a = 0.
 * The time order of the system is the highest of its element types.
 */
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <libmesh/libmesh_common.h>

/** @brief The step the caller asks for: end time and time scheme parameters (the Python time_step dict). */
struct TimeStepInput
{
    libMesh::Real time = 0.0;
    libMesh::Real theta = 1.0; ///< Theta method (first order).
    libMesh::Real beta = 0.25; ///< Newmark beta (second order).
    libMesh::Real gamma = 0.5; ///< Newmark gamma (second order).
    bool transient_effects = false;
};

/** @brief The time-step state of the equation systems. */
struct TimeStep
{
    libMesh::Real time = 0.0;
    libMesh::Real theta = 1.0;  ///< Theta method (first order).
    libMesh::Real beta = 0.25;  ///< Newmark beta (second order).
    libMesh::Real gamma = 0.5;  ///< Newmark gamma (second order).
    bool transient_effects = false;
    unsigned int time_step_idx = 0;
    libMesh::Real dt = 0.0;
    bool updated = true;
};

/** @brief What happened during the solve of one time step, filled by the systems; reported at its end. */
/** @brief Statistics of a MUMPS factorization (MUMPS INFOG / RINFOG), all ranks. */
struct FactorizationReport
{
    bool available = false;        ///< The preconditioner is a MUMPS factorization.
    int status = 0;                ///< INFOG(1): 0 ok, < 0 error, > 0 warning.
    int status_detail = 0;         ///< INFOG(2).
    double flops = 0.0;            ///< RINFOG(3): of the elimination.
    int memory_max_mb = 0;         ///< INFOG(21): memory used by the factorization, largest rank.
    int memory_total_mb = 0;       ///< INFOG(22): memory used by the factorization, all ranks.
    double matrix_entries = 0.0;   ///< Nonzeros of the factorized matrix.
    double factor_entries = 0.0;   ///< INFOG(29): entries of the factors.
    int null_pivots = -1;          ///< INFOG(28), when null pivot detection is on (ICNTL(24) = 1); else -1.
    int delayed_pivots = 0;        ///< INFOG(13): delayed pivots, a sign of numerical difficulties.
};

/** @brief One Newton iteration: the residual at the iterate and the step that led to it. */
struct NewtonIteration
{
    libMesh::Real residual = 0.0;      ///< l2 norm of the residual at the iterate.
    libMesh::Real step = -1.0;         ///< l2 norm of the step to the iterate; -1 for the initial guess.
    libMesh::Real relative_step = -1.0;
    libMesh::Real step_length = -1.0;  ///< Step taken / full Newton step (1 without line search); -1 if unknown.
};

/** @brief The largest input of a nonlinear law (e.g. |B| of a B-H curve) on the elements of one id. */
struct LawRange
{
    std::string law;
    std::string input;
    unsigned int id = 0;            ///< MATERIAL, CONSTANT or TYPE id, as the source of the law says.
    std::string source;             ///< "MATERIAL", "CONSTANT" or "TYPE".
    libMesh::Real max_input = 0.0;  ///< Largest element mean of |input|.
    libMesh::Real curve_max = 0.0;  ///< Largest x of the curve of the law (0 without curve).
    bool extrapolated = false;      ///< max_input is beyond the curve.
};

struct StepReport
{
    libMesh::Real time = 0.0;
    libMesh::Real dt = 0.0;
    unsigned int index = 0;
    bool transient = false;
    bool advanced = true;                   ///< A new time, or the same time solved again.
    bool element_matrices_computed = false; ///< Else reused from the previous steps.
    std::size_t elements = 0;               ///< Active elements (all ranks).
    bool system_matrix_updated = false;     ///< Else the factorization of the previous step was reused.
    bool nonlinear = false;
    unsigned int linear_iterations = 0;     ///< Of the linear solve, or of all the Newton iterations.
    libMesh::Real linear_residual = 0.0;
    std::string linear_reason;
    std::size_t unconstrained_dofs = 0;     ///< Unknowns of the reduced system (all ranks).
    std::size_t total_dofs = 0;
    unsigned int newton_iterations = 0;
    std::string newton_result;
    libMesh::Real true_residual = -1.0;     ///< ||b - A x|| / ||b|| of the (reduced) linear system; -1 if not computed.
    FactorizationReport factorization;
    std::vector<NewtonIteration> newton_history;
    std::vector<LawRange> law_ranges;
    double seconds_element_matrices = 0.0;
    double seconds_solve = 0.0; ///< Assembly and solve.
    double seconds_nodal_fields = 0.0;
    double seconds_element_fields = 0.0;
    double seconds_output = 0.0;
};

/** @brief The coefficients of v and a in u and in the state at the start of the step (see the file comment). */
struct TimeScheme
{
    libMesh::Real c_v = 0.0, hu_v = 0.0, hv_v = 0.0, ha_v = 0.0;
    libMesh::Real c_a = 0.0, hu_a = 0.0, hv_a = 0.0, ha_a = 0.0;

    /** @brief True when the step is not static: C and M enter the system. */
    bool transient() const { return c_v != 0.0 || c_a != 0.0; }

    static TimeScheme make(const TimeStep &step, unsigned int time_order);
};

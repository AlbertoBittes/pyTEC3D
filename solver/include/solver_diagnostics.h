/**
 * @file solver_diagnostics.h
 * @brief Quality checks of a linear solve: the true residual and the statistics of a MUMPS factorization.
 *
 * Both are collective (every rank calls them) and run once per solve, outside the assembly loops.
 */
#pragma once

#include <petscksp.h>

#include <libmesh/libmesh_common.h>

#include "time_step.h"

namespace pytec3d::diagnostics
{
    /** @brief ||b - A x|| / ||b|| (||b - A x|| when b = 0). */
    libMesh::Real relative_true_residual(Mat A, Vec x, Vec b);

    /** @brief The statistics of the factorization of @p ksp, when its preconditioner is a MUMPS LU or Cholesky. */
    FactorizationReport factorization_report(KSP ksp);

    /** @brief Warns about a failed or suspicious factorization and a large true residual. */
    void check_linear_solve(const FactorizationReport &factorization, libMesh::Real true_residual);
} // namespace pytec3d::diagnostics

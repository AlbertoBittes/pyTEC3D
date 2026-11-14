/**
 * @file solver_diagnostics.cpp
 * @brief True residual of a linear solve and MUMPS factorization statistics.
 */
#include "la.h"
#include "solver_diagnostics.h"

#include "logging.h"


#include <petscmat.h>

#include <cstring>

namespace pytec3d::diagnostics
{
    namespace
    {
        /** @brief True residuals above this are reported: a direct solve of a well-posed system is far below. */
        constexpr pytec3d::la::Real kTrueResidualWarning = 1e-6;

        /** @brief MUMPS gives large counts as minus the count in millions. */
        double mumps_count(PetscInt value)
        {
            return value < 0 ? -static_cast<double>(value) * 1e6 : static_cast<double>(value);
        }
    } // namespace


    FactorizationReport factorization_report(KSP ksp)
    {
        FactorizationReport report;
#if defined(PETSC_HAVE_MUMPS)
        if (!ksp)
            return report;
        PC pc = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, KSPGetPC(ksp, &pc));
        PetscBool factor = PETSC_FALSE;
        PetscCallAbort(PETSC_COMM_WORLD, PetscObjectTypeCompareAny(reinterpret_cast<PetscObject>(pc), &factor,
                                                                      PCLU, PCCHOLESKY, ""));
        if (!factor)
            return report;
        Mat F = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, PCFactorGetMatrix(pc, &F));
        MatSolverType type = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, MatFactorGetSolverType(F, &type));
        if (!type || std::strcmp(type, MATSOLVERMUMPS) != 0)
            return report;

        const auto infog = [F](PetscInt i)
        {
            PetscInt value = 0;
            PetscCallAbort(PETSC_COMM_WORLD, MatMumpsGetInfog(F, i, &value));
            return value;
        };
        report.available = true;
        report.status = static_cast<int>(infog(1));
        report.status_detail = static_cast<int>(infog(2));
        if (report.status < 0)
            return report;
        PetscReal flops = 0.0;
        PetscCallAbort(PETSC_COMM_WORLD, MatMumpsGetRinfog(F, 3, &flops));
        report.flops = flops;
        report.memory_max_mb = static_cast<int>(infog(21));
        report.memory_total_mb = static_cast<int>(infog(22));
        report.factor_entries = mumps_count(infog(29));
        report.delayed_pivots = static_cast<int>(infog(13));
        PetscInt null_pivot_detection = 0;
        PetscCallAbort(PETSC_COMM_WORLD, MatMumpsGetIcntl(F, 24, &null_pivot_detection));
        if (null_pivot_detection == 1)
            report.null_pivots = static_cast<int>(infog(28));

        Mat A = nullptr;
        PetscCallAbort(PETSC_COMM_WORLD, PCGetOperators(pc, &A, nullptr));
        MatInfo info;
        PetscCallAbort(PETSC_COMM_WORLD, MatGetInfo(A, MAT_GLOBAL_SUM, &info));
        report.matrix_entries = info.nz_used;
#else
        (void)ksp;
#endif
        return report;
    }

    void check_linear_solve(const FactorizationReport &factorization, pytec3d::la::Real true_residual)
    {
        if (factorization.available)
        {
            if (factorization.status < 0)
                logging::warn("The MUMPS factorization failed: INFOG(1) = {}, INFOG(2) = {} (see the MUMPS users' guide; "
                              "-9 or -13: not enough memory, -10: the matrix is singular).",
                              factorization.status, factorization.status_detail);
            else if (factorization.status > 0)
                logging::warn("MUMPS warning: INFOG(1) = {}, INFOG(2) = {}.", factorization.status,
                              factorization.status_detail);
            if (factorization.null_pivots > 0)
                logging::warn("MUMPS found {} null pivot(s): the matrix is singular (a missing constraint or a floating "
                              "region?).",
                              factorization.null_pivots);
        }
        if (true_residual > kTrueResidualWarning)
            logging::warn("Large true residual of the linear solve: ||b - Ax|| / ||b|| = {:.3e}; the system may be "
                          "ill-conditioned or singular.",
                          true_residual);
    }
} // namespace pytec3d::diagnostics

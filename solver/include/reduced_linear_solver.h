/**
 * @file reduced_linear_solver.h
 * @brief ReducedPetscLinearSolver: the PETSc linear solver of the pyTEC3D systems, on the unconstrained DOFs, with the
 *        constraints applied to the global (unconstrained) matrix and right-hand side.
 *
 * The constraints of the DofMap (its processed rows: each constrained DOF in terms of unconstrained ones, and their
 * values) write the DOFs as u = P w + g: w the unconstrained DOFs, P (rows: all the DOFs, columns: the unconstrained
 * ones) the identity on the unconstrained DOFs and the constraint rows on the others, g the values. A system A u = b
 * becomes (P^T A P) w = P^T (b - A g), u = P w + g; a Newton step (homogeneous constraints) delta = P w, with the
 * residual projected beforehand (project_residual()). No element matrix is constrained: the system matrix is the
 * global, unconstrained one, assembled once.
 *
 * The reduced matrix is kept between solves, so each solve does the least factorization work:
 *   - reused: the same matrix (reuse_preconditioner(true)): only the triangular solves;
 *   - numeric: new values on the same pattern (a Newton iteration, a new time step size or scheme): P^T A P is refilled
 *     in place (MatPtAP, MAT_REUSE_MATRIX), PETSc redoes only the numeric factorization;
 *   - full: other constraint rows (another P), another matrix or other solver settings (invalidate_pattern()): a new
 *     KSP, symbolic and numeric factorization.
 */
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <libmesh/dof_map.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/wrapped_petsc.h>

class ReducedPetscLinearSolver : public libMesh::PetscLinearSolver<libMesh::Number>
{
public:
    enum class Factorization
    {
        reused,
        numeric,
        full
    };

    /** @brief Factorizations done since the construction (all the solves of all the steps). */
    struct Counts
    {
        unsigned int full = 0;
        unsigned int numeric = 0;
        unsigned int reused = 0;
        double seconds_ptap = 0.0;          ///< P^T A P (the reduced matrix).
        double seconds_factorization = 0.0; ///< KSPSetUp: the preconditioner / the factorization (MUMPS).
        double seconds_linear_solve = 0.0;  ///< KSPSolve after the setup: the triangular solves or the iterations.
    };

    explicit ReducedPetscLinearSolver(const libMesh::Parallel::Communicator &comm);

    using libMesh::PetscLinearSolver<libMesh::Number>::solve;

    void clear() override;
    void init(const char *name = nullptr) override;

    /**
     * @brief The constraints of @p dof_map (after DofMap::process_constraints()): P and g. @p rows_changed: other
     * constrained DOFs or rows than at the last call (a new P, a full factorization); otherwise only the values.
     */
    void set_constraints(libMesh::DofMap &dof_map, bool rows_changed);

    /** @brief Newton steps: homogeneous constraints (delta = P w, the right-hand side already projected). */
    void set_homogeneous(bool homogeneous) { homogeneous_ = homogeneous; }

    /** @brief r <- P^T r at the unconstrained DOFs, 0 at the constrained ones: the residual of the reduced system. */
    void project_residual(libMesh::NumericVector<libMesh::Number> &r);

    std::pair<unsigned int, libMesh::Real> solve(libMesh::SparseMatrix<libMesh::Number> &matrix,
                                                 libMesh::SparseMatrix<libMesh::Number> &precond,
                                                 libMesh::NumericVector<libMesh::Number> &solution,
                                                 libMesh::NumericVector<libMesh::Number> &rhs,
                                                 const std::optional<double> tol = std::nullopt,
                                                 const std::optional<unsigned int> m_its = std::nullopt) override;

    /** @brief The next factorization is a full one (other settings). */
    void invalidate_pattern() { pattern_valid_ = false; }

    Factorization last_factorization() const { return last_; }
    const Counts &counts() const { return counts_; }

    /** @brief The reduced system of the last solve (for the true residual). */
    Mat reduced_matrix() const { return submat_; }
    Vec reduced_rhs() const { return subrhs_; }
    Vec reduced_solution() const { return subsolution_; }
    /** @brief Unconstrained DOFs of this rank. */
    PetscInt n_free_local() const { return n_free_local_; }

private:
    /** @brief A new KSP (the reduced matrix of a full factorization). */
    void reset_ksp();

    libMesh::WrappedPetsc<Mat> P_;            ///< Constraint prolongation: all the DOFs x the unconstrained ones.
    libMesh::WrappedPetsc<Vec> g_;            ///< The constraint values (0 at the unconstrained DOFs).
    libMesh::WrappedPetsc<VecScatter> free_;  ///< The unconstrained DOFs of a global vector -> a reduced vector.
    libMesh::WrappedPetsc<Mat> submat_;       ///< P^T A P.
    libMesh::WrappedPetsc<Vec> subrhs_;
    libMesh::WrappedPetsc<Vec> subsolution_;
    Mat parent_ = nullptr;                    ///< The matrix submat_ was computed from.
    PetscObjectState parent_state_ = -1;      ///< Its nonzero state then.
    PetscInt n_free_local_ = 0;
    bool homogeneous_ = false;
    bool pattern_valid_ = false;
    std::string prefix_;
    bool has_prefix_ = false;
    Factorization last_ = Factorization::full;
    Counts counts_;
};

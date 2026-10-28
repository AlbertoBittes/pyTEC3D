/**
 * @file reduced_linear_solver.h
 * @brief ReducedPetscLinearSolver: the PETSc linear solver of the pyTEC3D systems, restricted to the unconstrained
 *        DOFs, keeping between solves what the next factorization can reuse.
 */
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <libmesh/petsc_linear_solver.h>
#include <libmesh/wrapped_petsc.h>

/**
 * @brief The linear solver of the linear system and of the Newton iterations, on the unconstrained DOFs
 *        (restrict_solve_to, SUBSET_ZERO).
 *
 * The subset, the scatter, the reduced vectors and the reduced matrix are kept between solves, so each solve does
 * the least factorization work:
 *   - reused: the same matrix (reuse_preconditioner(true)): only the triangular solves;
 *   - numeric: new values on the same sparsity pattern (a Newton iteration, a new time step size or scheme): the
 *     reduced matrix is refilled in place (MAT_REUSE_MATRIX), PETSc sees the same matrix with the same nonzero state
 *     and redoes only the numeric factorization (ordering and analysis kept);
 *   - full: a new subset (constrained DOFs added or removed), other constraint rows (invalidate_pattern()) or
 *     another number of nonzeros of the matrix: a new KSP, symbolic and numeric factorization.
 * The restricted solve of libMesh extracts a new reduced matrix at every solve, which forces a full factorization.
 * Without a subset the solve is on all the local DOFs. The preconditioning matrix must be the system matrix.
 */
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
    };

    explicit ReducedPetscLinearSolver(const libMesh::Parallel::Communicator &comm);

    using libMesh::PetscLinearSolver<libMesh::Number>::solve;

    void clear() override;
    void init(const char *name = nullptr) override;
    void init_systems(const libMesh::System &system) override;

    /** @brief The DOFs of the solve (local, no duplicates); only SUBSET_ZERO. A new subset: full factorization. */
    void restrict_solve_to(const std::vector<unsigned int> *const dofs,
                           const libMesh::SubsetSolveMode subset_solve_mode = libMesh::SUBSET_ZERO) override;

    std::pair<unsigned int, libMesh::Real> solve(libMesh::SparseMatrix<libMesh::Number> &matrix,
                                                 libMesh::SparseMatrix<libMesh::Number> &precond,
                                                 libMesh::NumericVector<libMesh::Number> &solution,
                                                 libMesh::NumericVector<libMesh::Number> &rhs,
                                                 const std::optional<double> tol = std::nullopt,
                                                 const std::optional<unsigned int> m_its = std::nullopt) override;

    /** @brief The constraint rows changed (possibly another sparsity pattern): the next factorization is a full one. */
    void invalidate_pattern() { pattern_valid_ = false; }

    Factorization last_factorization() const { return last_; }
    const Counts &counts() const { return counts_; }

    /** @brief The reduced system of the last solve (for the true residual). */
    Mat reduced_matrix() const { return submat_; }
    Vec reduced_rhs() const { return subrhs_; }
    Vec reduced_solution() const { return subsolution_; }

private:
    /** @brief A new KSP and a new reduced system on the subset of @p parent. */
    void build(Mat parent, Vec rhs);
    void clear_reduced_system();

    std::vector<unsigned int> dofs_;
    bool has_subset_ = false;
    std::string prefix_;
    bool has_prefix_ = false;
    const libMesh::System *system_ = nullptr;

    libMesh::WrappedPetsc<IS> is_;
    libMesh::WrappedPetsc<VecScatter> scatter_;
    libMesh::WrappedPetsc<Mat> submat_;
    libMesh::WrappedPetsc<Vec> subrhs_;
    libMesh::WrappedPetsc<Vec> subsolution_;
    Mat parent_ = nullptr;          ///< Matrix the reduced matrix was extracted from.
    PetscLogDouble parent_nonzeros_ = -1.0; ///< Its nonzeros (all ranks) at the last full factorization.
    bool pattern_valid_ = false;

    Factorization last_ = Factorization::full;
    Counts counts_;
};

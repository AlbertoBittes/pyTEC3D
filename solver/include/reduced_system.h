/**
 * @file reduced_system.h
 * @brief ReducedSystem: the constrained linear systems of the DOF system, solved on the free DOFs.
 *
 * The constraints (constraints.h) write the unknowns as u = P w + g: w the free (unconstrained) DOFs, P the constraint
 * prolongation (a HypreParMatrix: identity on the free DOFs, the rows of the constrained ones resolved to free DOFs, the
 * chains of masters followed), g the resolved values. A system K u = f becomes (P^T K P) w = P^T (f - K g) (hypre RAP),
 * a Newton step on u the step on w with P^T J P and the residual P^T r: the constraints hold exactly.
 *
 * The reduced matrix goes to the PETSc linear solver (KSP: the solver settings, MUMPS for the direct factorizations).
 * Each solve does the least factorization work:
 *   - reused: the same matrix: only the triangular solves;
 *   - numeric: new values with the same pattern (a Newton iteration, a new time step size or scheme): the PETSc matrix is
 *     updated in place, MUMPS redoes only the numeric factorization;
 *   - full: another pattern (other constraint rows or DOFs), or new linear solver settings: a new KSP, symbolic and
 *     numeric factorization.
 */
#pragma once

#include <memory>
#include <vector>

#include <mfem.hpp>
#include <petscksp.h>

#include "constraints.h"
#include "dof_space.h"
#include "solver_configuration.h"

class ReducedSystem
{
public:
    enum class Factorization
    {
        reused,
        numeric,
        full
    };
    struct Counts
    {
        unsigned int full = 0;
        unsigned int numeric = 0;
        unsigned int reused = 0;
    };
    /** @brief What new constraints changed. */
    enum class Change
    {
        none,   ///< The same rows and values.
        values, ///< Only the values (g): the same reduced matrix.
        rows,   ///< The same constrained DOFs, other rows: another P.
        dofs    ///< Other constrained DOFs: another P and size.
    };

    explicit ReducedSystem(const pytec3d::fem::DofSpace &dofs);
    ~ReducedSystem();
    ReducedSystem(const ReducedSystem &) = delete;
    ReducedSystem &operator=(const ReducedSystem &) = delete;

    /** @brief Takes new constraints (the same on every rank); builds P and g. */
    Change set_constraints(const pytec3d::ConstraintSet &constraints);
    bool has_constraints() const { return P_ != nullptr; }

    /** @brief Free DOFs, all ranks / this rank. */
    HYPRE_BigInt n_free() const { return n_free_global_; }
    int n_free_local() const { return static_cast<int>(free_partition_[1] - free_partition_[0]); }

    /** @brief u = P w + g. */
    void expand(const mfem::Vector &w, mfem::Vector &u) const;
    /** @brief u = P w (a step: homogeneous constraints). */
    void prolong(const mfem::Vector &w, mfem::Vector &u) const { P_->Mult(w, u); }
    /** @brief The free DOFs of @p u. */
    void restrict_to_free(const mfem::Vector &u, mfem::Vector &w) const;
    /** @brief w = P^T f. */
    void project(const mfem::Vector &f, mfem::Vector &w) const { P_->MultTranspose(f, w); }
    /** @brief The values g of the constrained DOFs, zero elsewhere (T-vector). */
    const mfem::Vector &offset() const { return g_; }

    /** @brief The reduced matrix P^T K P of @p K; @p changed: K has other values than at the last call. */
    void set_matrix(const mfem::HypreParMatrix &K, bool changed);

    /** @brief The settings of the next solves; true when they changed (a full factorization next). */
    bool set_settings(const SolverSettings &settings, bool nonlinear);

    /** @brief Solves the reduced system; @p solution the initial guess when the settings say so. */
    std::pair<unsigned int, double> solve(const mfem::Vector &rhs, mfem::Vector &solution, double rtol);

    /** @brief ||rhs - Khat x|| / ||rhs|| (all ranks); -1 for a zero rhs. */
    double relative_true_residual(const mfem::Vector &x, const mfem::Vector &rhs) const;

    KSP ksp() const { return ksp_; }
    KSPConvergedReason converged_reason() const;
    Factorization last_factorization() const { return last_; }
    const Counts &counts() const { return counts_; }

private:
    void clear_ksp();
    /** @brief The PETSc matrix of Khat_: in place when the pattern is that of the previous one. */
    void update_petsc_matrix();

    const pytec3d::fem::DofSpace &dofs_;
    pytec3d::ConstraintSet constraints_;
    bool has_set_ = false;
    std::unique_ptr<mfem::HypreParMatrix> P_;
    mfem::Vector g_;
    std::vector<int> free_local_; ///< The owned free DOFs: their positions in the T-vectors.
    HYPRE_BigInt free_partition_[2] = {0, 0};
    HYPRE_BigInt n_free_global_ = 0;

    const mfem::HypreParMatrix *K_ = nullptr; ///< The matrix Khat_ was computed from.
    std::unique_ptr<mfem::HypreParMatrix> Khat_; ///< P^T K P while it is moved into A_.
    bool khat_new_ = true;     ///< Khat_ has a new pattern: full factorization.
    bool khat_changed_ = true; ///< Khat_ has new values since the last solve.

    Mat A_ = nullptr; ///< Khat_ as a PETSc matrix.
    std::vector<PetscInt> csr_i_, csr_j_;
    KSP ksp_ = nullptr;
    PetscSolverConfiguration configuration_;
    Factorization last_ = Factorization::full;
    Counts counts_;
};

/**
 * @file reduced_linear_solver.cpp
 * @brief ReducedPetscLinearSolver: constraints as a prolongation P, the reduced system P^T A P (reduced_linear_solver.h).
 */
#include <chrono>
#include "reduced_linear_solver.h"

#include <algorithm>
#include <map>
#include <set>

#include <libmesh/petsc_matrix_base.h>
#include <libmesh/solver_configuration.h>
#include <libmesh/petsc_solver_exception.h>
#include <libmesh/petsc_vector.h>

#include "logging.h"

using libMesh::libmesh_ignore; // used by the LibmeshPetscCall macros

namespace
{
    using Base = libMesh::PetscLinearSolver<libMesh::Number>;
    using libMesh::dof_id_type;
} // namespace

ReducedPetscLinearSolver::ReducedPetscLinearSolver(const libMesh::Parallel::Communicator &comm) : Base(comm)
{
}

void ReducedPetscLinearSolver::clear()
{
    P_.reset_to_zero();
    g_.reset_to_zero();
    free_.reset_to_zero();
    submat_.reset_to_zero();
    subrhs_.reset_to_zero();
    subsolution_.reset_to_zero();
    parent_ = nullptr;
    parent_state_ = -1;
    pattern_valid_ = false;
    Base::clear();
}

void ReducedPetscLinearSolver::init(const char *name)
{
    // The prefix of the first init, given again to the KSP of each full factorization.
    if (!this->initialized())
    {
        has_prefix_ = name != nullptr;
        prefix_ = name ? std::string(name) : std::string();
    }
    Base::init(name);
}

void ReducedPetscLinearSolver::reset_ksp()
{
    // A new KSP: the operators of the previous one may have other sizes. P and g are kept.
    Base::clear();
    Base::init(has_prefix_ ? prefix_.c_str() : nullptr);
}

void ReducedPetscLinearSolver::set_constraints(libMesh::DofMap &dof_map, bool rows_changed)
{
    const auto &comm = this->comm();
    const unsigned int n_ranks = comm.size();
    const dof_id_type first = dof_map.first_dof(), end = dof_map.end_dof();
    const libMesh::DofConstraints &rows = dof_map.get_dof_constraints();
    const libMesh::DofConstraintValueMap &values = dof_map.get_primal_constraint_values();

    // The constrained DOFs of every rank, on every rank: the index of any unconstrained DOF in the reduced system.
    std::vector<dof_id_type> constrained;
    for (const auto &[dof, row] : rows)
    {
        libMesh::libmesh_ignore(row);
        if (dof >= first && dof < end)
            constrained.push_back(dof);
    }
    std::sort(constrained.begin(), constrained.end());
    std::vector<dof_id_type> all_constrained = constrained;
    comm.allgather(all_constrained, false);
    std::sort(all_constrained.begin(), all_constrained.end());

    std::vector<dof_id_type> rank_first(n_ranks + 1, 0);
    for (unsigned int r = 0; r < n_ranks; ++r)
        rank_first[r + 1] = dof_map.end_dof(r);
    const auto constrained_between = [&](dof_id_type from, dof_id_type to)
    {
        return static_cast<PetscInt>(std::lower_bound(all_constrained.begin(), all_constrained.end(), to) -
                                     std::lower_bound(all_constrained.begin(), all_constrained.end(), from));
    };
    std::vector<PetscInt> free_first(n_ranks + 1, 0);
    for (unsigned int r = 0; r < n_ranks; ++r)
        free_first[r + 1] = free_first[r] + static_cast<PetscInt>(rank_first[r + 1] - rank_first[r]) -
                            constrained_between(rank_first[r], rank_first[r + 1]);
    const PetscInt rank = static_cast<PetscInt>(comm.rank());
    n_free_local_ = free_first[static_cast<std::size_t>(rank) + 1] - free_first[static_cast<std::size_t>(rank)];
    const PetscInt n_free = free_first.back();
    const auto free_index = [&](dof_id_type dof)
    {
        if (std::binary_search(all_constrained.begin(), all_constrained.end(), dof))
            pytec3d::logging::error_and_abort("Constraints: DOF " + std::to_string(dof) +
                                              " is the master of a constraint row but is constrained itself.");
        const auto r = static_cast<std::size_t>(std::upper_bound(rank_first.begin(), rank_first.end(), dof) -
                                                rank_first.begin()) - 1;
        return free_first[r] + static_cast<PetscInt>(dof - rank_first[r]) - constrained_between(rank_first[r], dof);
    };

    const PetscInt n_rows = static_cast<PetscInt>(end - first);
    const PetscInt n_global = static_cast<PetscInt>(dof_map.n_dofs());
    if (rows_changed || !P_)
    {
        // P: one row per owned DOF; an unconstrained DOF its column, a constrained one its masters.
        const PetscInt col0 = free_first[static_cast<std::size_t>(rank)], col1 = free_first[static_cast<std::size_t>(rank) + 1];
        std::vector<PetscInt> d_nnz(static_cast<std::size_t>(n_rows), 0), o_nnz(static_cast<std::size_t>(n_rows), 0);
        std::vector<PetscInt> free_dofs;
        for (PetscInt i = 0; i < n_rows; ++i)
        {
            const dof_id_type dof = first + static_cast<dof_id_type>(i);
            const auto row = rows.find(dof);
            if (row == rows.end())
            {
                d_nnz[static_cast<std::size_t>(i)] = 1;
                free_dofs.push_back(static_cast<PetscInt>(dof));
                continue;
            }
            for (const auto &[master, c] : row->second)
            {
                libMesh::libmesh_ignore(c);
                const PetscInt col = free_index(master);
                (col >= col0 && col < col1 ? d_nnz : o_nnz)[static_cast<std::size_t>(i)]++;
            }
        }
        P_.reset_to_zero();
        LibmeshPetscCall(MatCreateAIJ(comm.get(), n_rows, n_free_local_, n_global, n_free, 0, d_nnz.data(), 0, o_nnz.data(),
                                      P_.get()));
        for (PetscInt i = 0; i < n_rows; ++i)
        {
            const dof_id_type dof = first + static_cast<dof_id_type>(i);
            const auto row = rows.find(dof);
            const PetscInt r = static_cast<PetscInt>(dof);
            if (row == rows.end())
            {
                const PetscInt col = free_index(dof);
                LibmeshPetscCall(MatSetValue(P_, r, col, 1.0, INSERT_VALUES));
                continue;
            }
            for (const auto &[master, c] : row->second)
                if (c != 0.0)
                    LibmeshPetscCall(MatSetValue(P_, r, free_index(master), c, INSERT_VALUES));
        }
        LibmeshPetscCall(MatAssemblyBegin(P_, MAT_FINAL_ASSEMBLY));
        LibmeshPetscCall(MatAssemblyEnd(P_, MAT_FINAL_ASSEMBLY));

        // The scatter of the unconstrained DOFs: global vector <-> reduced vector, in order.
        libMesh::WrappedPetsc<IS> from, to;
        LibmeshPetscCall(ISCreateGeneral(comm.get(), static_cast<PetscInt>(free_dofs.size()), free_dofs.data(),
                                         PETSC_COPY_VALUES, from.get()));
        LibmeshPetscCall(ISCreateStride(comm.get(), n_free_local_, col0, 1, to.get()));
        libMesh::WrappedPetsc<Vec> u, w;
        LibmeshPetscCall(MatCreateVecs(P_, w.get(), u.get()));
        free_.reset_to_zero();
        LibmeshPetscCall(VecScatterCreate(u, from, w, to, free_.get()));

        // Another P: another reduced system, a full factorization.
        submat_.reset_to_zero();
        subrhs_.reset_to_zero();
        subsolution_.reset_to_zero();
        parent_ = nullptr;
        pattern_valid_ = false;
    }

    // g: the values of the owned constrained DOFs.
    g_.reset_to_zero();
    LibmeshPetscCall(MatCreateVecs(P_, nullptr, g_.get()));
    LibmeshPetscCall(VecZeroEntries(g_));
    for (const dof_id_type dof : constrained)
    {
        const auto it = values.find(dof);
        if (it != values.end() && it->second != 0.0)
            LibmeshPetscCall(VecSetValue(g_, static_cast<PetscInt>(dof), it->second, INSERT_VALUES));
    }
    LibmeshPetscCall(VecAssemblyBegin(g_));
    LibmeshPetscCall(VecAssemblyEnd(g_));
}

void ReducedPetscLinearSolver::project_residual(libMesh::NumericVector<libMesh::Number> &r_in)
{
    auto *r = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&r_in);
    r->close();
    libMesh::WrappedPetsc<Vec> w;
    LibmeshPetscCall(MatCreateVecs(P_, w.get(), nullptr));
    LibmeshPetscCall(MatMultTranspose(P_, r->vec(), w));
    LibmeshPetscCall(VecZeroEntries(r->vec()));
    libMesh::VecScatterBeginEnd(this->comm(), free_, w, r->vec(), INSERT_VALUES, SCATTER_REVERSE);
}

std::pair<unsigned int, libMesh::Real> ReducedPetscLinearSolver::solve(libMesh::SparseMatrix<libMesh::Number> &matrix_in,
                                                                       libMesh::SparseMatrix<libMesh::Number> &precond_in,
                                                                       libMesh::NumericVector<libMesh::Number> &solution_in,
                                                                       libMesh::NumericVector<libMesh::Number> &rhs_in,
                                                                       const std::optional<double> tol,
                                                                       const std::optional<unsigned int> m_its)
{
    if (&precond_in != &matrix_in)
        pytec3d::logging::error_and_abort("ReducedPetscLinearSolver::solve: the preconditioning matrix must be the system matrix.");
    if (!P_)
        pytec3d::logging::error_and_abort("ReducedPetscLinearSolver::solve: no constraints (set_constraints()).");
    auto *matrix = libMesh::cast_ptr<libMesh::PetscMatrixBase<libMesh::Number> *>(&matrix_in);
    auto *solution = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&solution_in);
    auto *rhs = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&rhs_in);
    this->init();
    matrix->close();
    solution->close();
    rhs->close();

    // The level, the same on all the ranks: full without a reduced system of this P and matrix (and its pattern);
    // else reused (the same matrix) or numeric.
    const Mat parent = matrix->mat();
    PetscObjectState state = 0;
    LibmeshPetscCall(MatGetNonzeroState(parent, &state));
    bool full = !submat_ || !pattern_valid_ || parent != parent_ || state != parent_state_;
    this->comm().max(full);
    const Factorization level = full ? Factorization::full
                                : this->same_preconditioner ? Factorization::reused
                                                            : Factorization::numeric;
    const auto t_ptap = std::chrono::steady_clock::now();
    if (level == Factorization::full)
    {
        reset_ksp();
        submat_.reset_to_zero();
        LibmeshPetscCall(MatPtAP(parent, P_, MAT_INITIAL_MATRIX, PETSC_DETERMINE, submat_.get()));
        subrhs_.reset_to_zero();
        subsolution_.reset_to_zero();
        LibmeshPetscCall(MatCreateVecs(submat_, subsolution_.get(), subrhs_.get()));
        parent_ = parent;
        parent_state_ = state;
        pattern_valid_ = true;
    }
    else if (level == Factorization::numeric)
        // The values of the same nonzeros: PETSc keeps the ordering and the analysis of the factorization.
        LibmeshPetscCall(MatPtAP(parent, P_, MAT_REUSE_MATRIX, PETSC_DETERMINE, submat_.get()));
    counts_.seconds_ptap += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ptap).count();

    // The reduced right-hand side: P^T (b - A g); a Newton step: b, projected already (project_residual()).
    if (homogeneous_)
        libMesh::VecScatterBeginEnd(this->comm(), free_, rhs->vec(), subrhs_, INSERT_VALUES, SCATTER_FORWARD);
    else
    {
        libMesh::WrappedPetsc<Vec> t;
        LibmeshPetscCall(VecDuplicate(rhs->vec(), t.get()));
        LibmeshPetscCall(MatMult(parent, g_, t));
        LibmeshPetscCall(VecAYPX(t, -1.0, rhs->vec()));
        LibmeshPetscCall(MatMultTranspose(P_, t, subrhs_));
    }
    libMesh::VecScatterBeginEnd(this->comm(), free_, solution->vec(), subsolution_, INSERT_VALUES, SCATTER_FORWARD);

    const KSP ksp = this->ksp();
    LibmeshPetscCall(KSPSetOperators(ksp, submat_, submat_));
    LibmeshPetscCall(KSPSetReusePreconditioner(ksp, level == Factorization::reused ? PETSC_TRUE : PETSC_FALSE));
    // As PetscLinearSolver::solve: the tolerances of the caller or of the settings, the command line, then the
    // pyTEC3D settings.
    const double rel_tol = this->get_real_solver_setting("rel_tol", tol);
    const double abs_tol = this->get_real_solver_setting("abs_tol", std::nullopt, static_cast<double>(PETSC_DEFAULT));
    const double max_its = this->get_int_solver_setting("max_its", m_its);
    LibmeshPetscCall(KSPSetTolerances(ksp, rel_tol, abs_tol, PETSC_DEFAULT, static_cast<PetscInt>(max_its)));
    LibmeshPetscCall(KSPSetFromOptions(ksp));
    if (this->_solver_configuration)
        this->_solver_configuration->configure_solver();
    // the setup apart (the factorization of a direct solver), then the solve
    const auto t_setup = std::chrono::steady_clock::now();
    LibmeshPetscCall(KSPSetUp(ksp));
    const auto t_solve = std::chrono::steady_clock::now();
    LibmeshPetscCall(KSPSolve(ksp, subrhs_, subsolution_));
    counts_.seconds_factorization += std::chrono::duration<double>(t_solve - t_setup).count();
    counts_.seconds_linear_solve += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_solve).count();

    PetscInt its = 0;
    PetscReal residual = 0.0;
    LibmeshPetscCall(KSPGetIterationNumber(ksp, &its));
    LibmeshPetscCall(KSPGetResidualNorm(ksp, &residual));

    // u = P w (+ g).
    LibmeshPetscCall(MatMult(P_, subsolution_, solution->vec()));
    if (!homogeneous_)
        LibmeshPetscCall(VecAXPY(solution->vec(), 1.0, g_));

    last_ = level;
    switch (level)
    {
    case Factorization::full:
        ++counts_.full;
        break;
    case Factorization::numeric:
        ++counts_.numeric;
        break;
    case Factorization::reused:
        ++counts_.reused;
        break;
    }
    return {static_cast<unsigned int>(its), static_cast<libMesh::Real>(residual)};
}

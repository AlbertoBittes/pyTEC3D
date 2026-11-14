/**
 * @file reduced_system.cpp
 * @brief ReducedSystem: constraint prolongation, reduced matrix (hypre RAP) and PETSc linear solver.
 */
#include "reduced_system.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

#include "logging.h"

namespace
{
    void check(PetscErrorCode ierr, const char *what)
    {
        if (ierr)
            pytec3d::logging::error_and_abort(std::string("PETSc error in ") + what + ".");
    }

    /** @brief A row resolved to free DOFs: u = sum c_m u_m + value, the masters free. */
    struct Resolved
    {
        std::map<PetscInt, double> masters;
        double value = 0.0;
    };
} // namespace

ReducedSystem::ReducedSystem(const pytec3d::fem::DofSpace &dofs) : dofs_(dofs)
{
}

ReducedSystem::~ReducedSystem()
{
    clear_ksp();
    if (A_)
        MatDestroy(&A_);
}

void ReducedSystem::clear_ksp()
{
    if (ksp_)
        KSPDestroy(&ksp_);
}

ReducedSystem::Change ReducedSystem::set_constraints(const pytec3d::ConstraintSet &constraints)
{
    Change change = Change::dofs;
    if (has_set_)
    {
        bool same_dofs = constraints.size() == constraints_.size();
        bool same_rows = same_dofs, same_values = same_dofs;
        if (same_dofs)
            for (auto a = constraints.cbegin(), b = constraints_.cbegin(); a != constraints.cend(); ++a, ++b)
            {
                if (a->first != b->first)
                {
                    same_dofs = same_rows = same_values = false;
                    break;
                }
                same_rows = same_rows && a->second.same_row(b->second);
                same_values = same_values && a->second.value == b->second.value;
            }
        change = !same_dofs ? Change::dofs : !same_rows ? Change::rows : !same_values ? Change::values : Change::none;
    }
    if (change == Change::none)
        return change;
    if (has_set_ && change == Change::dofs)
        pytec3d::logging::info("Boundary conditions update required: constrained DOF set changed");
    else if (change == Change::rows)
        pytec3d::logging::info("Boundary conditions update required: constraint rows changed");
    constraints_ = constraints;
    has_set_ = true;

    const MPI_Comm comm = dofs_.mesh().comm();
    const int rank = dofs_.mesh().rank();
    const std::vector<PetscInt> &owned = dofs_.owned_by_rank();
    const std::size_t n_ranks = owned.size();
    std::vector<PetscInt> first(n_ranks + 1, 0);
    for (std::size_t r = 0; r < n_ranks; ++r)
        first[r + 1] = first[r] + owned[r];

    // The free DOFs of each rank are numbered in order, after those of the lower ranks: every rank knows all the
    // constrained DOFs, hence the free index of any DOF.
    std::vector<PetscInt> keys;
    keys.reserve(constraints.size());
    for (const auto &[dof, row] : constraints)
        keys.push_back(dof);
    const auto constrained_between = [&](PetscInt from, PetscInt to)
    {
        return static_cast<PetscInt>(std::lower_bound(keys.begin(), keys.end(), to) - std::lower_bound(keys.begin(), keys.end(), from));
    };
    std::vector<PetscInt> free_first(n_ranks + 1, 0);
    for (std::size_t r = 0; r < n_ranks; ++r)
        free_first[r + 1] = free_first[r] + owned[r] - constrained_between(first[r], first[r + 1]);
    n_free_global_ = static_cast<HYPRE_BigInt>(free_first.back());
    free_partition_[0] = static_cast<HYPRE_BigInt>(free_first[static_cast<std::size_t>(rank)]);
    free_partition_[1] = static_cast<HYPRE_BigInt>(free_first[static_cast<std::size_t>(rank) + 1]);
    const auto free_index = [&](PetscInt g)
    {
        const auto r = static_cast<std::size_t>(std::upper_bound(first.begin(), first.end(), g) - first.begin()) - 1;
        return free_first[r] + (g - first[r]) - constrained_between(first[r], g);
    };

    // The rows resolved to free DOFs, chains of constrained masters followed.
    std::unordered_map<PetscInt, Resolved> resolved;
    std::set<PetscInt> visiting;
    std::function<const Resolved &(PetscInt)> resolve = [&](PetscInt dof) -> const Resolved &
    {
        if (const auto done = resolved.find(dof); done != resolved.end())
            return done->second;
        if (!visiting.insert(dof).second)
            pytec3d::logging::error_and_abort("Constraints: the constraint rows of DOF " + std::to_string(dof) +
                                              " form a cycle (couplings constraining each other).");
        Resolved out;
        const pytec3d::ConstraintRow &row = constraints.at(dof);
        out.value = row.value;
        for (const auto &[master, coefficient] : row.masters)
        {
            if (!constraints.count(master))
            {
                out.masters[master] += coefficient;
                continue;
            }
            const Resolved &m = resolve(master);
            out.value += coefficient * m.value;
            for (const auto &[free, c] : m.masters)
                out.masters[free] += coefficient * c;
        }
        visiting.erase(dof);
        return resolved.emplace(dof, std::move(out)).first->second;
    };

    // P: rows the owned unknowns, columns the free DOFs; g: the values.
    const PetscInt row0 = dofs_.first_owned();
    const int n_rows = static_cast<int>(dofs_.n_owned());
    std::vector<int> I(static_cast<std::size_t>(n_rows) + 1, 0);
    std::vector<HYPRE_BigInt> J;
    std::vector<mfem::real_t> data;
    g_.SetSize(n_rows);
    g_ = 0.0;
    free_local_.clear();
    for (int i = 0; i < n_rows; ++i)
    {
        const PetscInt g = row0 + i;
        if (!constraints.count(g))
        {
            J.push_back(static_cast<HYPRE_BigInt>(free_index(g)));
            data.push_back(1.0);
            free_local_.push_back(i);
        }
        else
        {
            const Resolved &r = resolve(g);
            for (const auto &[master, c] : r.masters)
                if (c != 0.0)
                {
                    J.push_back(static_cast<HYPRE_BigInt>(free_index(master)));
                    data.push_back(c);
                }
            g_(i) = r.value;
        }
        I[static_cast<std::size_t>(i) + 1] = static_cast<int>(J.size());
    }
    P_ = std::make_unique<mfem::HypreParMatrix>(comm, n_rows, static_cast<HYPRE_BigInt>(dofs_.n_global()), n_free_global_,
                                                I.data(), J.data(), data.data(), dofs_.true_partition(), free_partition_);

    // Another P: another reduced matrix and pattern, a full factorization.
    if (change == Change::dofs || change == Change::rows)
    {
        if (A_)
            MatDestroy(&A_);
        K_ = nullptr;
        khat_new_ = true;
        khat_changed_ = true;
    }
    return change;
}

void ReducedSystem::expand(const mfem::Vector &w, mfem::Vector &u) const
{
    P_->Mult(w, u);
    u += g_;
}

void ReducedSystem::restrict_to_free(const mfem::Vector &u, mfem::Vector &w) const
{
    w.SetSize(n_free_local());
    for (std::size_t k = 0; k < free_local_.size(); ++k)
        w(static_cast<int>(k)) = u(free_local_[k]);
}

void ReducedSystem::set_matrix(const mfem::HypreParMatrix &K, bool changed)
{
    if (A_ && &K == K_ && !changed && !khat_new_)
        return;
    // P^T K P, moved into the PETSc matrix of the solver at once: only one copy of the reduced matrix is kept.
    Khat_.reset(mfem::RAP(&K, P_.get()));
    K_ = &K;
    khat_changed_ = true;
    update_petsc_matrix();
    Khat_.reset();
}

void ReducedSystem::update_petsc_matrix()
{
    // Local rows with global columns, as PETSc's MPIAIJ takes them.
    mfem::SparseMatrix merged;
    Khat_->MergeDiagAndOffd(merged);
    merged.SortColumnIndices();
    const int m = merged.Height();
    const int *I = merged.GetI();
    const int *J = merged.GetJ();
    const double *values = merged.GetData();
    const bool same = A_ && !khat_new_ && static_cast<int>(csr_i_.size()) == m + 1 && std::equal(I, I + m + 1, csr_i_.begin()) &&
                      static_cast<int>(csr_j_.size()) == I[m] && std::equal(J, J + I[m], csr_j_.begin());
    if (same)
    {
        check(MatUpdateMPIAIJWithArray(A_, values), "MatUpdateMPIAIJWithArray");
        return;
    }
    if (A_)
        MatDestroy(&A_);
    csr_i_.assign(I, I + m + 1);
    csr_j_.assign(J, J + I[m]);
    check(MatCreateMPIAIJWithArrays(dofs_.mesh().comm(), m, m, static_cast<PetscInt>(n_free_global_),
                                    static_cast<PetscInt>(n_free_global_), csr_i_.data(), csr_j_.data(), values, &A_),
          "MatCreateMPIAIJWithArrays");
    khat_new_ = true;
}

bool ReducedSystem::set_settings(const SolverSettings &settings, bool nonlinear)
{
    if (!configuration_.set_settings(settings, nonlinear))
        return false;
    // Another KSP or PC type: a full factorization.
    pytec3d::logging::info("New linear solver settings: updating system matrix");
    clear_ksp();
    khat_new_ = true;
    khat_changed_ = true;
    return true;
}

std::pair<unsigned int, double> ReducedSystem::solve(const mfem::Vector &rhs, mfem::Vector &solution, double rtol)
{
    if (!A_)
        pytec3d::logging::error_and_abort("ReducedSystem::solve: no matrix.");
    if (n_free_global_ == 0)
        pytec3d::logging::error_and_abort("ReducedSystem::solve: every DOF is constrained.");

    if (!ksp_ || khat_new_)
    {
        clear_ksp();
        check(KSPCreate(dofs_.mesh().comm(), &ksp_), "KSPCreate");
        last_ = Factorization::full;
        ++counts_.full;
    }
    else if (khat_changed_)
    {
        last_ = Factorization::numeric;
        ++counts_.numeric;
    }
    else
    {
        last_ = Factorization::reused;
        ++counts_.reused;
    }
    check(KSPSetOperators(ksp_, A_, A_), "KSPSetOperators");
    check(KSPSetReusePreconditioner(ksp_, last_ == Factorization::reused ? PETSC_TRUE : PETSC_FALSE),
          "KSPSetReusePreconditioner");
    configuration_.configure(ksp_, rtol);

    // The mfem vectors as PETSc vectors, without copies.
    const MPI_Comm comm = dofs_.mesh().comm();
    solution.SetSize(n_free_local());
    Vec b = nullptr, x = nullptr;
    check(VecCreateMPIWithArray(comm, 1, n_free_local(), static_cast<PetscInt>(n_free_global_), rhs.GetData(), &b),
          "VecCreateMPIWithArray");
    check(VecCreateMPIWithArray(comm, 1, n_free_local(), static_cast<PetscInt>(n_free_global_), solution.GetData(), &x),
          "VecCreateMPIWithArray");
    check(KSPSolve(ksp_, b, x), "KSPSolve");
    VecDestroy(&b);
    VecDestroy(&x);
    khat_new_ = false;
    khat_changed_ = false;

    PetscInt its = 0;
    PetscReal norm = 0.0;
    check(KSPGetIterationNumber(ksp_, &its), "KSPGetIterationNumber");
    check(KSPGetResidualNorm(ksp_, &norm), "KSPGetResidualNorm");
    return {static_cast<unsigned int>(its), static_cast<double>(norm)};
}

double ReducedSystem::relative_true_residual(const mfem::Vector &x, const mfem::Vector &rhs) const
{
    const MPI_Comm comm = dofs_.mesh().comm();
    mfem::Vector r(rhs.Size());
    Vec px = nullptr, pr = nullptr;
    check(VecCreateMPIWithArray(comm, 1, n_free_local(), static_cast<PetscInt>(n_free_global_), x.GetData(), &px),
          "VecCreateMPIWithArray");
    check(VecCreateMPIWithArray(comm, 1, n_free_local(), static_cast<PetscInt>(n_free_global_), r.GetData(), &pr),
          "VecCreateMPIWithArray");
    check(MatMult(A_, px, pr), "MatMult");
    VecDestroy(&px);
    VecDestroy(&pr);
    r.Neg();
    r += rhs;
    const double rn = mfem::InnerProduct(comm, r, r), bn = mfem::InnerProduct(comm, rhs, rhs);
    return bn > 0.0 ? std::sqrt(rn / bn) : -1.0;
}

KSPConvergedReason ReducedSystem::converged_reason() const
{
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    if (ksp_)
        check(KSPGetConvergedReason(ksp_, &reason), "KSPGetConvergedReason");
    return reason;
}

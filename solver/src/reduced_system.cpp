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
#include "report_tables.h"

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
    u.SetSize(P_->Height());
    u = 0.0; // hypre forms beta * y with beta = 0: no garbage of a new allocation may survive
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
    if ((A_ || mfem_A_) && &K == K_ && !changed && !khat_new_)
        return;
    // P^T K P, moved into the PETSc matrix of the solver at once: only one copy of the reduced matrix is kept.
    Khat_.reset(mfem::RAP(&K, P_.get()));
    K_ = &K;
    khat_changed_ = true;
    if (mfem_backend())
    {
        // The pattern of the previous one: a numeric factorization (MUMPS reuses its analysis).
        const bool same = mfem_A_ && !khat_new_ && mfem_A_->NNZ() == Khat_->NNZ() &&
                          mfem_A_->Height() == Khat_->Height();
        if (!same)
            khat_new_ = true;
        mfem_A_ = std::move(Khat_);
        return;
    }
    mfem_A_.reset();
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
    mfem_solver_.reset();
    mfem_prec_.reset();
    khat_new_ = true;
    khat_changed_ = true;
    return true;
}

std::pair<unsigned int, double> ReducedSystem::solve(const mfem::Vector &rhs, mfem::Vector &solution, double rtol)
{
    if (mfem_backend())
        return solve_mfem(rhs, solution, rtol);
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
    r = 0.0;
    if (mfem_backend() && mfem_A_)
    {
        // A singular (ungauged) system: the residual of the consistent right-hand side it solved.
        const mfem::Vector &b = projected_ && consistent_rhs_.Size() == rhs.Size() ? consistent_rhs_ : rhs;
        mfem_A_->Mult(x, r);
        r.Neg();
        r += b;
        const double rn = mfem::InnerProduct(comm, r, r), bn = mfem::InnerProduct(comm, b, b);
        return bn > 0.0 ? std::sqrt(rn / bn) : -1.0;
    }
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
    if (mfem_backend())
        return mfem_converged_ ? KSP_CONVERGED_RTOL : KSP_DIVERGED_ITS;
    if (ksp_)
        check(KSPGetConvergedReason(ksp_, &reason), "KSPGetConvergedReason");
    return reason;
}

// ---- MFEM backend ---------------------------------------------------------------------------------------------------

namespace
{
    /** @brief A Krylov solver of MFEM by its PETSc name. */
    std::unique_ptr<mfem::IterativeSolver> krylov(const std::string &type, MPI_Comm comm)
    {
        if (type == "cg")
            return std::make_unique<mfem::CGSolver>(comm);
        if (type == "gmres")
            return std::make_unique<mfem::GMRESSolver>(comm);
        if (type == "fgmres")
            return std::make_unique<mfem::FGMRESSolver>(comm);
        if (type == "minres")
            return std::make_unique<mfem::MINRESSolver>(comm);
        if (type == "bcgs" || type == "bicgstab")
            return std::make_unique<mfem::BiCGSTABSolver>(comm);
        pytec3d::logging::error_and_abort("Linear solver (backend mfem): unknown ksp_type '" + type +
                                          "' (preonly, cg, gmres, fgmres, minres, bcgs).");
        return nullptr;
    }
} // namespace

std::unique_ptr<mfem::Solver> ReducedSystem::make_ams()
{
    // The one first-order Nedelec variable that holds every free DOF.
    int edge_var = -1;
    for (unsigned int v = 0; v < dofs_.n_vars(); ++v)
        if (dofs_.n_global(v) > 0)
        {
            if (!dofs_.variable(v).edge() || dofs_.variable(v).order != 1 || edge_var >= 0)
                pytec3d::logging::error_and_abort("Linear solver (backend mfem): pc_type 'ams' needs the DOFs of one "
                                                  "first-order Nedelec field only (e.g. AVA3D_edge without conductors).");
            edge_var = static_cast<int>(v);
        }
    if (edge_var < 0)
        pytec3d::logging::error_and_abort("Linear solver (backend mfem): pc_type 'ams' without a Nedelec field.");
    const unsigned int v = static_cast<unsigned int>(edge_var);
    auto &nd = const_cast<mfem::ParFiniteElementSpace &>(dofs_.space(v));
    mfem::ParMesh &pmesh = *nd.GetParMesh();
    const int dim = pmesh.Dimension();

    // The vertices: a first-order H1 space (its DOFs are the mesh vertices), and their coordinates.
    mfem::H1_FECollection h1_collection(1, dim);
    mfem::ParFiniteElementSpace h1(&pmesh, &h1_collection);
    mfem::ParDiscreteLinearOperator gradient(&h1, &nd);
    gradient.AddDomainInterpolator(new mfem::GradientInterpolator);
    gradient.Assemble();
    gradient.Finalize();
    std::unique_ptr<mfem::HypreParMatrix> G_mfem(gradient.ParallelAssemble());
    std::array<std::unique_ptr<mfem::HypreParVector>, 3> xyz;
    for (int c = 0; c < dim; ++c)
    {
        mfem::ParGridFunction coordinate(&h1);
        for (int vertex = 0; vertex < pmesh.GetNV(); ++vertex)
            coordinate(vertex) = pmesh.GetVertex(vertex)[c];
        xyz[static_cast<std::size_t>(c)].reset(coordinate.ParallelProject());
    }

    // S: the reduced rows (the owned free DOFs) of the MFEM true edge DOFs, with the orientation of each DOF.
    std::unordered_map<PetscInt, std::size_t> local_of_global;
    for (std::size_t l = 0; l < dofs_.n_local(); ++l)
        if (dofs_.owned(l))
            local_of_global.emplace(dofs_.global(l), l);
    const int rows = n_free_local();
    std::vector<int> I(static_cast<std::size_t>(rows) + 1, 0);
    std::vector<HYPRE_BigInt> J(static_cast<std::size_t>(rows));
    std::vector<double> data(static_cast<std::size_t>(rows));
    const std::size_t offset = dofs_.local_index(v, 0);
    for (int i = 0; i < rows; ++i)
    {
        const PetscInt global = dofs_.first_owned() + free_local_[static_cast<std::size_t>(i)];
        const std::size_t l = local_of_global.at(global);
        const int vdof = static_cast<int>(l - offset);
        J[static_cast<std::size_t>(i)] = nd.GetGlobalTDofNumber(vdof);
        data[static_cast<std::size_t>(i)] = dofs_.edge_sign(l);
        I[static_cast<std::size_t>(i) + 1] = i + 1;
    }
    HYPRE_BigInt row_starts[2] = {free_partition_[0], free_partition_[1]};
    mfem::HypreParMatrix S(dofs_.mesh().comm(), rows, n_free_global_, nd.GlobalTrueVSize(), I.data(), J.data(),
                           data.data(), row_starts, nd.GetTrueDofOffsets());
    mfem::HypreParMatrix *G = mfem::ParMult(&S, G_mfem.get(), true);

    // An ungauged curl-curl system (no conductivity, no tree gauge) has the discrete gradients in its kernel: K G = 0.
    // A random probe tells: |K G x| / |G x| against |K y| / |y|. AMS then drops its Poisson correction (singular).
    const MPI_Comm comm = dofs_.mesh().comm();
    // The probe on the interior vertices only: those whose edges are all free (|G| counts the edges of each vertex,
    // in the whole space and in the reduced one).
    const auto abs_copy = [](const mfem::HypreParMatrix &A)
    {
        auto copy = std::make_unique<mfem::HypreParMatrix>(A);
        hypre_ParCSRMatrix *h = *copy;
        for (hypre_CSRMatrix *block : {hypre_ParCSRMatrixDiag(h), hypre_ParCSRMatrixOffd(h)})
        {
            HYPRE_Complex *values = hypre_CSRMatrixData(block);
            for (HYPRE_Int k = 0; k < hypre_CSRMatrixNumNonzeros(block); ++k)
                values[k] = std::abs(values[k]);
        }
        return copy;
    };
    mfem::Vector all_edges(G_mfem->Height()), free_edges(G->Height()), count_all(G->Width()), count_free(G->Width());
    count_all = 0.0;
    count_free = 0.0;
    all_edges = 1.0;
    free_edges = 1.0;
    abs_copy(*G_mfem)->MultTranspose(all_edges, count_all);
    abs_copy(*G)->MultTranspose(free_edges, count_free);
    mfem::Vector x(G->Width()), Gx(G->Height()), KGx(G->Height()), y(G->Height()), Ky(G->Height());
    Gx = 0.0;
    KGx = 0.0;
    Ky = 0.0;
    x.Randomize(1 + dofs_.mesh().rank());
    for (int i = 0; i < x.Size(); ++i)
        if (count_free(i) < count_all(i) - 0.5)
            x(i) = 0.0;
    y.Randomize(7 + dofs_.mesh().rank());
    G->Mult(x, Gx);
    mfem_A_->Mult(Gx, KGx);
    mfem_A_->Mult(y, Ky);
    const auto norm = [comm](const mfem::Vector &a) { return std::sqrt(mfem::InnerProduct(comm, a, a)); };
    const double gradient_gain = norm(KGx) / std::max(norm(Gx), 1e-300);
    const double typical_gain = norm(Ky) / std::max(norm(y), 1e-300);
    const bool singular = gradient_gain <= 1e-8 * typical_gain;

    kernel_G_.reset();
    kernel_L_.reset();
    kernel_amg_.reset();
    kernel_cg_.reset();
    if (singular)
    {
        mfem::Vector interior(G->Width());
        for (int i = 0; i < interior.Size(); ++i)
            interior(i) = count_free(i) < count_all(i) - 0.5 ? 0.0 : 1.0;
        std::unique_ptr<mfem::HypreParMatrix> Gt(G->Transpose());
        Gt->ScaleRows(interior);
        kernel_G_.reset(Gt->Transpose());
        Gt.reset(kernel_G_->Transpose());
        std::unique_ptr<mfem::HypreParMatrix> L(mfem::ParMult(Gt.get(), kernel_G_.get(), true));
        // The boundary vertices (no column of G): identity rows, the Laplacian stays definite.
        mfem::SparseMatrix boundary(interior.Size());
        for (int i = 0; i < interior.Size(); ++i)
            boundary.Set(i, i, 1.0 - interior(i));
        boundary.Finalize();
        mfem::HypreParMatrix B(comm, G->GetGlobalNumCols(), G->ColPart(), &boundary);
        kernel_L_.reset(mfem::Add(1.0, *L, 1.0, B));
        kernel_amg_ = std::make_unique<mfem::HypreBoomerAMG>(*kernel_L_);
        kernel_amg_->SetPrintLevel(0);
        kernel_cg_ = std::make_unique<mfem::CGSolver>(comm);
        kernel_cg_->SetRelTol(1e-12);
        kernel_cg_->SetAbsTol(0.0);
        kernel_cg_->SetMaxIter(500);
        kernel_cg_->SetOperator(*kernel_L_);
        kernel_cg_->SetPreconditioner(*kernel_amg_);
    }
    auto ams = std::make_unique<mfem::HypreAMS>(*mfem_A_, G, xyz[0].release(), xyz[1].release(),
                                                dim == 3 ? xyz[2].release() : nullptr);
    ams->SetPrintLevel(0);
    if (singular)
        ams->SetSingularProblem();
    pytec3d::logging::info("AMS preconditioner: {} edge DOFs, {} vertices{}.", n_free_global_, G->GetGlobalNumCols(),
                           singular ? ", singular curl-curl (no gauge): no Poisson correction" : "");
    return ams;
}

void ReducedSystem::setup_mfem_solver(double rtol)
{
    const LinearSolverSettings &linear = configuration_.settings().linear;
    const MPI_Comm comm = dofs_.mesh().comm();
    const bool direct_pc = linear.pc_type == "lu" || linear.pc_type == "cholesky";
    const bool numeric = mfem_solver_ && !khat_new_;

    auto make_mumps = [&]() -> std::unique_ptr<mfem::Solver>
    {
#ifdef MFEM_USE_MUMPS
        auto mumps = std::make_unique<mfem::MUMPSSolver>(comm);
        mumps->SetPrintLevel(0);
        mumps->SetMatrixSymType(linear.pc_type == "cholesky" ? mfem::MUMPSSolver::MatType::SYMMETRIC_INDEFINITE
                                                             : mfem::MUMPSSolver::MatType::UNSYMMETRIC);
        mumps->SetReorderingReuse(true);
        return mumps;
#else
        pytec3d::logging::error_and_abort("Linear solver (backend mfem): this MFEM has no MUMPS (MFEM_USE_MUMPS); use "
                                          "backend 'petsc'.");
        return nullptr;
#endif
    };

    if (linear.ksp_type == "preonly")
    {
        if (!direct_pc)
            pytec3d::logging::error_and_abort("Linear solver (backend mfem): ksp_type 'preonly' needs pc_type lu or "
                                              "cholesky.");
        if (!numeric)
            mfem_solver_ = make_mumps();
        mfem_solver_->SetOperator(*mfem_A_); // MUMPS: analysis reused on the same pattern
        return;
    }

    // A Krylov solver: a new preconditioner on new values (MUMPS: numeric only).
    if (direct_pc)
    {
        if (!numeric || !mfem_prec_)
            mfem_prec_ = make_mumps();
        mfem_prec_->SetOperator(*mfem_A_);
    }
    else if (linear.pc_type == "hypre" || linear.pc_type == "boomeramg")
    {
        auto amg = std::make_unique<mfem::HypreBoomerAMG>(*mfem_A_);
        amg->SetPrintLevel(0);
        mfem_prec_ = std::move(amg);
    }
    else if (linear.pc_type == "ams")
        mfem_prec_ = make_ams();
    else if (linear.pc_type == "jacobi")
        mfem_prec_ = std::make_unique<mfem::HypreSmoother>(*mfem_A_, mfem::HypreSmoother::Jacobi);
    else if (linear.pc_type == "none")
        mfem_prec_.reset();
    else
        pytec3d::logging::error_and_abort("Linear solver (backend mfem): unknown pc_type '" + linear.pc_type +
                                          "' (lu, cholesky, boomeramg, ams, jacobi, none).");
    auto solver = krylov(linear.ksp_type, comm);
    solver->SetRelTol(rtol);
    solver->SetAbsTol(linear.absolute_tolerance);
    solver->SetMaxIter(linear.max_iterations);
    solver->SetPrintLevel(mfem::IterativeSolver::PrintLevel().None());
    solver->iterative_mode = linear.initial_guess_nonzero;
    solver->SetOperator(*mfem_A_);
    if (mfem_prec_)
        solver->SetPreconditioner(*mfem_prec_);
    mfem_solver_ = std::move(solver);
}

std::pair<unsigned int, double> ReducedSystem::solve_mfem(const mfem::Vector &rhs, mfem::Vector &solution, double rtol)
{
    if (!mfem_A_)
        pytec3d::logging::error_and_abort("ReducedSystem::solve: no matrix.");
    if (n_free_global_ == 0)
        pytec3d::logging::error_and_abort("ReducedSystem::solve: every DOF is constrained.");
    if (configuration_.take_report())
    {
        const LinearSolverSettings &linear = configuration_.settings().linear;
        const bool direct = linear.pc_type == "lu" || linear.pc_type == "cholesky";
        pytec3d::logging::table_root(pytec3d::report::format_solver_configuration_table(
            configuration_.settings(), configuration_.nonlinear(), linear.ksp_type, linear.pc_type,
            direct ? "mumps (MFEM MUMPSSolver)" : "-"));
    }

    if (!mfem_solver_ || khat_new_)
    {
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
    if (last_ != Factorization::reused)
        setup_mfem_solver(rtol);
    else if (auto *iterative = dynamic_cast<mfem::IterativeSolver *>(mfem_solver_.get()))
        iterative->SetRelTol(rtol);
    khat_new_ = false;
    khat_changed_ = false;

    solution.SetSize(n_free_local());
    if (!configuration_.settings().linear.initial_guess_nonzero)
        solution = 0.0;
    if (kernel_cg_ && dynamic_cast<mfem::IterativeSolver *>(mfem_solver_.get()))
    {
        // Ungauged: f - G (G^T G)^-1 G^T f, the right-hand side without its part in the kernel (the gradients).
        mfem::Vector Gtf(kernel_G_->Width()), phi(kernel_G_->Width()), consistent(rhs);
        Gtf = 0.0;
        kernel_G_->MultTranspose(rhs, Gtf);
        phi = 0.0;
        kernel_cg_->Mult(Gtf, phi);
        kernel_G_->AddMult(phi, consistent, -1.0);
        mfem_solver_->Mult(consistent, solution);
        consistent_rhs_ = consistent;
        projected_ = true;
    }
    else
    {
        mfem_solver_->Mult(rhs, solution);
        projected_ = false;
    }
    if (auto *iterative = dynamic_cast<mfem::IterativeSolver *>(mfem_solver_.get()))
    {
        mfem_converged_ = iterative->GetConverged();
        return {static_cast<unsigned int>(iterative->GetNumIterations()), iterative->GetFinalNorm()};
    }
    mfem_converged_ = true;
    return {1u, 0.0};
}

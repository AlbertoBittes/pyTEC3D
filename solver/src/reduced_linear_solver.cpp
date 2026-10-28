/**
 * @file reduced_linear_solver.cpp
 * @brief ReducedPetscLinearSolver: restricted PETSc solves keeping the reduced system and its factorization.
 */
#include "reduced_linear_solver.h"

#include <libmesh/petsc_macro.h>
#include <libmesh/petsc_matrix_base.h>
#include <libmesh/petsc_solver_exception.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/solver_configuration.h>

#include "logging.h"

using libMesh::libmesh_ignore;
using Base = libMesh::PetscLinearSolver<libMesh::Number>;

ReducedPetscLinearSolver::ReducedPetscLinearSolver(const libMesh::Parallel::Communicator &comm)
    : Base(comm)
{
}

void ReducedPetscLinearSolver::clear()
{
    clear_reduced_system();
    Base::clear();
}

void ReducedPetscLinearSolver::init(const char *name)
{
    // The prefix of the first init, given again to the KSP of each full factorization (ksp() calls init() without).
    if (!this->initialized())
    {
        has_prefix_ = name != nullptr;
        prefix_ = name ? std::string(name) : std::string();
    }
    Base::init(name);
}

void ReducedPetscLinearSolver::init_systems(const libMesh::System &system)
{
    system_ = &system;
    Base::init_systems(system);
}

void ReducedPetscLinearSolver::restrict_solve_to(const std::vector<unsigned int> *const dofs,
                                                 const libMesh::SubsetSolveMode subset_solve_mode)
{
    if (subset_solve_mode != libMesh::SUBSET_ZERO)
        pytec3d::logging::error_and_abort("ReducedPetscLinearSolver::restrict_solve_to: only SUBSET_ZERO is supported.");
    const bool has_subset = dofs != nullptr;
    if (has_subset == has_subset_ && (!has_subset || *dofs == dofs_))
        return; // the same subset: the reduced system is kept
    has_subset_ = has_subset;
    dofs_ = has_subset ? *dofs : std::vector<unsigned int>();
    clear_reduced_system();
}

void ReducedPetscLinearSolver::clear_reduced_system()
{
    is_.reset_to_zero();
    scatter_.reset_to_zero();
    submat_.reset_to_zero();
    subrhs_.reset_to_zero();
    subsolution_.reset_to_zero();
    parent_ = nullptr;
    parent_nonzeros_ = -1.0;
    pattern_valid_ = false;
}

void ReducedPetscLinearSolver::build(Mat parent, Vec rhs)
{
    // A new KSP: the operators of the previous one may have other sizes.
    const std::string prefix = prefix_;
    const bool has_prefix = has_prefix_;
    this->clear();
    this->init(has_prefix ? prefix.c_str() : nullptr);
    if (system_)
        Base::init_systems(*system_);

    PetscInt n_local = 0;
    PetscInt *petsc_dofs = nullptr;
    if (has_subset_)
    {
        n_local = static_cast<PetscInt>(dofs_.size());
        LibmeshPetscCall(PetscMalloc1(n_local, &petsc_dofs));
        for (PetscInt i = 0; i < n_local; ++i)
            petsc_dofs[i] = static_cast<PetscInt>(dofs_[i]);
    }
    else
    {
        PetscInt first = 0, last = 0;
        LibmeshPetscCall(VecGetOwnershipRange(rhs, &first, &last));
        n_local = last - first;
        LibmeshPetscCall(PetscMalloc1(n_local, &petsc_dofs));
        for (PetscInt i = 0; i < n_local; ++i)
            petsc_dofs[i] = first + i;
    }
    LibmeshPetscCall(ISCreateGeneral(this->comm().get(), n_local, petsc_dofs, PETSC_OWN_POINTER, is_.get()));

    LibmeshPetscCall(VecCreate(this->comm().get(), subrhs_.get()));
    LibmeshPetscCall(VecSetSizes(subrhs_, n_local, PETSC_DECIDE));
    LibmeshPetscCall(VecSetFromOptions(subrhs_));
    LibmeshPetscCall(VecDuplicate(subrhs_, subsolution_.get()));
    LibmeshPetscCall(VecScatterCreate(rhs, is_, subrhs_, nullptr, scatter_.get()));
    LibmeshPetscCall(LibMeshCreateSubMatrix(parent, is_, is_, MAT_INITIAL_MATRIX, submat_.get()));

    MatInfo info;
    LibmeshPetscCall(MatGetInfo(parent, MAT_GLOBAL_SUM, &info));
    parent_ = parent;
    parent_nonzeros_ = info.nz_used;
    pattern_valid_ = true;
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
    auto *matrix = libMesh::cast_ptr<libMesh::PetscMatrixBase<libMesh::Number> *>(&matrix_in);
    auto *solution = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&solution_in);
    auto *rhs = libMesh::cast_ptr<libMesh::PetscVector<libMesh::Number> *>(&rhs_in);
    this->init();
    matrix->close();
    solution->close();
    rhs->close();

    // The level, the same on all the ranks: full without a reduced system of this subset and pattern, or when the
    // matrix has another number of nonzeros; else reused (the same matrix) or numeric.
    const Mat parent = matrix->mat();
    MatInfo info;
    LibmeshPetscCall(MatGetInfo(parent, MAT_GLOBAL_SUM, &info));
    bool full = !submat_ || !pattern_valid_ || parent != parent_ || info.nz_used != parent_nonzeros_;
    this->comm().max(full);
    const Factorization level = full ? Factorization::full
                                : this->same_preconditioner ? Factorization::reused
                                                            : Factorization::numeric;
    if (level == Factorization::full)
        build(parent, rhs->vec());
    else if (level == Factorization::numeric)
        // The values of the same nonzeros: PETSc keeps the ordering and the analysis of the factorization.
        LibmeshPetscCall(LibMeshCreateSubMatrix(parent, is_, is_, MAT_REUSE_MATRIX, submat_.get()));

    const KSP ksp = this->ksp();
    LibmeshPetscCall(KSPSetOperators(ksp, submat_, submat_));
    LibmeshPetscCall(KSPSetReusePreconditioner(ksp, level == Factorization::reused ? PETSC_TRUE : PETSC_FALSE));

    libMesh::VecScatterBeginEnd(this->comm(), scatter_, rhs->vec(), subrhs_, INSERT_VALUES, SCATTER_FORWARD);
    libMesh::VecScatterBeginEnd(this->comm(), scatter_, solution->vec(), subsolution_, INSERT_VALUES, SCATTER_FORWARD);

    // As PetscLinearSolver::solve: the tolerances of the caller or of the settings, the command line, then the
    // pyTEC3D settings.
    const double rel_tol = this->get_real_solver_setting("rel_tol", tol);
    const double abs_tol = this->get_real_solver_setting("abs_tol", std::nullopt, static_cast<double>(PETSC_DEFAULT));
    const double max_its = this->get_int_solver_setting("max_its", m_its);
    LibmeshPetscCall(KSPSetTolerances(ksp, rel_tol, abs_tol, PETSC_DEFAULT, static_cast<PetscInt>(max_its)));
    LibmeshPetscCall(KSPSetFromOptions(ksp));
    if (this->_solver_configuration)
        this->_solver_configuration->configure_solver();
    LibmeshPetscCall(KSPSolve(ksp, subrhs_, subsolution_));

    PetscInt its = 0;
    PetscReal residual = 0.0;
    LibmeshPetscCall(KSPGetIterationNumber(ksp, &its));
    LibmeshPetscCall(KSPGetResidualNorm(ksp, &residual));

    LibmeshPetscCall(VecZeroEntries(solution->vec()));
    libMesh::VecScatterBeginEnd(this->comm(), scatter_, subsolution_, solution->vec(), INSERT_VALUES, SCATTER_REVERSE);

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

/**
 * @file dof_system.h
 * @brief DofSystem: the DOF system of a model, its time scheme and its solve (linear, or Newton).
 *
 * The DOF system solves M a + C v + K u = F at the end of each step (time_step.h) on the unknowns of fem::DofSpace,
 * with the constraints of the step (constraints.h) through the reduced system (reduced_system.h). It is assembled the
 * way an MFEM code assembles its forms:
 *   - K, C and M once per model: the element matrices of the element types into an L-level SparseMatrix, then the
 *     T-level HypreParMatrix P_L^T A P_L; the operators of the static nodal and element fields likewise, as sparse
 *     operators from the L-vectors of the solution to the field values of the elements;
 *   - K_eff = K + c_v C + c_a M (hypre Add) when the step size or the scheme changes; the history of the step
 *     F - C h_v - M h_a as products of T-vectors;
 *   - the right-hand side at every step (it follows the loads), an L-vector assembled to a T-vector;
 *   - nonlinear (laws on some elements): Newton on R(u) = K_eff u - F + R_laws(u), with the Jacobian K_eff + J_laws(u);
 *     only the elements with laws are evaluated at each iteration. The algorithm, line search and convergence tests are
 *     those of libMesh's NewtonSolver, which the "nonlinear" solver settings describe.
 *
 * Vectors: the solution u, its velocity v and acceleration a at the end and at the start of the step, T-vectors;
 * their L-vectors (DofSpace::gather) for the elements.
 */
#pragma once

#include <array>
#include <memory>
#include <unordered_map>
#include <vector>

#include <mfem.hpp>

#include "constraints.h"
#include "element_builder.h"
#include "reduced_system.h"
#include "time_step.h"

class pyTEC3DEquationSystems;

class DofSystem
{
public:
    enum VectorId
    {
        SOLUTION,
        OLD_SOLUTION,
        VELOCITY,
        ACCELERATION,
        OLD_VELOCITY,
        OLD_ACCELERATION,
        N_VECTORS
    };

    DofSystem(pyTEC3DEquationSystems &equation_systems, bool nonlinear);
    ~DofSystem();
    DofSystem(const DofSystem &) = delete;
    DofSystem &operator=(const DofSystem &) = delete;

    bool nonlinear() const { return nonlinear_; }

    /** @brief The solver settings of the next solve; a change of the linear solver settings refactorizes. */
    void apply_solver_settings(const SolverSettings &settings);
    /** @brief The step: time-step state, start-of-step vectors, constraints. */
    void init_time_step(const TimeStepInput &time_step);
    /** @brief K, C, M and the static field operators, once per model. */
    void assemble_operators();
    /** @brief Assembles and solves the step; then the time derivatives (and the rates of the law outputs). */
    void solve();

    /** @brief The nodal fields at the points / element fields at the cells of the local elements: per variable of the
     * layout, the sums of the values (n_components per point or cell) and the number of elements that gave one. */
    void accumulate_fields(pytec3d::FieldKind kind, std::vector<std::vector<double>> &sums,
                           std::vector<std::vector<double>> &counts) const;
    /** @brief The same for the DOF fields: the solution at the points of the elements. */
    void accumulate_dof_fields(std::vector<std::vector<double>> &sums, std::vector<std::vector<double>> &counts) const;

    /** @brief The L-vector of a vector (gathered on demand). */
    const std::vector<double> &local(VectorId id) const;
    const mfem::Vector &vector(VectorId id) const { return vectors_[id]; }

    const ReducedSystem &reduced() const { return *reduced_; }

private:
    void invalidate(VectorId id) { gathered_[id] = false; }
    LocalVectors local_vectors() const;

    /** @brief Adds an element matrix / vector, in the element's shape functions, to an L-level matrix / vector. */
    void add_element_matrix(mfem::SparseMatrix &A, int e, const ElementBuilder::Matrix &m) const;
    /** @brief A zero, finalized L-level matrix with the sparsity @p pattern (DOF-DOF connectivity). */
    std::unique_ptr<mfem::SparseMatrix> matrix_with_pattern(const mfem::Table &pattern) const;
    /** @brief The DOF-DOF connectivity of @p elements (all the local elements if null). */
    mfem::Table connectivity(const std::vector<int> *elements) const;
    void add_element_vector(mfem::Vector &b, int e, const ElementBuilder::Vector &v) const;
    /** @brief The same, safe on several threads (atomic additions). */
    void add_element_vector_atomic(mfem::Vector &b, int e, const ElementBuilder::Vector &v) const;
    /** @brief The nodal loads of the DOF fields (point arrays) times @p scale, added to the T-vector @p b. */
    void add_nodal_loads(mfem::Vector &b, double scale) const;

    /** @brief True when K_eff must be formed again for this step (logs why). */
    bool check_update_system_matrix();
    /** @brief K_eff = K + c_v C + c_a M of the current scheme. */
    void update_effective_matrix();
    /** @brief F of the step: element right-hand sides, nodal loads, history - C h_v - M h_a. */
    void assemble_rhs(mfem::Vector &F) const;
    void update_time_derivatives();

    void solve_linear();
    /** @brief The law terms of the elements with laws at the current solution: their residual (T-vector) and Jacobian. */
    void assemble_laws(bool jacobian, mfem::Vector *residual);
    /** @brief R(u) = K_eff u - F + R_laws(u) at the current solution, and ||P^T R||. */
    double assemble_residual();
    void solve_nonlinear();
    double line_search(double tol, double last_residual, double &current_residual, const mfem::Vector &du);
    bool test_convergence(double residual, double step, bool linear_solve_finished);
    void record_newton_iterate(double residual);
    void update_law_rates();

    double norm(const mfem::Vector &x) const;
    /** @brief The Jacobians at the mesh points of every element (Nedelec outputs), computed at the first use. */
    const pytec3d::fem::QuadratureData &point_quadrature() const;

    pyTEC3DEquationSystems &es_;
    const pytec3d::fem::DofSpace &dofs_;
    bool nonlinear_;
    std::unique_ptr<ReducedSystem> reduced_;
    pytec3d::CouplingPairs coupling_pairs_;

    std::array<mfem::Vector, N_VECTORS> vectors_;
    mutable std::array<std::vector<double>, N_VECTORS> local_{};
    mutable std::array<bool, N_VECTORS> gathered_{};

    // The operators of the model.
    std::unique_ptr<mfem::HypreParMatrix> K_, C_, M_;
    std::unique_ptr<mfem::HypreParMatrix> Keff_owned_; ///< K + c_v C + c_a M of a transient step.
    const mfem::HypreParMatrix *Keff_ = nullptr;        ///< Of the scheme of the current step (K_ on a static one).
    std::unique_ptr<mfem::HypreParMatrix> jacobian_; ///< Nonlinear: K_eff + J_laws at the current iterate.
    /// The static field operators, by kind (nodal, element) and dt() order: per element, a dense block (field values of
    /// the element x its DOFs, row-major, in the element's shape functions), at block_offsets_[kind][order][e] of
    /// field_blocks_ (an empty block: equal offsets). The columns are the element DOFs (DofSpace::element_dofs): only the
    /// values are stored.
    std::array<std::array<std::vector<double>, 3>, 2> field_blocks_;
    std::array<std::array<std::vector<std::size_t>, 3>, 2> block_offsets_;
    std::array<std::vector<int>, 2> field_offsets_; ///< First field value of each element.
    /// The sparsity of K, C, M (all the elements) and of the Jacobian of the laws (the elements with laws): CSR on the
    /// L-DOFs, from the element-DOF connectivity, as MFEM's forms precompute it.
    mfem::Table pattern_, law_pattern_;
    std::array<std::vector<int>, 2> evaluated_field_elements_; ///< Elements with fields evaluated at every output.
    std::vector<int> law_elements_;                             ///< Elements with laws selected.
    mutable std::unique_ptr<pytec3d::fem::QuadratureData> point_quadrature_;
    bool has_C_ = false, has_M_ = false;
    bool operators_assembled_ = false;

    mfem::Vector F_;        ///< Right-hand side of the step.
    mfem::Vector residual_; ///< Nonlinear: R at the current iterate.
    mfem::Vector reduced_residual_;
    bool settings_changed_ = false;
    bool constraints_changed_ = false;
    bool first_solve_ = true;
    SolverSettings settings_;

    // Newton (libMesh NewtonSolver): its state across the steps of the model.
    unsigned int solve_result_ = 0;
    unsigned int inner_iterations_ = 0;
    unsigned int outer_iterations_ = 0;
    double max_residual_norm_ = 0.0;
    double max_solution_norm_ = 0.0;
    bool has_previous_iterate_ = false;
    mfem::Vector previous_iterate_;
    double full_newton_step_ = -1.0;

    /** @brief The first and second time derivatives of the law outputs under dt() of one element. */
    struct LawRates
    {
        std::array<ElementBuilder::Vector, 2> v; ///< q'
        std::array<ElementBuilder::Vector, 2> a; ///< q''
    };
    std::unordered_map<int, LawRates> law_rates_;     ///< At the end of the step, by local element.
    std::unordered_map<int, LawRates> old_law_rates_; ///< At the start of the step.
};

/**
 * @file dof_system.cpp
 * @brief DofSystem: operators, right-hand side, linear and Newton solves, time derivatives, fields (dof_system.h).
 */
#include "dof_system.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>

#include <omp.h>

#include "element_type.h"
#include "logging.h"
#include "model_element_types.h"
#include "pytec3d_equation_systems.h"
#include "pytec3d_model.h"
#include "registries.h"
#include "report_tables.h"
#include "solver_diagnostics.h"

using pytec3d::la::Number;
using pytec3d::la::Real;

namespace
{
    /** @brief The results of the Newton solver: the flags of libMesh's DiffSolver::SolveResult. */
    enum NewtonResult : unsigned int
    {
        INVALID_SOLVE_RESULT = 0,
        CONVERGED_ABSOLUTE_RESIDUAL = 2,
        CONVERGED_RELATIVE_RESIDUAL = 4,
        CONVERGED_ABSOLUTE_STEP = 8,
        CONVERGED_RELATIVE_STEP = 16,
        DIVERGED_NO_REASON = 32,
        DIVERGED_MAX_NONLINEAR_ITERATIONS = 64,
        DIVERGED_BACKTRACKING_FAILURE = 128,
        DIVERGED_LINEAR_SOLVER_FAILURE = 256,
    };

    /** @brief libMesh's TOLERANCE in double precision: the line search runs to sqrt(TOLERANCE). */
    constexpr double kTolerance = 1.e-8;
    /** @brief libMesh's DiffSolver::max_linear_iterations: a linear solve of as many iterations is unfinished. */
    constexpr unsigned int kMaxLinearIterations = 1000;

    void accumulate(std::vector<double> &sums, std::vector<double> &counts, std::int64_t at, unsigned int n_components,
                    const double *values)
    {
        double *sum = sums.data() + static_cast<std::size_t>(at) * n_components;
        for (unsigned int c = 0; c < n_components; ++c)
        {
#pragma omp atomic
            sum[c] += values[c];
        }
#pragma omp atomic
        counts[static_cast<std::size_t>(at)] += 1.0;
    }

    /**
     * @brief Adds the values at the points of element @p e of the field of variable @p var, whose coefficients in the
     * element's shape functions are @p coefficients, to (sums, counts). The shape functions at the points come from
     * the reference tables of the cell type's points; a Nedelec field also needs the Jacobians there (@p points).
     */
    void accumulate_at_points(const pytec3d::fem::FEMesh &mesh, const pytec3d::fem::FECatalog &catalog,
                              const pytec3d::fem::ReferenceTables &tables, const pytec3d::fem::QuadratureData *points_data,
                              int e, const pytec3d::fem::Variable &var, const double *coefficients, std::vector<double> &sums,
                              std::vector<double> &counts)
    {
        int n_points = 0;
        const std::int64_t *points = mesh.element_points(e, n_points);
        const mfem::IntegrationRule &rule = mesh.point_rule(e);
        const mfem::FiniteElement *fe = catalog.element(var.family, var.order, mesh.mesh().GetElementGeometry(e));
        const pytec3d::fem::ReferenceTable &table = tables.get(*fe, rule);
        const int n = table.n_dofs;
        const int dim = table.dim;
        double values[3] = {0.0, 0.0, 0.0};
        if (var.edge())
        {
            // phi = J^-T phi_ref at each point.
            const double *J = points_data->J(e);
            const double *det = points_data->det(e);
            for (int k = 0; k < n_points; ++k)
            {
                double reference[3] = {0.0, 0.0, 0.0};
                for (int i = 0; i < n; ++i)
                    for (int c = 0; c < dim; ++c)
                        reference[c] += coefficients[i] * table.B[static_cast<std::size_t>((k * n + i) * dim + c)];
                const double *Jk = J + 9 * k;
                const double d = det[k];
                if (dim == 3)
                {
                    const double K[9] = {(Jk[4] * Jk[8] - Jk[5] * Jk[7]) / d, (Jk[5] * Jk[6] - Jk[3] * Jk[8]) / d,
                                         (Jk[3] * Jk[7] - Jk[4] * Jk[6]) / d, (Jk[2] * Jk[7] - Jk[1] * Jk[8]) / d,
                                         (Jk[0] * Jk[8] - Jk[2] * Jk[6]) / d, (Jk[1] * Jk[6] - Jk[0] * Jk[7]) / d,
                                         (Jk[1] * Jk[5] - Jk[2] * Jk[4]) / d, (Jk[2] * Jk[3] - Jk[0] * Jk[5]) / d,
                                         (Jk[0] * Jk[4] - Jk[1] * Jk[3]) / d};
                    for (int r = 0; r < 3; ++r)
                        values[r] = K[3 * r] * reference[0] + K[3 * r + 1] * reference[1] + K[3 * r + 2] * reference[2];
                }
                else
                {
                    values[0] = (Jk[4] * reference[0] - Jk[3] * reference[1]) / d;
                    values[1] = (-Jk[1] * reference[0] + Jk[0] * reference[1]) / d;
                }
                accumulate(sums, counts, points[k], var.n_components, values);
            }
            return;
        }
        for (int k = 0; k < n_points; ++k)
        {
            const double *B = table.B.data() + static_cast<std::size_t>(k * n);
            for (unsigned int c = 0; c < var.n_components; ++c)
            {
                values[c] = 0.0;
                for (int i = 0; i < n; ++i)
                    values[c] += coefficients[static_cast<std::size_t>(c) * static_cast<std::size_t>(n) + i] * B[i];
            }
            accumulate(sums, counts, points[k], var.n_components, values);
        }
    }

    /** @brief Elements computed in parallel, then inserted in order: the work of a chunk at a time in memory. */
    constexpr int kChunk = 4096;

    /**
     * @brief fn(builder, i) for i in [begin, end), on the OpenMP threads, each with its own builder. An exception of a
     * thread is rethrown after the loop (an abort of the solver reaches the caller as from a serial loop).
     */
    template <typename Fn>
    void parallel_for(const pyTEC3DEquationSystems &es, const LocalVectors *vectors, int begin, int end, Fn &&fn)
    {
        std::exception_ptr failure;
#pragma omp parallel
        {
            ElementBuilder builder(es);
            if (vectors)
                builder.set_vectors(*vectors);
#pragma omp for schedule(dynamic, 16)
            for (int i = begin; i < end; ++i)
            {
                try
                {
                    fn(builder, i);
                }
                catch (...)
                {
#pragma omp critical(pytec3d_parallel_failure)
                    if (!failure)
                        failure = std::current_exception();
                }
            }
        }
        if (failure)
            std::rethrow_exception(failure);
    }

    std::unique_ptr<mfem::HypreParMatrix> add(double a, const mfem::HypreParMatrix &A, double b, const mfem::HypreParMatrix &B)
    {
        return std::unique_ptr<mfem::HypreParMatrix>(mfem::Add(a, A, b, B));
    }
} // namespace

DofSystem::DofSystem(pyTEC3DEquationSystems &equation_systems, bool nonlinear)
    : es_(equation_systems), dofs_(equation_systems.dof_space()), nonlinear_(nonlinear),
      reduced_(std::make_unique<ReducedSystem>(equation_systems.dof_space()))
{
    for (mfem::Vector &x : vectors_)
    {
        x.SetSize(static_cast<int>(dofs_.n_owned()));
        x = 0.0;
    }
    gathered_.fill(false);
    coupling_pairs_ = pytec3d::coupling_pairs(es_.model().arrays(), es_.mesh().n_points());
}

DofSystem::~DofSystem() = default;

double DofSystem::norm(const mfem::Vector &x) const
{
    return std::sqrt(mfem::InnerProduct(es_.mesh().comm(), x, x));
}

const std::vector<double> &DofSystem::local(VectorId id) const
{
    if (!gathered_[id])
    {
        dofs_.gather(vectors_[id], local_[id]);
        gathered_[id] = true;
    }
    return local_[id];
}

LocalVectors DofSystem::local_vectors() const
{
    LocalVectors v;
    v.solution = &local(SOLUTION);
    v.old_solution = &local(OLD_SOLUTION);
    v.velocity = &local(VELOCITY);
    v.acceleration = &local(ACCELERATION);
    v.old_velocity = &local(OLD_VELOCITY);
    v.old_acceleration = &local(OLD_ACCELERATION);
    return v;
}

// ---- settings, time step, constraints --------------------------------------------------------------------------

void DofSystem::apply_solver_settings(const SolverSettings &settings)
{
    settings_ = settings;
    if (reduced_->set_settings(settings, nonlinear_))
        settings_changed_ = true;
}

void DofSystem::init_time_step(const TimeStepInput &time_step_input)
{
    es_.prepare_time_step_solve(time_step_input);
    const TimeStep *step = es_.current_time_step();
    if (step->updated)
    {
        // The end of the previous step is the start of this one.
        vectors_[OLD_SOLUTION] = vectors_[SOLUTION];
        vectors_[OLD_VELOCITY] = vectors_[VELOCITY];
        vectors_[OLD_ACCELERATION] = vectors_[ACCELERATION];
        invalidate(OLD_SOLUTION);
        invalidate(OLD_VELOCITY);
        invalidate(OLD_ACCELERATION);
        if (nonlinear_)
            old_law_rates_ = law_rates_;
        pytec3d::logging::info("Advanced to new time step");
        pytec3d::logging::table_root(pytec3d::report::format_time_step_table(*step, es_.time_order()));
    }

    // The constraints of a new step; the same step solved again keeps them.
    if (step->updated || !reduced_->has_constraints())
    {
        pytec3d::ConstraintInputs inputs;
        inputs.mesh = &es_.mesh();
        inputs.dofs = &dofs_;
        inputs.arrays = &es_.model().arrays();
        inputs.element_types = &es_.element_types();
        inputs.couplings = &es_.get_model_data().dof_couplings();
        inputs.pairs = &coupling_pairs_;
        inputs.node_rotations = &es_.node_rotations();
        inputs.static_step = !step->transient_effects;
        const ReducedSystem::Change change = reduced_->set_constraints(pytec3d::build_constraints(inputs));
        constraints_changed_ = constraints_changed_ || change == ReducedSystem::Change::rows || change == ReducedSystem::Change::dofs;
    }
    else
        pytec3d::logging::info("Reusing cached DOF constraints (time step unchanged)");
}

bool DofSystem::check_update_system_matrix()
{
    if (constraints_changed_ || settings_changed_)
    {
        pytec3d::logging::info(constraints_changed_ ? "New constraints: updating system matrix"
                                                    : "New linear solver settings: updating system matrix");
        constraints_changed_ = false;
        settings_changed_ = false;
        return true;
    }
    bool update = false;
    const TimeStep *current = es_.current_time_step();
    const TimeStep *old = es_.previous_time_step();
    if (first_solve_ || !Keff_ || !current)
        update = true;
    else if (current->updated)
    {
        if (!old)
            update = true;
        else
        {
            const auto differs = [](Real a, Real b) { return std::abs(a - b) > std::numeric_limits<Real>::epsilon(); };
            update = differs(current->dt, old->dt) || differs(current->theta, old->theta) ||
                     differs(current->beta, old->beta) || differs(current->gamma, old->gamma) ||
                     current->transient_effects != old->transient_effects;
        }
    }
    pytec3d::logging::info(update ? "Updating system matrix for current time-step" : "Reusing system matrix from previous time step");
    return update;
}

// ---- operators -------------------------------------------------------------------------------------------------

void DofSystem::add_element_matrix(mfem::SparseMatrix &A, int e, const ElementBuilder::Matrix &m) const
{
    const int *local = nullptr;
    const double *sign = nullptr;
    unsigned int n = 0;
    dofs_.element_dofs(e, local, sign, n);
    thread_local mfem::Array<int> indices;
    thread_local mfem::DenseMatrix block;
    indices.SetSize(static_cast<int>(n));
    block.SetSize(static_cast<int>(n));
    for (unsigned int i = 0; i < n; ++i)
    {
        indices[static_cast<int>(i)] = local[i];
        for (unsigned int j = 0; j < n; ++j)
            block(static_cast<int>(i), static_cast<int>(j)) = sign[i] * sign[j] * m(i, j);
    }
    // Into the precomputed pattern of a finalized matrix.
    A.AddSubMatrix(indices, indices, block, 0);
}

mfem::Table DofSystem::connectivity(const std::vector<int> *elements) const
{
    // element -> DOFs, then DOF -> elements -> DOFs: the DOFs coupled by an element (MFEM's Table algebra).
    const int n_elements = elements ? static_cast<int>(elements->size()) : es_.mesh().n_local_elements();
    mfem::Table element_dofs;
    element_dofs.MakeI(n_elements);
    for (int pass = 0; pass < 2; ++pass)
    {
        for (int k = 0; k < n_elements; ++k)
        {
            const int e = elements ? (*elements)[static_cast<std::size_t>(k)] : k;
            const int *local = nullptr;
            const double *sign = nullptr;
            unsigned int n = 0;
            dofs_.element_dofs(e, local, sign, n);
            for (unsigned int i = 0; i < n; ++i)
            {
                if (pass == 0)
                    element_dofs.AddAColumnInRow(k);
                else
                    element_dofs.AddConnection(k, local[i]);
            }
        }
        if (pass == 0)
            element_dofs.MakeJ();
    }
    element_dofs.ShiftUpI();
    element_dofs.Finalize();
    mfem::Table dof_elements;
    mfem::Transpose(element_dofs, dof_elements, static_cast<int>(dofs_.n_local()));
    mfem::Table out;
    mfem::Mult(dof_elements, element_dofs, out);
    return out;
}

std::unique_ptr<mfem::SparseMatrix> DofSystem::matrix_with_pattern(const mfem::Table &pattern) const
{
    const int n = static_cast<int>(dofs_.n_local());
    const int nnz = pattern.Size_of_connections();
    // new[]: the matrix owns and deletes the arrays.
    int *I = new int[static_cast<std::size_t>(n) + 1];
    int *J = new int[static_cast<std::size_t>(nnz)];
    std::copy(pattern.GetI(), pattern.GetI() + n + 1, I);
    std::copy(pattern.GetJ(), pattern.GetJ() + nnz, J);
    auto *data = new mfem::real_t[static_cast<std::size_t>(nnz)]();
    auto A = std::make_unique<mfem::SparseMatrix>(I, J, data, n, n, true, true, false);
    A->SortColumnIndices();
    return A;
}

void DofSystem::add_element_vector(mfem::Vector &b, int e, const ElementBuilder::Vector &v) const
{
    const int *local = nullptr;
    const double *sign = nullptr;
    unsigned int n = 0;
    dofs_.element_dofs(e, local, sign, n);
    for (unsigned int i = 0; i < n; ++i)
        b(local[i]) += sign[i] * v(i);
}

void DofSystem::add_element_vector_atomic(mfem::Vector &b, int e, const ElementBuilder::Vector &v) const
{
    const int *local = nullptr;
    const double *sign = nullptr;
    unsigned int n = 0;
    dofs_.element_dofs(e, local, sign, n);
    double *data = b.GetData();
    for (unsigned int i = 0; i < n; ++i)
    {
        const double value = sign[i] * v(i);
#pragma omp atomic
        data[local[i]] += value;
    }
}

void DofSystem::assemble_operators()
{
    StepReport &report = es_.step_report();
    if (operators_assembled_)
    {
        pytec3d::logging::info("Reusing the operators of the model");
        return;
    }
    report.element_matrices_computed = true;
    const pytec3d::fem::FEMesh &mesh = es_.mesh();
    pytec3d::logging::ScopedTimer timer("Operators assembled (" + std::to_string(mesh.n_local_elements()) + " local elements)",
                                        &report.seconds_element_matrices);
    // K, C, M on the sparsity of the element-DOF connectivity, computed once.
    pattern_ = connectivity(nullptr);
    const std::unique_ptr<mfem::SparseMatrix> K_local = matrix_with_pattern(pattern_);
    std::unique_ptr<mfem::SparseMatrix> C_local, M_local;
    mfem::SparseMatrix &K = *K_local;

    // The field values of the elements, element by element.
    const std::array<const pytec3d::fem::SystemLayout *, 2> layouts = {&es_.nodal_layout(), &es_.element_layout()};
    for (std::size_t k = 0; k < 2; ++k)
    {
        field_offsets_[k].assign(1, 0);
        for (int e = 0; e < mesh.n_local_elements(); ++e)
            field_offsets_[k].push_back(field_offsets_[k].back() + static_cast<int>(layouts[k]->n_element_dofs(e)));
        for (std::size_t order = 0; order < 3; ++order)
        {
            field_blocks_[k][order].clear();
            block_offsets_[k][order].assign(1, 0);
        }
        evaluated_field_elements_[k].clear();
    }
    law_elements_.clear();

    bool local_C = false, local_M = false;
    // The element matrices of a chunk, computed in parallel, then added in element order.
    struct Built
    {
        std::vector<ElementBuilder::Matrix> dof;
        std::array<std::vector<ElementBuilder::Matrix>, 2> fields;
        bool laws = false;
        std::array<bool, 2> evaluated = {false, false};
    };
    std::vector<Built> built;
    for (int begin = 0; begin < mesh.n_local_elements(); begin += kChunk)
    {
        const int end = std::min(mesh.n_local_elements(), begin + kChunk);
        built.assign(static_cast<std::size_t>(end - begin), Built{});
        parallel_for(es_, nullptr, begin, end, [&](ElementBuilder &builder, int e)
                     {
                         builder.set_element(e);
                         builder.build_element();
                         Built &out = built[static_cast<std::size_t>(e - begin)];
                         out.dof = std::move(builder.dof_matrices());
                         out.fields[0] = std::move(builder.field_matrices(pytec3d::FieldKind::Nodal));
                         out.fields[1] = std::move(builder.field_matrices(pytec3d::FieldKind::Element));
                         out.laws = builder.nonlinear_element();
                         out.evaluated = {builder.has_evaluated_fields(pytec3d::FieldKind::Nodal),
                                          builder.has_evaluated_fields(pytec3d::FieldKind::Element)}; });
        for (int e = begin; e < end; ++e)
        {
            const Built &out = built[static_cast<std::size_t>(e - begin)];
            const std::vector<ElementBuilder::Matrix> &m = out.dof;
            add_element_matrix(K, e, m[0]);
            if (m[1].m() > 0)
            {
                if (!C_local)
                    C_local = matrix_with_pattern(pattern_);
                add_element_matrix(*C_local, e, m[1]);
                local_C = true;
            }
            if (m.size() > 2 && m[2].m() > 0)
            {
                if (!M_local)
                    M_local = matrix_with_pattern(pattern_);
                add_element_matrix(*M_local, e, m[2]);
                local_M = true;
            }
            if (out.laws)
                law_elements_.push_back(e);

            for (std::size_t kind = 0; kind < 2; ++kind)
            {
                const std::vector<ElementBuilder::Matrix> &F = out.fields[kind];
                for (std::size_t order = 0; order < 3; ++order)
                {
                    std::vector<double> &blocks = field_blocks_[kind][order];
                    if (order < F.size() && F[order].m() > 0)
                        blocks.insert(blocks.end(), F[order].get_values().begin(), F[order].get_values().end());
                    block_offsets_[kind][order].push_back(blocks.size());
                }
                if (out.evaluated[kind])
                    evaluated_field_elements_[kind].push_back(e);
            }
        }
    }
    for (auto &kind : field_blocks_)
        for (auto &blocks : kind)
            blocks.shrink_to_fit();

    int flags[2] = {local_C ? 1 : 0, local_M ? 1 : 0};
    MPI_Allreduce(MPI_IN_PLACE, flags, 2, MPI_INT, MPI_MAX, mesh.comm());
    has_C_ = flags[0] != 0;
    has_M_ = flags[1] != 0;
    // The sizes of what is kept: the local operators (here) and the field operators (until the end of the session).
    long long field_entries = 0;
    for (const auto &kind : field_blocks_)
        for (const auto &blocks : kind)
            field_entries += static_cast<long long>(blocks.size());
    pytec3d::logging::info("Operators: K {} local entries; field operators {} values ({:.1f} MB on this rank).",
                           K.NumNonZeroElems(), field_entries, 8.0 * static_cast<double>(field_entries) / 1048576.0);
    // Collective: every rank assembles C (M) when one has it.
    if (has_C_ && !C_local)
        C_local = matrix_with_pattern(pattern_);
    if (has_M_ && !M_local)
        M_local = matrix_with_pattern(pattern_);
    K_ = dofs_.assemble(K);
    if (has_C_)
        C_ = dofs_.assemble(*C_local);
    if (has_M_)
        M_ = dofs_.assemble(*M_local);
    law_pattern_ = connectivity(&law_elements_);
    operators_assembled_ = true;
}

void DofSystem::update_effective_matrix()
{
    const TimeScheme scheme = es_.time_scheme();
    // A new K_eff object: the reduced system recomputes P^T K_eff P (the old one is freed first).
    Keff_owned_.reset();
    const bool c = has_C_ && scheme.transient() && scheme.c_v != 0.0;
    const bool m = has_M_ && scheme.transient() && scheme.c_a != 0.0;
    if (!c && !m)
    {
        Keff_ = K_.get(); // A static step: K itself, no copy.
        return;
    }
    std::unique_ptr<mfem::HypreParMatrix> keff;
    if (c)
        keff = add(1.0, *K_, scheme.c_v, *C_);
    if (m)
        keff = add(1.0, keff ? *keff : *K_, scheme.c_a, *M_);
    Keff_owned_ = std::move(keff);
    Keff_ = Keff_owned_.get();
}

void DofSystem::add_nodal_loads(mfem::Vector &b, double scale) const
{
    // Nodal loads are point values on the DOFs, not integrals, added once by the rank that owns the DOF. A scalar load
    // goes to the DOF of the point (an edge field: the DOF of the midside point, in its orientation). A vector load of a
    // vector field is given in the node frame (node_angles, like the Dirichlet values _F): f = R f_node in global.
    const pytec3d::ModelArrays &arrays = es_.model().arrays();
    const auto &rotations = es_.node_rotations();
    const auto add = [&](long long l, double value)
    {
        if (l >= 0 && dofs_.owned(static_cast<std::size_t>(l)) && value != 0.0)
            b(static_cast<int>(dofs_.global(static_cast<std::size_t>(l)) - dofs_.first_owned())) += scale * value;
    };
    for (unsigned int v = 0; v < dofs_.n_vars(); ++v)
    {
        const pytec3d::fem::Variable &variable = dofs_.variable(v);
        const std::string_view load = pytec3d::registry::info(variable.id).nodal_load;
        if (load.empty())
            continue;
        const pytec3d::RealArray *array = arrays.point_field(std::string(load));
        if (!array)
            continue;
        for (const auto &[point, dof] : dofs_.point_dofs(v))
        {
            (void)dof;
            const auto p = static_cast<std::size_t>(point);
            if (variable.vector_lagrange())
            {
                pytec3d::la::RealVectorValue f((*array)(p, 0), (*array)(p, 1), (*array)(p, 2));
                if (!rotations.empty())
                    f = rotations[p] * f;
                for (unsigned int c = 0; c < static_cast<unsigned int>(variable.vdim); ++c)
                    add(dofs_.point_dof(v, point, c), f(c));
                continue;
            }
            add(dofs_.point_dof(v, point, 0), array->n_components() == 1 ? (*array)(p) : (*array)(p, variable.components[0]));
        }
    }
}

void DofSystem::assemble_rhs(mfem::Vector &F) const
{
    const pytec3d::fem::FEMesh &mesh = es_.mesh();
    mfem::Vector b(static_cast<int>(dofs_.n_local()));
    b = 0.0;
    const LocalVectors vectors = local_vectors();
    parallel_for(es_, &vectors, 0, mesh.n_local_elements(), [&](ElementBuilder &builder, int e)
                 {
                     builder.set_element(e);
                     add_element_vector_atomic(b, e, builder.build_element_rhs()); });
    F.SetSize(static_cast<int>(dofs_.n_owned()));
    dofs_.assemble(b, F);
    add_nodal_loads(F, 1.0);

    // The history of the step: - C h_v - M h_a, h = hu u_n + hv v_n + ha a_n.
    const TimeScheme scheme = es_.time_scheme();
    if (!scheme.transient())
        return;
    const auto history = [&](const mfem::HypreParMatrix &A, Real hu, Real hv, Real ha)
    {
        if (hu == 0.0 && hv == 0.0 && ha == 0.0)
            return;
        mfem::Vector h(vectors_[OLD_SOLUTION]);
        h *= hu;
        h.Add(hv, vectors_[OLD_VELOCITY]);
        h.Add(ha, vectors_[OLD_ACCELERATION]);
        A.Mult(-1.0, h, 1.0, F);
    };
    if (has_C_)
        history(*C_, scheme.hu_v, scheme.hv_v, scheme.ha_v);
    if (has_M_)
        history(*M_, scheme.hu_a, scheme.hv_a, scheme.ha_a);
}

// ---- solve -----------------------------------------------------------------------------------------------------

void DofSystem::solve()
{
    if (nonlinear_)
        solve_nonlinear();
    else
        solve_linear();
    first_solve_ = false;
}

void DofSystem::solve_linear()
{
    StepReport &report = es_.step_report();
    pytec3d::logging::ScopedTimer timer("Linear system assembled and solved", &report.seconds_solve);
    const bool update_matrix = check_update_system_matrix();
    report.system_matrix_updated |= update_matrix;
    if (update_matrix)
        update_effective_matrix();
    pytec3d::logging::info("Building the right-hand side");
    assemble_rhs(F_);

    // (P^T K P) w = P^T (F - K g); u = P w + g.
    reduced_->set_matrix(*Keff_, update_matrix);
    mfem::Vector f(F_);
    Keff_->Mult(-1.0, reduced_->offset(), 1.0, f);
    mfem::Vector fhat(reduced_->n_free_local()), w(reduced_->n_free_local());
    reduced_->project(f, fhat);
    w = 0.0;
    if (settings_.linear.initial_guess_nonzero)
        reduced_->restrict_to_free(vectors_[SOLUTION], w);
    const ReducedSystem::Counts before = reduced_->counts();
    const auto [iterations, residual] = reduced_->solve(fhat, w, settings_.linear.relative_tolerance);
    report.full_factorizations += reduced_->counts().full - before.full;
    report.numeric_factorizations += reduced_->counts().numeric - before.numeric;
    report.true_residual = reduced_->relative_true_residual(w, fhat);
    report.factorization = pytec3d::diagnostics::factorization_report(reduced_->ksp());
    reduced_->expand(w, vectors_[SOLUTION]);
    invalidate(SOLUTION);
    update_time_derivatives();

    const KSPConvergedReason reason = reduced_->converged_reason();
    report.linear_iterations = iterations;
    report.linear_residual = residual;
    report.linear_reason = KSPConvergedReasons[reason];
    report.unconstrained_dofs = static_cast<std::size_t>(reduced_->n_free());
    pytec3d::logging::info("Linear solve: {} iteration(s), residual {:.3e}, true residual {:.3e}, {} ({} unconstrained "
                           "DOFs on this rank).",
                           iterations, residual, report.true_residual, report.linear_reason, reduced_->n_free_local());
    pytec3d::diagnostics::check_linear_solve(report.factorization, report.true_residual);
    if (reason < 0)
        pytec3d::logging::warn("The linear solver did not converge: {}.", report.linear_reason);
}

void DofSystem::update_time_derivatives()
{
    // v = c_v u + hu_v u_n + hv_v v_n + ha_v a_n, a likewise (time_step.h); zero on a static step.
    const TimeScheme scheme = es_.time_scheme();
    const auto update = [&](mfem::Vector &x, Real c, Real hu, Real hv, Real ha)
    {
        x = 0.0;
        if (!scheme.transient())
            return;
        x.Add(c, vectors_[SOLUTION]);
        x.Add(hu, vectors_[OLD_SOLUTION]);
        x.Add(hv, vectors_[OLD_VELOCITY]);
        x.Add(ha, vectors_[OLD_ACCELERATION]);
    };
    update(vectors_[VELOCITY], scheme.c_v, scheme.hu_v, scheme.hv_v, scheme.ha_v);
    update(vectors_[ACCELERATION], scheme.c_a, scheme.hu_a, scheme.hv_a, scheme.ha_a);
    invalidate(VELOCITY);
    invalidate(ACCELERATION);
}

// ---- nonlinear ---------------------------------------------------------------------------------------------------

void DofSystem::assemble_laws(bool jacobian, mfem::Vector *residual)
{
    const TimeScheme scheme = es_.time_scheme();
    const int n_local = static_cast<int>(dofs_.n_local());
    const std::unique_ptr<mfem::SparseMatrix> J_local = matrix_with_pattern(law_pattern_);
    mfem::SparseMatrix &J = *J_local;
    mfem::Vector r(n_local);
    r = 0.0;
    const LocalVectors vectors = local_vectors();
    const bool picard = settings_.nonlinear.picard;
    std::vector<ElementBuilder::Matrix> tangents;
    const int n_laws = static_cast<int>(law_elements_.size());
    for (int begin = 0; begin < n_laws; begin += kChunk)
    {
        const int end = std::min(n_laws, begin + kChunk);
        tangents.assign(static_cast<std::size_t>(end - begin), ElementBuilder::Matrix());
        parallel_for(es_, &vectors, begin, end, [&](ElementBuilder &builder, int i)
                     {
                         const int e = law_elements_[static_cast<std::size_t>(i)];
                         builder.set_element(e);
                         ElementBuilder::Vector law_residual;
                         builder.build_element_tangent_stiffness_residual(scheme, picard, tangents[static_cast<std::size_t>(i - begin)],
                                                                          law_residual);
                         if (!residual)
                             return;
                         // The stored rates of the law outputs under dt(): hv_d q'_n + ha_d q''_n (the builder adds the rest).
                         if (scheme.transient())
                             if (const auto it = old_law_rates_.find(e); it != old_law_rates_.end())
                             {
                                 const std::array<std::array<Real, 2>, 2> h = {{{scheme.hv_v, scheme.ha_v}, {scheme.hv_a, scheme.ha_a}}};
                                 for (unsigned int d = 0; d < 2; ++d)
                                 {
                                     if (it->second.v[d].size() == law_residual.size())
                                         law_residual.add(h[d][0], it->second.v[d]);
                                     if (it->second.a[d].size() == law_residual.size())
                                         law_residual.add(h[d][1], it->second.a[d]);
                                 }
                             }
                         add_element_vector_atomic(r, e, law_residual); });
        if (jacobian)
            for (int i = begin; i < end; ++i)
                add_element_matrix(J, law_elements_[static_cast<std::size_t>(i)], tangents[static_cast<std::size_t>(i - begin)]);
    }
    if (residual)
    {
        residual->SetSize(static_cast<int>(dofs_.n_owned()));
        dofs_.assemble(r, *residual);
    }
    if (jacobian)
    {
        const std::unique_ptr<mfem::HypreParMatrix> laws = dofs_.assemble(J);
        jacobian_ = add(1.0, *Keff_, 1.0, *laws);
    }
}

double DofSystem::assemble_residual()
{
    // R = K_eff u - F + R_laws.
    mfem::Vector laws;
    assemble_laws(false, &laws);
    residual_.SetSize(static_cast<int>(dofs_.n_owned()));
    residual_ = 0.0; // hypre forms beta * y with beta = 0
    Keff_->Mult(vectors_[SOLUTION], residual_);
    residual_ -= F_;
    residual_ += laws;
    reduced_residual_.SetSize(reduced_->n_free_local());
    reduced_->project(residual_, reduced_residual_);
    return norm(reduced_residual_);
}

void DofSystem::record_newton_iterate(double residual)
{
    std::vector<NewtonIteration> &history = es_.step_report().newton_history;
    NewtonIteration iteration;
    iteration.residual = residual;
    if (!has_previous_iterate_)
    {
        previous_iterate_ = vectors_[SOLUTION];
        has_previous_iterate_ = true;
    }
    else
    {
        // The step from the previous iterate, line search included.
        previous_iterate_ -= vectors_[SOLUTION];
        iteration.step = norm(previous_iterate_);
        const double n = norm(vectors_[SOLUTION]);
        iteration.relative_step = n > 0.0 ? iteration.step / n : iteration.step;
        if (full_newton_step_ > 0.0)
            iteration.step_length = iteration.step / full_newton_step_;
        previous_iterate_ = vectors_[SOLUTION];
    }
    full_newton_step_ = -1.0;
    history.push_back(iteration);
}

bool DofSystem::test_convergence(double residual, double step, bool linear_solve_finished)
{
    const NewtonSettings &newton = settings_.nonlinear;
    bool converged = false;
    if (residual < newton.absolute_residual_tolerance)
    {
        solve_result_ |= CONVERGED_ABSOLUTE_RESIDUAL;
        converged = true;
    }
    if (residual / max_residual_norm_ < newton.relative_residual_tolerance)
    {
        solve_result_ |= CONVERGED_RELATIVE_RESIDUAL;
        converged = true;
    }
    // For incomplete linear solves, it is not safe to test step sizes.
    if (!linear_solve_finished)
        return converged;
    if (step < newton.absolute_step_tolerance)
    {
        solve_result_ |= CONVERGED_ABSOLUTE_STEP;
        converged = true;
    }
    if (max_solution_norm_ != 0.0 && step / max_solution_norm_ < newton.relative_step_tolerance)
    {
        solve_result_ |= CONVERGED_RELATIVE_STEP;
        converged = true;
    }
    return converged;
}

double DofSystem::line_search(double tol, double last_residual, double &current_residual, const mfem::Vector &du)
{
    // libMesh's NewtonSolver::line_search(): u holds u_prev - du; the step shrinks to bx du, u = u_prev - bx du.
    const NewtonSettings &newton = settings_.nonlinear;
    if (current_residual < last_residual ||
        (!newton.require_residual_reduction && (!newton.require_finite_residual || !std::isnan(current_residual))))
        return 1.0;

    mfem::Vector &u = vectors_[SOLUTION];
    const auto evaluate = [&]()
    {
        invalidate(SOLUTION);
        current_residual = assemble_residual();
        if (newton.verbose)
            pytec3d::logging::info("  Current Residual: {}", current_residual);
        return current_residual;
    };

    double ax = 0.0; // First abscissa: no negative steps.
    double cx = 1.0; // Second abscissa: no extrapolated steps.
    double bx = 1.0;
    while (std::isnan(current_residual) || (current_residual > last_residual && newton.require_residual_reduction))
    {
        double substep = 0.5;
        if (newton.brent_line_search && !std::isnan(current_residual))
            substep = std::max(std::min(0.5, last_residual / current_residual), tol * 2.0);
        u.Add(bx * (1.0 - substep), du);
        bx *= substep;
        if (newton.verbose)
            pytec3d::logging::info("  Shrinking Newton step to {}", bx);
        evaluate();
        if (bx / 2.0 < newton.minimum_step_length && (std::isnan(current_residual) || current_residual > last_residual))
        {
            pytec3d::logging::info("Inexact Newton step FAILED at step {}", outer_iterations_);
            if (!newton.continue_after_backtrack_failure)
                pytec3d::logging::error_and_abort("The Newton solver failed: the line search could not reduce the residual "
                                                  "(continue_after_backtrack_failure is off).");
            pytec3d::logging::info("Continuing anyway ...");
            solve_result_ = DIVERGED_BACKTRACKING_FAILURE;
            return bx;
        }
    }
    if (!newton.brent_line_search)
        return bx;

    // Brent's method on the residual norm along the step.
    double e = 0.0;
    double x = bx, w = bx, v = bx;
    double fx = current_residual, fw = current_residual, fv = current_residual;
    const unsigned int max_i = 20;
    const double golden_ratio = 1.0 - (std::sqrt(5.0) - 1.0) / 2.0;
    const auto sign = [](double a, double b) { return b >= 0.0 ? std::abs(a) : -std::abs(a); };
    for (unsigned int i = 1; i <= max_i; ++i)
    {
        const double xm = (ax + cx) * 0.5;
        const double tol1 = tol * std::abs(x) + tol * tol;
        const double tol2 = 2.0 * tol1;
        if (std::abs(x - xm) <= (tol2 - 0.5 * (cx - ax)))
            return x;
        double d = 0.0;
        if (std::abs(e) > tol1)
        {
            const double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            else
                q = std::abs(q);
            if (std::abs(p) >= std::abs(0.5 * q * e) || p <= q * (ax - x) || p >= q * (cx - x))
            {
                e = x >= xm ? ax - x : cx - x;
                d = golden_ratio * e;
            }
            else
            {
                d = p / q;
                if (x + d - ax < tol2 || cx - (x + d) < tol2)
                    d = sign(tol1, xm - x);
            }
        }
        else
        {
            e = x >= xm ? ax - x : cx - x;
            d = golden_ratio * e;
        }
        const double uu = std::abs(d) >= tol1 ? x + d : x + sign(tol1, d);
        u.Add(bx - uu, du);
        bx = uu;
        if (newton.verbose)
            pytec3d::logging::info("  Shrinking Newton step to {}", bx);
        const double fu = evaluate();
        if (fu <= fx)
        {
            if (uu >= x)
                ax = x;
            else
                cx = x;
            v = w;
            w = x;
            x = uu;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if (uu < x)
                ax = uu;
            else
                cx = uu;
            if (fu <= fw || w == x)
            {
                v = w;
                w = uu;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = uu;
                fv = fu;
            }
        }
    }
    pytec3d::logging::warn("Too many iterations used in the Brent line search.");
    return bx;
}

void DofSystem::solve_nonlinear()
{
    StepReport &report = es_.step_report();
    pytec3d::logging::ScopedTimer timer("Nonlinear system solved", &report.seconds_solve);
    const NewtonSettings &newton = settings_.nonlinear;
    has_previous_iterate_ = false;
    full_newton_step_ = -1.0;
    solve_result_ = INVALID_SOLVE_RESULT;
    inner_iterations_ = 0;
    const ReducedSystem::Counts before = reduced_->counts();

    const bool update_matrix = check_update_system_matrix();
    report.system_matrix_updated |= update_matrix;
    if (update_matrix)
        update_effective_matrix();
    assemble_rhs(F_);

    // The constraints exactly: u = P w + g, w the free DOFs of u.
    mfem::Vector w(reduced_->n_free_local());
    reduced_->restrict_to_free(vectors_[SOLUTION], w);
    reduced_->expand(w, vectors_[SOLUTION]);
    invalidate(SOLUTION);
    mfem::Vector delta(reduced_->n_free_local()), du(static_cast<int>(dofs_.n_owned()));

    for (outer_iterations_ = 0; outer_iterations_ < newton.max_iterations; ++outer_iterations_)
    {
        if (newton.verbose)
            pytec3d::logging::info("Assembling the System");
        // R and J at the iterate: the law terms once for both.
        mfem::Vector laws;
        assemble_laws(true, &laws);
        residual_.SetSize(static_cast<int>(dofs_.n_owned()));
        residual_ = 0.0; // hypre forms beta * y with beta = 0
        Keff_->Mult(vectors_[SOLUTION], residual_);
        residual_ -= F_;
        residual_ += laws;
        reduced_residual_.SetSize(reduced_->n_free_local());
        reduced_->project(residual_, reduced_residual_);
        double current_residual = norm(reduced_residual_);
        record_newton_iterate(current_residual);
        if (std::isnan(current_residual))
            pytec3d::logging::error_and_abort("Nonlinear solver DIVERGED at step " + std::to_string(outer_iterations_) +
                                              " with residual Not-a-Number.");
        if (current_residual <= newton.absolute_residual_tolerance)
        {
            if (newton.verbose)
                pytec3d::logging::info("Linear solve unnecessary; residual {} meets absolute tolerance {}", current_residual,
                                       newton.absolute_residual_tolerance);
            solve_result_ |= CONVERGED_ABSOLUTE_RESIDUAL;
            if (current_residual == 0.0)
            {
                if (newton.relative_residual_tolerance > 0)
                    solve_result_ |= CONVERGED_RELATIVE_RESIDUAL;
                if (newton.absolute_step_tolerance > 0)
                    solve_result_ |= CONVERGED_ABSOLUTE_STEP;
                if (newton.relative_step_tolerance > 0)
                    solve_result_ |= CONVERGED_RELATIVE_STEP;
            }
            break;
        }
        const double last_residual = current_residual;
        max_residual_norm_ = std::max(current_residual, max_residual_norm_);
        double norm_total = norm(vectors_[SOLUTION]);
        max_solution_norm_ = std::max(max_solution_norm_, norm_total);
        if (newton.verbose)
            pytec3d::logging::info("Nonlinear Residual: {}", current_residual);

        // The Newton step on the free DOFs: (P^T J P) delta = P^T R, du = P delta.
        reduced_->set_matrix(*jacobian_, true);
        delta = 0.0;
        const auto [iterations, linear_residual] =
            reduced_->solve(reduced_residual_, delta, settings_.linear.relative_tolerance);
        (void)linear_residual;
        reduced_->prolong(delta, du);
        inner_iterations_ += iterations;
        const bool linear_solve_finished = iterations != kMaxLinearIterations;
        double norm_delta = norm(du);
        full_newton_step_ = norm_delta;
        if (newton.verbose)
            pytec3d::logging::info("Linear solve finished, step {}; trying full Newton step", iterations);
        vectors_[SOLUTION] -= du;
        invalidate(SOLUTION);

        if (newton.require_residual_reduction || newton.require_finite_residual ||
            outer_iterations_ + 1 < newton.max_iterations || !newton.continue_after_max_iterations)
        {
            current_residual = assemble_residual();
            if (newton.verbose)
                pytec3d::logging::info("  Current Residual: {}", current_residual);
            if (test_convergence(current_residual, norm_delta, linear_solve_finished && current_residual <= last_residual))
            {
                ++outer_iterations_;
                break;
            }
        }

        const double step_length = line_search(std::sqrt(kTolerance), last_residual, current_residual, du);
        norm_delta *= step_length;
        if (solve_result_ == DIVERGED_BACKTRACKING_FAILURE)
        {
            ++outer_iterations_;
            break;
        }
        if (outer_iterations_ + 1 >= newton.max_iterations)
        {
            pytec3d::logging::info("  Nonlinear solver reached maximum step {}, latest evaluated residual {}",
                                   newton.max_iterations, current_residual);
            if (!newton.continue_after_max_iterations)
                pytec3d::logging::error_and_abort("The Newton solver reached its maximum iterations "
                                                  "(continue_after_max_iterations is off).");
            solve_result_ = DIVERGED_MAX_NONLINEAR_ITERATIONS;
            pytec3d::logging::info("  Continuing...");
            continue;
        }
        norm_total = norm(vectors_[SOLUTION]);
        max_solution_norm_ = std::max(max_solution_norm_, norm_total);
        if (newton.verbose)
            pytec3d::logging::info("  Nonlinear step: |du|/|u| = {}, |du| = {}", norm_delta / norm_total, norm_delta);
        if (test_convergence(current_residual, norm_delta / step_length, linear_solve_finished))
        {
            ++outer_iterations_;
            break;
        }
    }

    report.full_factorizations += reduced_->counts().full - before.full;
    report.numeric_factorizations += reduced_->counts().numeric - before.numeric;
    // The last iterate: the residual of the last step, without the Jacobian.
    if (has_previous_iterate_)
    {
        record_newton_iterate(norm(reduced_residual_));
        if (report.newton_history.size() >= 2 && report.newton_history.back().step == 0.0)
            report.newton_history.pop_back(); // Converged at the start of an iteration: no new iterate.
        has_previous_iterate_ = false;
    }
    update_time_derivatives();
    update_law_rates();
    report.factorization = pytec3d::diagnostics::factorization_report(reduced_->ksp());
    pytec3d::diagnostics::check_linear_solve(report.factorization, -1.0);

    std::string reasons;
    for (const auto &[flag, text] : std::initializer_list<std::pair<unsigned int, const char *>>{
             {CONVERGED_ABSOLUTE_RESIDUAL, "absolute residual"},
             {CONVERGED_RELATIVE_RESIDUAL, "relative residual"},
             {CONVERGED_ABSOLUTE_STEP, "absolute step"},
             {CONVERGED_RELATIVE_STEP, "relative step"},
             {DIVERGED_MAX_NONLINEAR_ITERATIONS, "maximum iterations reached"},
             {DIVERGED_BACKTRACKING_FAILURE, "line search failed"},
             {DIVERGED_LINEAR_SOLVER_FAILURE, "linear solver failed"},
             {DIVERGED_NO_REASON, "diverged"}})
        if (solve_result_ & flag)
            reasons += (reasons.empty() ? "" : ", ") + std::string(text);
    report.newton_iterations = outer_iterations_;
    report.linear_iterations = inner_iterations_;
    report.unconstrained_dofs = static_cast<std::size_t>(reduced_->n_free());
    report.newton_result = reasons.empty() ? "no convergence test met" : reasons;
    pytec3d::logging::info("Newton: {} iteration(s), {} linear iteration(s); {}.", report.newton_iterations,
                           report.linear_iterations, report.newton_result);
    if (solve_result_ >= DIVERGED_NO_REASON)
        pytec3d::logging::warn("The Newton solver did not converge ({}); the solution is the last iterate.", reasons);
}

void DofSystem::update_law_rates()
{
    // q' = c_v q(u) + hu_v q(u_n) + hv_v q'_n + ha_v q''_n and q'' likewise (as the DOFs): the rates of the law outputs
    // under dt(), per element; zero on a static step.
    const TimeScheme scheme = es_.time_scheme();
    ElementBuilder builder(es_);
    law_rates_.clear();
    for (const int e : law_elements_)
    {
        builder.set_element(e);
        if (!builder.rate_residual_element())
            continue;
        LawRates rates;
        const auto old = old_law_rates_.find(e);
        for (unsigned int d = 0; d < 2; ++d)
        {
            const ElementBuilder::Vector q = builder.law_output(local(SOLUTION), d + 1);
            const ElementBuilder::Vector q_n = builder.law_output(local(OLD_SOLUTION), d + 1);
            const auto rate = [&](Real c, Real hu, Real hv, Real ha)
            {
                ElementBuilder::Vector x(q.size());
                if (!scheme.transient())
                    return x;
                x.add(c, q);
                x.add(hu, q_n);
                if (old != old_law_rates_.end() && old->second.v[d].size() == x.size())
                {
                    x.add(hv, old->second.v[d]);
                    x.add(ha, old->second.a[d]);
                }
                return x;
            };
            rates.v[d] = rate(scheme.c_v, scheme.hu_v, scheme.hv_v, scheme.ha_v);
            rates.a[d] = rate(scheme.c_a, scheme.hu_a, scheme.hv_a, scheme.ha_a);
        }
        law_rates_.emplace(e, std::move(rates));
    }
}

// ---- fields ------------------------------------------------------------------------------------------------------

const pytec3d::fem::QuadratureData &DofSystem::point_quadrature() const
{
    if (!point_quadrature_)
    {
        std::vector<const mfem::IntegrationRule *> rules;
        for (int e = 0; e < es_.mesh().n_local_elements(); ++e)
            rules.push_back(&es_.mesh().point_rule(e));
        point_quadrature_ = std::make_unique<pytec3d::fem::QuadratureData>(es_.mesh(), rules);
    }
    return *point_quadrature_;
}

void DofSystem::accumulate_dof_fields(std::vector<std::vector<double>> &sums, std::vector<std::vector<double>> &counts) const
{
    const pytec3d::fem::FEMesh &mesh = es_.mesh();
    const std::vector<double> &u = local(SOLUTION);
    bool edges = false;
    for (unsigned int v = 0; v < dofs_.n_vars(); ++v)
        edges = edges || dofs_.variable(v).edge();
    const pytec3d::fem::QuadratureData *points_data = edges ? &point_quadrature() : nullptr;
#pragma omp parallel
    {
        std::vector<double> values;
#pragma omp for schedule(dynamic, 64)
        for (int e = 0; e < mesh.n_local_elements(); ++e)
        {
            dofs_.element_values(e, u, values);
            std::size_t offset = 0;
            for (unsigned int v = 0; v < dofs_.n_vars(); ++v)
            {
                const unsigned int n = dofs_.n_element_dofs(v, e);
                if (n == 0)
                    continue;
                accumulate_at_points(mesh, es_.catalog(), es_.reference_tables(), points_data, e, dofs_.variable(v),
                                     values.data() + offset, sums[v], counts[v]);
                offset += n;
            }
        }
    }
}

void DofSystem::accumulate_fields(pytec3d::FieldKind kind, std::vector<std::vector<double>> &sums,
                                  std::vector<std::vector<double>> &counts) const
{
    const std::size_t k = kind == pytec3d::FieldKind::Nodal ? 0 : 1;
    const pytec3d::fem::FEMesh &mesh = es_.mesh();
    const pytec3d::fem::SystemLayout &layout = k == 0 ? es_.nodal_layout() : es_.element_layout();
    StepReport &report = es_.step_report();
    pytec3d::logging::ScopedTimer timer(k == 0 ? "Nodal fields computed" : "Element fields computed",
                                        k == 0 ? &report.seconds_nodal_fields : &report.seconds_element_fields);

    // The static operators applied to the solution and its rates: F0 u + F1 du/dt + F2 d2u/dt2.
    mfem::Vector values(field_offsets_[k].back());
    values = 0.0;
    const std::array<VectorId, 3> rates = {SOLUTION, VELOCITY, ACCELERATION};
    for (std::size_t order = 0; order < 3; ++order)
    {
        const std::vector<double> &blocks = field_blocks_[k][order];
        if (blocks.empty())
            continue;
        const std::vector<double> &l = local(rates[order]);
        const std::vector<std::size_t> &offsets = block_offsets_[k][order];
#pragma omp parallel for schedule(static)
        for (int e = 0; e < mesh.n_local_elements(); ++e)
        {
            const std::size_t begin = offsets[static_cast<std::size_t>(e)];
            if (offsets[static_cast<std::size_t>(e) + 1] == begin)
                continue;
            const int *local_dofs = nullptr;
            const double *sign = nullptr;
            unsigned int n = 0;
            dofs_.element_dofs(e, local_dofs, sign, n);
            const int row0 = field_offsets_[k][static_cast<std::size_t>(e)];
            const int rows = field_offsets_[k][static_cast<std::size_t>(e) + 1] - row0;
            const double *F = blocks.data() + begin;
            for (int r = 0; r < rows; ++r)
            {
                double sum = 0.0;
                for (unsigned int j = 0; j < n; ++j)
                    sum += F[static_cast<std::size_t>(r) * n + j] * sign[j] * l[static_cast<std::size_t>(local_dofs[j])];
                values(row0 + r) += sum;
            }
        }
    }
    // The terms evaluated at every output: dynamic fields (laws), loads alone (e.g. Jt = Js).
    // On the threads: each element updates only its own block of values.
    const std::vector<int> &evaluated = evaluated_field_elements_[k];
    const LocalVectors vectors = local_vectors();
    parallel_for(es_, &vectors, 0, static_cast<int>(evaluated.size()), [&](ElementBuilder &builder, int j)
                 {
                     const int e = evaluated[static_cast<std::size_t>(j)];
                     builder.set_element(e);
                     const int row0 = field_offsets_[k][static_cast<std::size_t>(e)];
                     const int n = field_offsets_[k][static_cast<std::size_t>(e) + 1] - row0;
                     ElementBuilder::Vector block(static_cast<unsigned int>(n));
                     for (int i = 0; i < n; ++i)
                         block(static_cast<unsigned int>(i)) = values(row0 + i);
                     builder.add_evaluated_fields(kind, layout, block);
                     for (int i = 0; i < n; ++i)
                         values(row0 + i) = block(static_cast<unsigned int>(i)); });

#pragma omp parallel for schedule(dynamic, 64)
    for (int e = 0; e < mesh.n_local_elements(); ++e)
    {
        std::size_t offset = static_cast<std::size_t>(field_offsets_[k][static_cast<std::size_t>(e)]);
        for (unsigned int v = 0; v < layout.n_vars(); ++v)
        {
            const unsigned int n = layout.n_element_dofs(v, e);
            if (n == 0)
                continue;
            const pytec3d::fem::Variable &var = layout.variable(v);
            const double *coefficients = values.GetData() + offset;
            if (var.monomial())
                accumulate(sums[v], counts[v], mesh.cell(e), var.n_components, coefficients);
            else
                accumulate_at_points(mesh, es_.catalog(), es_.reference_tables(), nullptr, e, var, coefficients, sums[v],
                                     counts[v]);
            offset += n;
        }
    }
}

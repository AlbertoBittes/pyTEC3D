/**
 * @file element_kernels.cpp
 * @brief ElementKernels: kernel-driven element matrices, tangent/residual and post-processed fields.
 */
#include "element_kernels.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

#include <libmesh/elem.h>
#include <libmesh/vector_value.h>

#include "bh_curve.h"
#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using libMesh::DenseMatrix;
        using libMesh::DenseVector;
        using libMesh::Real;
        using libMesh::TensorValue;
        using Vec3 = std::array<Real, 3>;

        const Real kMu0 = 4.0 * libMesh::pi * 1e-7;
        const char *const kBHCoefficient = "nu[][]"; // The coefficient made nonlinear by a BH curve.

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementKernels: " + message);
        }

        struct Coefficient
        {
            bool is_tensor = false;
            bool is_vector = false;
            Real scalar = 0.0;
            TensorValue<Real> tensor;
            Vec3 vector{0.0, 0.0, 0.0};
        };

        // "X[]" but not "X[][]".
        bool is_vector_coefficient(const std::string &name)
        {
            const auto ends_with = [&name](const char *suffix)
            {
                const std::string s(suffix);
                return name.size() >= s.size() && name.compare(name.size() - s.size(), s.size(), s) == 0;
            };
            return ends_with("[]") && !ends_with("[][]");
        }

        // Vector coefficients other than the material magnetization are element loads.
        const char *const kMagnetization = "M[]";

        Real conductivity(Real resistivity)
        {
            return (resistivity > 0.0 && std::isfinite(resistivity)) ? 1.0 / resistivity : 0.0;
        }

        Coefficient evaluate_coefficient(const std::string &name,
                                         const ElementParameters &p,
                                         const ElementTypeInstance &element_type,
                                         const ElementLoads *loads = nullptr)
        {
            const TensorValue<Real> &R = p.local_to_global;
            const auto rotate_diagonal = [&R](Real a, Real b, Real c)
            {
                const TensorValue<Real> diagonal(a, 0.0, 0.0, 0.0, b, 0.0, 0.0, 0.0, c);
                return TensorValue<Real>(R * diagonal * R.transpose());
            };

            Coefficient out;
            if (name == "nu[][]")
            {
                const auto &mur = p.material.relative_permeability;
                out.is_tensor = true;
                out.tensor = rotate_diagonal(1.0 / (kMu0 * mur[0]), 1.0 / (kMu0 * mur[1]), 1.0 / (kMu0 * mur[2]));
            }
            else if (name == "sigma[][]")
            {
                const auto &rho = p.material.resistivity;
                out.is_tensor = true;
                out.tensor = rotate_diagonal(conductivity(rho[0]), conductivity(rho[1]), conductivity(rho[2]));
            }
            else if (name == "cbrt(abs(nu))")
            {
                const auto &mur = p.material.relative_permeability;
                out.scalar = std::cbrt(std::abs(1.0 / (kMu0 * mur[0]) / (kMu0 * mur[1]) / (kMu0 * mur[2])));
            }
            else if (name.rfind("coil_", 0) == 0)
            {
                const RealConstants &c = p.constants;
                const libMesh::RealVectorValue t = R * libMesh::RealVectorValue(c.current_direction[0],
                                                                                c.current_direction[1],
                                                                                c.current_direction[2]);
                const Real turns_density = c.turns / (c.cross_section * c.resistance);
                const Real coil_length = element_type.dimension() == 2
                                             ? (element_type.options().axisymmetric ? 2.0 * libMesh::pi * c.radius : 1.0) *
                                                   c.cross_section * c.thickness
                                             : c.volume;
                const Real kvv = 1.0 / (c.resistance * coil_length);

                if (name == "coil_av[][]" || name == "coil_j[][]")
                {
                    // Column 0 couples a scalar trial to the vector test.
                    const Real scale = name == "coil_j[][]" ? turns_density : -turns_density;
                    out.is_tensor = true;
                    for (unsigned int a = 0; a < 3; ++a)
                        out.tensor(a, 0) = scale * t(a);
                }
                else if (name == "coil_va[][]" || name == "coil_ea[][]")
                {
                    // Row 0 couples a vector trial to the scalar test.
                    const Real scale = name == "coil_va[][]" ? -c.symmetry * turns_density : turns_density;
                    out.is_tensor = true;
                    for (unsigned int b = 0; b < 3; ++b)
                        out.tensor(0, b) = scale * t(b);
                }
                else if (name == "coil_vv")
                    out.scalar = c.symmetry * kvv;
                else if (name == "coil_ee")
                    out.scalar = -kvv;
                else
                    fail("unknown coefficient '" + name + "'.");
            }
            else if (is_vector_coefficient(name))
            {
                std::array<Real, 3> local{0.0, 0.0, 0.0};
                if (name == kMagnetization)
                    local = p.material.magnetization;
                else if (loads)
                {
                    const auto load = loads->element.find(name.substr(0, name.size() - 2));
                    if (load != loads->element.end())
                        local = load->second;
                }
                const libMesh::RealVectorValue global = R * libMesh::RealVectorValue(local[0], local[1], local[2]);
                out.is_vector = true;
                out.vector = {global(0), global(1), global(2)};
            }
            else
                fail("unknown coefficient '" + name + "'.");
            return out;
        }

        // Isotropic H(B) and its Jacobian dH/dB from a BH curve.
        void evaluate_bh(const BHCurve &curve, const Real B[3], unsigned int dim, Real H[3], Real jac[3][3])
        {
            Real b_sq = 0.0;
            for (unsigned int k = 0; k < dim; ++k)
                b_sq += B[k] * B[k];
            const Real b_mag = std::sqrt(b_sq);
            const Real nu = curve.evaluate_reluctivity_from_b(b_mag);
            const Real dhdb = curve.slope_h_vs_b(b_mag);

            for (unsigned int k = 0; k < 3; ++k)
            {
                H[k] = 0.0;
                for (unsigned int l = 0; l < 3; ++l)
                    jac[k][l] = 0.0;
            }
            constexpr Real kZeroB = 1e-14;
            for (unsigned int k = 0; k < dim; ++k)
            {
                H[k] = nu * B[k];
                jac[k][k] += nu;
                if (b_mag > kZeroB)
                    for (unsigned int l = 0; l < dim; ++l)
                        jac[k][l] += (dhdb - nu) * B[k] * B[l] / b_sq;
            }
        }

        // Coefficient times a point value; for nonlinear elements the BH coefficient gives H(B).
        Vec3 apply_coefficient(const std::string &name, const Vec3 &value, const ElementParameters &p, const ElementTypeInstance &et)
        {
            if (name.empty())
                return value;
            if (p.bh_curve && name == kBHCoefficient)
            {
                Real H[3];
                Real jac[3][3];
                evaluate_bh(*p.bh_curve, value.data(), et.dimension(), H, jac);
                return {H[0], H[1], H[2]};
            }

            const Coefficient c = evaluate_coefficient(name, p, et);
            if (c.is_vector)
                fail("vector coefficient '" + name + "' cannot scale a field expression.");
            Vec3 out{0.0, 0.0, 0.0};
            for (unsigned int a = 0; a < 3; ++a)
            {
                if (!c.is_tensor)
                    out[a] = c.scalar * value[a];
                else
                    for (unsigned int b = 0; b < 3; ++b)
                        out[a] += c.tensor(a, b) * value[b];
            }
            return out;
        }

        // Rows of the largest operator: the 3x3 gradient of a 3D vector.
        constexpr unsigned int kMaxRows = 9;

        // integral of test^T * coefficient * trial: the (test DOFs x trial DOFs) block of one kernel.
        DenseMatrix<Real> weak_form_block(const Nabla::Operand &test,
                                          const Coefficient &coefficient,
                                          const Nabla::Operand &trial,
                                          const std::vector<Real> &weights,
                                          Real scale)
        {
            if (test.values.size() != trial.values.size() || test.values.size() != weights.size())
                fail("test, trial and weights must have one entry per quadrature point.");
            if (test.values.empty())
                return DenseMatrix<Real>();
            if (coefficient.is_vector)
                fail("a vector coefficient needs a source kernel without trial.");

            const unsigned int rows_test = test.values[0].m();
            const unsigned int rows_trial = trial.values[0].m();
            if (rows_test > kMaxRows || rows_trial > kMaxRows)
                fail("operators with more than 9 rows are not supported.");
            if (coefficient.is_tensor)
            {
                for (unsigned int c : test.components)
                    if (c > 2)
                        fail("a tensor coefficient needs physical components (x, y, z) on the test side.");
                for (unsigned int c : trial.components)
                    if (c > 2)
                        fail("a tensor coefficient needs physical components (x, y, z) on the trial side.");
            }
            else if (rows_test != rows_trial)
                fail("a scalar coefficient needs test and trial with the same number of rows.");

            DenseMatrix<Real> block(test.values[0].n(), trial.values[0].n());
            for (std::size_t qp = 0; qp < test.values.size(); ++qp)
            {
                const DenseMatrix<Real> &T = test.values[qp];
                const DenseMatrix<Real> &U = trial.values[qp];
                const Real w = scale * weights[qp];

                for (unsigned int j = 0; j < U.n(); ++j)
                {
                    Real cu[kMaxRows] = {};
                    for (unsigned int a = 0; a < rows_test; ++a)
                    {
                        if (coefficient.is_tensor)
                            for (unsigned int b = 0; b < rows_trial; ++b)
                                cu[a] += coefficient.tensor(test.components[a], trial.components[b]) * U(b, j);
                        else
                            cu[a] = coefficient.scalar * U(a, j);
                    }
                    for (unsigned int i = 0; i < T.n(); ++i)
                    {
                        Real sum = 0.0;
                        for (unsigned int a = 0; a < rows_test; ++a)
                            sum += T(a, i) * cu[a];
                        block(i, j) += w * sum;
                    }
                }
            }
            return block;
        }

        // integral of test^T * vector: the source vector of one kernel over the test DOFs.
        DenseVector<Real> weak_form_source(const Nabla::Operand &test,
                                           const Vec3 &vector,
                                           const std::vector<Real> &weights,
                                           Real scale)
        {
            if (test.values.size() != weights.size())
                fail("test and weights must have one entry per quadrature point.");
            for (unsigned int c : test.components)
                if (c > 2)
                    fail("a source kernel needs physical components (x, y, z) on the test side.");

            DenseVector<Real> source(test.values.empty() ? 0 : test.values[0].n());
            for (std::size_t qp = 0; qp < test.values.size(); ++qp)
            {
                const DenseMatrix<Real> &T = test.values[qp];
                for (unsigned int i = 0; i < T.n(); ++i)
                {
                    Real sum = 0.0;
                    for (unsigned int a = 0; a < T.m(); ++a)
                        sum += T(a, i) * vector[test.components[a]];
                    source(i) += scale * weights[qp] * sum;
                }
            }
            return source;
        }

        // Adds @p block to K at the matrix indices of the test and trial variables' columns.
        void add_block(DenseMatrix<Real> &K,
                       const std::vector<unsigned int> &test_columns,
                       const std::vector<unsigned int> &trial_columns,
                       const DenseMatrix<Real> &block)
        {
            if (block.m() != test_columns.size() || block.n() != trial_columns.size())
                fail("a kernel block does not match the DOF columns of its variables.");
            for (unsigned int i = 0; i < block.m(); ++i)
                for (unsigned int j = 0; j < block.n(); ++j)
                    K(test_columns[i], trial_columns[j]) += block(i, j);
        }

        // Operator evaluations are shared between kernels and terms of one element.
        class OperandCache
        {
        public:
            OperandCache(const ElementTypeInstance &element_type, const libMesh::Elem &elem, Real thickness)
                : element_type_(element_type), elem_(elem), thickness_(thickness) {}

            const Nabla::Operand &get(NablaOperator op, const std::string &variable)
            {
                const auto key = std::make_pair(op, variable);
                auto it = cache_.find(key);
                if (it == cache_.end())
                    it = cache_.emplace(key, Nabla::operand(element_type_, variable, op, elem_, thickness_)).first;
                return it->second;
            }

        private:
            const ElementTypeInstance &element_type_;
            const libMesh::Elem &elem_;
            Real thickness_;
            std::map<std::pair<NablaOperator, std::string>, Nabla::Operand> cache_;
        };

        // Everything one element evaluation needs.
        struct Context
        {
            Context(const ElementTypeInstance &element_type, const libMesh::Elem &elem, const ElementParameters &parameters)
                : thickness(element_type.dimension() == 2 ? parameters.constants.thickness : 1.0),
                  cache(element_type, elem, thickness),
                  dofs(parameters.dofs),
                  weights(Nabla::integration_weights(element_type, elem, thickness))
            {
            }

            Real thickness;
            OperandCache cache;
            const ElementDofs &dofs;
            std::vector<Real> weights;
        };

        const std::vector<unsigned int> &columns_of(const ElementDofs &dofs, const std::string &variable)
        {
            const auto it = dofs.columns.find(variable);
            if (it == dofs.columns.end())
                fail("variable '" + variable + "' has no DOF columns on this element.");
            return it->second;
        }

        // Physical 3-vector of an operand applied to a local solution, at each quadrature point.
        std::vector<Vec3> sample_operand(const Nabla::Operand &operand,
                                         const std::vector<unsigned int> &columns,
                                         const DenseVector<Real> &solution)
        {
            for (unsigned int c : operand.components)
                if (c > 2)
                    fail("operator result has no physical components to evaluate fields from.");
            for (unsigned int column : columns)
                if (column >= solution.size())
                    fail("solution vector is smaller than the element DOF count.");

            std::vector<Vec3> samples(operand.values.size(), Vec3{0.0, 0.0, 0.0});
            for (std::size_t qp = 0; qp < samples.size(); ++qp)
            {
                const DenseMatrix<Real> &m = operand.values[qp];
                if (m.n() != columns.size())
                    fail("an operator does not match the DOF columns of its variable.");
                for (unsigned int r = 0; r < m.m(); ++r)
                    for (unsigned int j = 0; j < m.n(); ++j)
                        samples[qp][operand.components[r]] += m(r, j) * solution(columns[j]);
            }
            return samples;
        }
    } // namespace

    ElementKernels::ElementKernels(const ElementTypeInstance &element_type) : element_type_(element_type)
    {
        const BaseElementType &base = element_type.base();

        const auto is_dof = [&element_type](const std::string &name)
        {
            const VariableSpec *variable = element_type.find_variable(name);
            return variable && variable->kind == FieldKind::Dof;
        };

        for (const KernelSpec &spec : base.kernels)
        {
            if (!spec.domains.empty() &&
                std::find(spec.domains.begin(), spec.domains.end(), element_type.domain_name()) == spec.domains.end())
                continue;

            Kernel kernel;
            kernel.name = spec.name;
            kernel.has_trial = !spec.trial.empty();
            if (kernel.has_trial)
                kernel.trial = OperatorExpression::parse(spec.trial);
            kernel.test = OperatorExpression::parse(spec.test);
            if (!is_dof(kernel.test.variable) || (kernel.has_trial && !is_dof(kernel.trial.variable)))
                continue; // A variable of the kernel does not exist in this domain.

            kernel.sign = spec.coeff[0] == '-' ? -1.0 : 1.0;
            kernel.coefficient = spec.coeff[0] == '-' ? spec.coeff.substr(1) : spec.coeff;
            if (!kernel.has_trial)
            {
                if (!is_vector_coefficient(kernel.coefficient))
                    fail("kernel '" + kernel.name + "' has no trial, so its coefficient must be a vector (X[]).");
                const std::string load = kernel.coefficient.substr(0, kernel.coefficient.size() - 2);
                if (kernel.coefficient != kMagnetization &&
                    std::find(element_load_names_.begin(), element_load_names_.end(), load) == element_load_names_.end())
                    element_load_names_.push_back(load);
            }
            else if (is_vector_coefficient(kernel.coefficient))
                fail("kernel '" + kernel.name + "' has a trial, so its coefficient cannot be a vector.");
            kernels_.push_back(std::move(kernel));
        }

        for (const VariableSpec &variable : element_type.variables())
        {
            if (variable.kind == FieldKind::Dof)
                continue;

            const auto spec = std::find_if(base.fields.begin(), base.fields.end(), [&variable](const FieldSpec &f)
                                           { return f.name == variable.name && f.kind == variable.kind; });
            if (spec == base.fields.end())
                fail("field '" + variable.name + "' not found in the base element type.");

            Field field;
            field.variable = &variable;
            for (ExpressionTerm &term : ExpressionTerm::parse_sum(spec->expression, base.coefficients))
                if (is_dof(term.operand.variable))
                    field.terms.push_back(std::move(term));
            (variable.kind == FieldKind::Nodal ? nodal_ : element_).push_back(std::move(field));
        }
    }

    void ElementKernels::build_matrices(const libMesh::Elem &elem,
                                        const ElementParameters &parameters,
                                        DenseMatrix<Real> &Ke,
                                        DenseMatrix<Real> &Ce) const
    {
        Context ctx(element_type_, elem, parameters);
        Ke.resize(ctx.dofs.n_dofs, ctx.dofs.n_dofs);
        Ce.resize(ctx.dofs.n_dofs, ctx.dofs.n_dofs);
        const bool nonlinear = parameters.bh_curve != nullptr;

        for (const Kernel &kernel : kernels_)
        {
            if (!kernel.has_trial || (nonlinear && kernel.coefficient == kBHCoefficient))
                continue;

            const Coefficient coefficient = evaluate_coefficient(kernel.coefficient, parameters, element_type_);
            const Nabla::Operand &test = ctx.cache.get(kernel.test.op, kernel.test.variable);
            const Nabla::Operand &trial = ctx.cache.get(kernel.trial.op, kernel.trial.variable);
            DenseMatrix<Real> &target = kernel.trial.time_derivative ? Ce : Ke;

            const DenseMatrix<Real> block = weak_form_block(test, coefficient, trial, ctx.weights, kernel.sign);
            add_block(target, columns_of(ctx.dofs, kernel.test.variable), columns_of(ctx.dofs, kernel.trial.variable), block);
        }
    }

    void ElementKernels::build_rhs(const libMesh::Elem &elem,
                                   const ElementParameters &parameters,
                                   const ElementLoads &loads,
                                   DenseVector<Real> &rhs) const
    {
        Context ctx(element_type_, elem, parameters);
        rhs.resize(ctx.dofs.n_dofs);

        for (const Kernel &kernel : kernels_)
        {
            if (kernel.has_trial)
                continue;

            const Coefficient coefficient = evaluate_coefficient(kernel.coefficient, parameters, element_type_, &loads);
            const Nabla::Operand &test = ctx.cache.get(kernel.test.op, kernel.test.variable);
            const DenseVector<Real> source = weak_form_source(test, coefficient.vector, ctx.weights, kernel.sign);
            const std::vector<unsigned int> &columns = columns_of(ctx.dofs, kernel.test.variable);
            if (source.size() != columns.size())
                fail("the source of kernel '" + kernel.name + "' does not match the DOF columns of its variable.");
            for (unsigned int i = 0; i < source.size(); ++i)
                rhs(columns[i]) += source(i);
        }

        // Nodal loads are point values on the nodes, not integrals.
        for (const VariableSpec &variable : element_type_.variables())
        {
            if (variable.kind != FieldKind::Dof || variable.nodal_loads.empty())
                continue;
            const auto load = loads.nodal.find(variable.name);
            if (load == loads.nodal.end())
                continue;

            const std::vector<unsigned int> &columns = columns_of(ctx.dofs, variable.name);
            const std::size_t n = std::min(columns.size(), load->second.values.size());
            for (std::size_t i = 0; i < n; ++i)
                if (load->second.inject[i])
                    rhs(columns[i]) += load->second.values[i];
        }
    }

    void ElementKernels::build_tangent_residual(const libMesh::Elem &elem,
                                                const ElementParameters &parameters,
                                                const DenseVector<Real> &u,
                                                DenseMatrix<Real> &tangent,
                                                DenseVector<Real> &residual) const
    {
        Context ctx(element_type_, elem, parameters);
        tangent.resize(ctx.dofs.n_dofs, ctx.dofs.n_dofs);
        residual.resize(ctx.dofs.n_dofs);
        if (!parameters.bh_curve)
            return;
        if (u.size() < ctx.dofs.n_dofs)
            fail("solution vector is smaller than the element DOF count.");

        const unsigned int dim = element_type_.dimension();
        for (const Kernel &kernel : kernels_)
        {
            if (!kernel.has_trial || kernel.coefficient != kBHCoefficient)
                continue;

            const Nabla::Operand &test = ctx.cache.get(kernel.test.op, kernel.test.variable);
            const Nabla::Operand &trial = ctx.cache.get(kernel.trial.op, kernel.trial.variable);
            const std::vector<unsigned int> &test_columns = columns_of(ctx.dofs, kernel.test.variable);
            const std::vector<unsigned int> &trial_columns = columns_of(ctx.dofs, kernel.trial.variable);
            for (unsigned int c : test.components)
                if (c > 2)
                    fail("kernel '" + kernel.name + "': the BH coefficient needs physical components.");
            for (unsigned int c : trial.components)
                if (c > 2)
                    fail("kernel '" + kernel.name + "': the BH coefficient needs physical components.");

            for (std::size_t qp = 0; qp < ctx.weights.size(); ++qp)
            {
                const DenseMatrix<Real> &T = test.values[qp];
                const DenseMatrix<Real> &U = trial.values[qp];
                const Real w = kernel.sign * ctx.weights[qp];

                Real B[3] = {0.0, 0.0, 0.0};
                for (unsigned int r = 0; r < U.m(); ++r)
                    for (unsigned int j = 0; j < U.n(); ++j)
                        B[trial.components[r]] += U(r, j) * u(trial_columns[j]);
                Real H[3];
                Real jac[3][3];
                evaluate_bh(*parameters.bh_curve, B, dim, H, jac);

                for (unsigned int i = 0; i < T.n(); ++i)
                {
                    Real h_term = 0.0;
                    for (unsigned int a = 0; a < T.m(); ++a)
                        h_term += T(a, i) * H[test.components[a]];
                    residual(test_columns[i]) += w * h_term;

                    for (unsigned int j = 0; j < U.n(); ++j)
                    {
                        Real tangent_term = 0.0;
                        for (unsigned int a = 0; a < T.m(); ++a)
                            for (unsigned int b = 0; b < U.m(); ++b)
                                tangent_term += T(a, i) * jac[test.components[a]][trial.components[b]] * U(b, j);
                        tangent(test_columns[i], trial_columns[j]) += w * tangent_term;
                    }
                }
            }
        }
    }

    std::vector<DenseVector<Real>> ElementKernels::element_fields(const libMesh::Elem &elem,
                                                                  const ElementParameters &parameters,
                                                                  const DenseVector<Real> &u,
                                                                  const DenseVector<Real> &u_dt) const
    {
        Context ctx(element_type_, elem, parameters);

        std::vector<DenseVector<Real>> out;
        for (const Field &field : element_)
        {
            Vec3 total{0.0, 0.0, 0.0};
            for (const ExpressionTerm &term : field.terms)
            {
                // Each operand is averaged over the quadrature points, then scaled and combined.
                const Nabla::Operand &operand = ctx.cache.get(term.operand.op, term.operand.variable);
                const auto samples = sample_operand(operand, columns_of(ctx.dofs, term.operand.variable),
                                                    term.operand.time_derivative ? u_dt : u);

                Vec3 mean{0.0, 0.0, 0.0};
                for (const Vec3 &sample : samples)
                    for (unsigned int a = 0; a < 3; ++a)
                        mean[a] += sample[a] / samples.size();

                const Vec3 value = apply_coefficient(term.coefficient, mean, parameters, element_type_);
                for (unsigned int a = 0; a < 3; ++a)
                    total[a] += term.sign * value[a];
            }

            const auto &components = field.variable->active_components;
            DenseVector<Real> values(static_cast<unsigned int>(components.size()));
            for (unsigned int k = 0; k < components.size(); ++k)
                values(k) = total[components[k]];
            out.push_back(std::move(values));
        }
        return out;
    }

    std::vector<DenseVector<Real>> ElementKernels::nodal_fields(const libMesh::Elem &elem,
                                                                const ElementParameters &parameters,
                                                                const DenseVector<Real> &u,
                                                                const DenseVector<Real> &u_dt) const
    {
        Context ctx(element_type_, elem, parameters);

        std::vector<DenseVector<Real>> out;
        for (const Field &field : nodal_)
        {
            const Nabla::Operand &basis = ctx.cache.get(NablaOperator::Identity, field.variable->name);
            const unsigned int n = basis.values[0].n();
            const unsigned int rows = basis.values[0].m();

            // Each basis row is a component; its k-th column belongs to node k.
            std::vector<std::vector<unsigned int>> node_columns(rows);
            for (unsigned int i = 0; i < n; ++i)
            {
                unsigned int best_row = 0;
                Real best = -1.0;
                for (unsigned int r = 0; r < rows; ++r)
                {
                    Real magnitude = 0.0;
                    for (const auto &m : basis.values)
                        magnitude += std::abs(m(r, i));
                    if (magnitude > best)
                    {
                        best = magnitude;
                        best_row = r;
                    }
                }
                node_columns[best_row].push_back(i);
            }
            const std::size_t n_nodes = node_columns[0].size();
            for (const auto &columns : node_columns)
                if (columns.size() != n_nodes)
                    fail("field '" + field.variable->name + "': components have different node counts.");

            // L2 projection: M x = integral of basis^T * value.
            Coefficient identity;
            identity.scalar = 1.0;
            const DenseMatrix<Real> mass = weak_form_block(basis, identity, basis, ctx.weights, 1.0);

            DenseVector<Real> total(n);
            for (const ExpressionTerm &term : field.terms)
            {
                // Each operand is projected onto the nodal basis, then scaled node by node and combined.
                const Nabla::Operand &operand = ctx.cache.get(term.operand.op, term.operand.variable);
                const auto samples = sample_operand(operand, columns_of(ctx.dofs, term.operand.variable),
                                                    term.operand.time_derivative ? u_dt : u);

                DenseVector<Real> rhs(n);
                for (std::size_t qp = 0; qp < samples.size(); ++qp)
                    for (unsigned int r = 0; r < rows; ++r)
                        for (unsigned int i = 0; i < n; ++i)
                            rhs(i) += ctx.weights[qp] * basis.values[qp](r, i) * samples[qp][basis.components[r]];
                DenseVector<Real> projected;
                mass.cholesky_solve(rhs, projected);

                for (std::size_t k = 0; k < n_nodes; ++k)
                {
                    Vec3 nodal{0.0, 0.0, 0.0};
                    for (unsigned int r = 0; r < rows; ++r)
                        nodal[basis.components[r]] = projected(node_columns[r][k]);
                    const Vec3 value = apply_coefficient(term.coefficient, nodal, parameters, element_type_);
                    for (unsigned int r = 0; r < rows; ++r)
                        total(node_columns[r][k]) += term.sign * value[basis.components[r]];
                }
            }
            out.push_back(std::move(total));
        }
        return out;
    }
} // namespace pytec3d

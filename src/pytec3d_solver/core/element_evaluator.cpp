/**
 * @file element_evaluator.cpp
 * @brief ElementEvaluator: one implementation per element type, compiled from its spec.
 *
 * Evaluator<Id> walks the constexpr tables of formulation::spec(Id) with compile-time loops (static_for): the
 * index of every term, kernel and evaluation is a template argument, so the value types of each kernel (test,
 * coefficient, trial) are known and the kernel templates of kernels.h are called
 * directly. What remains at run time is per element: the domain of the element, the laws selected there, the values.
 */
#include "element_evaluator.h"

#include <array>
#include <cmath>
#include <string>
#include <type_traits>
#include <utility>

#include "formulations.h"
#include "kernels.h"
#include "logging.h"
#include "model_data.h"

namespace pytec3d
{
    using formulation::CoefficientKind;
    using formulation::EvaluationKind;
    using formulation::Stage;
    using formulation::ValueType;
    using libMesh::Real;
    using libMesh::RealTensorValue;
    using libMesh::RealVectorValue;

    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementEvaluator: " + message);
        }

        template <std::size_t... I, typename F>
        void for_each_index(std::index_sequence<I...>, F &&f)
        {
            (f(std::integral_constant<std::size_t, I>{}), ...);
        }

        /** @brief f(integral_constant<0>) ... f(integral_constant<N - 1>). */
        template <std::size_t N, typename F>
        void static_for(F &&f)
        {
            for_each_index(std::make_index_sequence<N>{}, std::forward<F>(f));
        }

        /** @brief C++ type of a value type. */
        template <ValueType Type>
        struct cpp_type
        {
            using type = Real;
        };
        template <>
        struct cpp_type<ValueType::Vector>
        {
            using type = RealVectorValue;
        };
        template <>
        struct cpp_type<ValueType::Tensor>
        {
            using type = RealTensorValue;
        };
        template <ValueType Type>
        using cpp_t = typename cpp_type<Type>::type;

        /** @brief The values of a kernel coefficient: one, or one per quadrature point. */
        template <typename C>
        struct TypedCoefficient
        {
            std::vector<C> values;
            bool per_qp = false;
            kernel::Coefficient<C> view() const { return {values.data(), per_qp}; }
        };

        /** @brief The evaluation a term belongs to. */
        constexpr std::size_t evaluation_of(const formulation::ElementTypeSpec &spec, std::size_t term)
        {
            for (std::size_t e = 0; e < spec.evaluations.size(); ++e)
                if (term >= spec.evaluations[e].first_term && term < spec.evaluations[e].first_term + spec.evaluations[e].n_terms)
                    return e;
            return spec.evaluations.size();
        }

        constexpr std::size_t law_position(const formulation::ElementTypeSpec &spec, registry::LawId id)
        {
            for (std::size_t l = 0; l < spec.laws.size(); ++l)
                if (spec.laws[l].id == id)
                    return l;
            return spec.laws.size();
        }

        constexpr std::size_t output_position(const formulation::LawSpec &law, registry::LawOutputId id)
        {
            for (std::size_t o = 0; o < law.outputs.size(); ++o)
                if (law.outputs[o] == id)
                    return o;
            return law.outputs.size();
        }

        /** @brief A kernel on the element, in the shape its spec gives (kernels.h). */
        struct KernelResult
        {
            bool done = false;
            libMesh::DenseMatrix<Real> matrix;  ///< Test and trial.
            libMesh::DenseVector<Real> vector;  ///< Test only.
            ValueTable table;                   ///< Trial only.
            std::vector<ParameterValue> values; ///< Neither: per quadrature point.
        };

        struct LawState
        {
            const CurveData *curve = nullptr; ///< Of the model data row selecting the law, if the law has a curve.
            bool done = false;
            std::vector<std::vector<ParameterValue>> outputs; ///< [output][qp].
        };

        template <registry::ElementTypeId Id>
        class Evaluator final : public ElementEvaluator::Implementation
        {
            using Matrix = ElementEvaluator::Matrix;
            using Vector = ElementEvaluator::Vector;
            using Inputs = ElementEvaluator::Inputs;

            static constexpr const formulation::ElementTypeSpec &spec = formulation::kElementTypes[static_cast<std::size_t>(Id)];
            static constexpr std::size_t kTerms = spec.terms.size();
            static constexpr std::size_t kKernels = spec.kernels.size();
            static constexpr std::size_t kEvaluations = spec.evaluations.size();
            static constexpr std::size_t kLaws = spec.laws.size();
            static constexpr std::size_t kNablas = spec.nablas.size();

        public:
            explicit Evaluator(const Inputs &inputs) : in_(inputs), type_(*inputs.type)
            {
                const libMesh::FEAbstract *geometry = type_.reinit(*in_.elem, in_.domain, in_.operation, in_.laws);
                if (geometry)
                    weights_ = geometry->get_JxW();
                nq_ = weights_.size();
                length_.assign(nq_, 1.0);
                if (geometry && spec.dimension == 2)
                {
                    Real thickness = 1.0;
                    if constexpr (spec.out_of_plane_length >= 0)
                        thickness = in_.parameters->symbol(static_cast<registry::ParameterId>(spec.out_of_plane_length)).scalar();
                    const auto &xyz = geometry->get_xyz();
                    for (std::size_t qp = 0; qp < nq_; ++qp)
                    {
                        // Out-of-plane length: the declared length, times the circumference of the axisymmetric element.
                        const Real r = std::sqrt(xyz[qp](0) * xyz[qp](0) + xyz[qp](2) * xyz[qp](2));
                        length_[qp] = type_.options().axisymmetric ? thickness * 2.0 * libMesh::pi * r : thickness;
                        weights_[qp] *= length_[qp];
                    }
                }

                // The curves of the laws selected, keyed by the row that selects each law.
                for (std::size_t l = 0; l < kLaws; ++l)
                {
                    const registry::LawId id = spec.laws[l].id;
                    const registry::LawInfo &info = registry::info(id);
                    if (!(in_.laws & registry::bit(id)) || info.curve == registry::CurveId::COUNT)
                        continue;
                    const unsigned int row = info.source == registry::LawSource::Material   ? in_.material_id
                                             : info.source == registry::LawSource::Constant ? in_.constants_id
                                                                                            : in_.type_id;
                    laws_[l].curve = in_.model_data ? in_.model_data->curve(info.curve, row) : nullptr;
                    if (!laws_[l].curve)
                    {
                        const auto &source = registry::kCurveSources[static_cast<std::size_t>(info.curve)];
                        fail("law '" + std::string(registry::name(id)) + "' is selected by row " + std::to_string(row) +
                             " but the model data has no '" + std::string(source.section) + "' curve for it.");
                    }
                }
            }

            // ---- evaluations -------------------------------------------------------------------------------

            void add_matrices(Matrix &K, Matrix &C, Matrix &M) override
            {
                for_terms(Operation::Matrices, [&](auto term_index)
                          {
                    if constexpr (spec.evaluations[evaluation_of(spec, term_index)].kind == EvaluationKind::Block)
                    {
                        constexpr formulation::TermSpec term = spec.terms[term_index];
                        Matrix &target = term.trial_dt == 0 ? K : term.trial_dt == 1 ? C : M;
                        scatter(kernel<term.kernel>().matrix, term, target);
                    } });
            }

            void field_operators(FieldKind kind, const ElementType::FieldIndices &indices, std::array<Matrix, 3> &F,
                                 Vector &f) override
            {
                if (in_.operation != Operation::Matrices)
                    fail("the evaluator was built for another operation.");
                // Only the operators the spec can fill; the others stay empty (0 x 0).
                const unsigned int operators = type_.field_operators(kind, in_.domain);
                for (std::size_t k = 0; k < F.size(); ++k)
                    F[k].resize(operators & (1u << k) ? indices.n_dofs : 0, operators & (1u << k) ? in_.dofs->n_dofs : 0);
                f.resize(operators & 8u ? indices.n_dofs : 0);
                if (nq_ == 0)
                    return;

                static_for<kEvaluations>([&](auto e)
                                         {
                    constexpr formulation::EvaluationSpec evaluation = spec.evaluations[e];
                    if constexpr (evaluation.kind == EvaluationKind::Average || evaluation.kind == EvaluationKind::L2Projection)
                    {
                        constexpr formulation::FieldSpec field = spec.fields[static_cast<std::size_t>(evaluation.field)];
                        if (field.kind != kind || !ElementType::active(field, in_.domain) ||
                            type_.dynamic(evaluation, in_.domain, in_.laws))
                            return;
                        const auto &ids = indices.field[static_cast<std::size_t>(evaluation.field)];
                        if constexpr (evaluation.kind == EvaluationKind::Average)
                            average_operator<e>(ids, F, f);
                        else
                            projection_operator<e>(ids, F, f);
                    } });
            }

            void add_rhs(Vector &F) override
            {
                for_terms(Operation::Rhs, [&](auto term_index)
                          {
                    constexpr formulation::TermSpec term = spec.terms[term_index];
                    scatter(kernel<term.kernel>().vector, term, F, 1.0); });
            }

            void add_nonlinear(Matrix &J, Vector &R) override
            {
                for_terms(Operation::Nonlinear, [&](auto term_index)
                          {
                    constexpr formulation::TermSpec term = spec.terms[term_index];
                    if constexpr (spec.evaluations[evaluation_of(spec, term_index)].kind == EvaluationKind::Rhs)
                        scatter(kernel<term.kernel>().vector, term, R, 1.0);
                    else if constexpr (term.stage == Stage::Tangent)
                    {
                        // d(dt^d u)/du = c_d of the time scheme; nothing on a static step.
                        const Real scale = in_.jacobian_scale[term.trial_dt];
                        if (scale != 0.0)
                            scatter(kernel<term.kernel>().matrix, term, J, scale);
                    }
                    else
                    {
                        // The residual of the term at the current solution (its rates for a dt() trial).
                        if (term.trial_dt > 0 && in_.jacobian_scale[term.trial_dt] == 0.0)
                            return;
                        const Matrix &m = kernel<term.kernel>().matrix;
                        const Vector u = trial_solution(term);
                        const auto &rows = in_.dofs->field[static_cast<std::size_t>(term.test_field)];
                        for (unsigned int i = 0; i < m.m(); ++i)
                        {
                            Real sum = 0.0;
                            for (unsigned int j = 0; j < m.n(); ++j)
                                sum += m(i, j) * u(j);
                            R(rows[i]) += sum;
                        }
                    } });
            }

            void field_values(FieldKind kind, const ElementType::FieldIndices &indices, Vector &values) override
            {
                const Operation operation = kind == FieldKind::Nodal ? Operation::NodalFields : Operation::ElementFields;
                if (operation != in_.operation)
                    fail("the evaluator was built for another operation.");
                if (nq_ == 0)
                    return; // Nothing to evaluate in this domain.

                static_for<kEvaluations>([&](auto e)
                                         {
                    constexpr formulation::EvaluationSpec evaluation = spec.evaluations[e];
                    if constexpr (evaluation.kind == EvaluationKind::Average || evaluation.kind == EvaluationKind::L2Projection)
                    {
                        constexpr formulation::FieldSpec field = spec.fields[static_cast<std::size_t>(evaluation.field)];
                        if (field.kind != kind || !ElementType::active(field, in_.domain) ||
                            !type_.dynamic(evaluation, in_.domain, in_.laws))
                            return;
                        const auto &ids = indices.field[static_cast<std::size_t>(evaluation.field)];
                        if constexpr (evaluation.kind == EvaluationKind::Average)
                            average<e>(ids, values);
                        else
                            project<e>(ids, values);
                    } });
            }

            void check_parameters() override
            {
                static_for<kTerms>([&](auto term_index)
                                   {
                    constexpr formulation::TermSpec term = spec.terms[term_index];
                    constexpr formulation::CoefficientSpec c = spec.kernels[term.kernel].coefficient;
                    const auto &evaluation = spec.evaluations[evaluation_of(spec, term_index)];
                    if (!type_.uses(evaluation, term, Operation::Parameters, in_.laws, in_.domain))
                        return;
                    if constexpr (c.kind == CoefficientKind::Parameter)
                        in_.parameters->symbol(static_cast<registry::ParameterId>(c.index)); });
            }

        private:
            // ---- terms ----------------------------------------------------------------------------------------

            /** @brief Calls fn(term index) for each term @p operation computes on the element. */
            template <typename Fn>
            void for_terms(Operation operation, Fn &&fn)
            {
                // The FEs reinit()ed are those of the operation of the inputs.
                if (operation != in_.operation)
                    fail("the evaluator was built for another operation.");
                static_for<kTerms>([&](auto term_index)
                                   {
                    const auto &evaluation = spec.evaluations[evaluation_of(spec, term_index)];
                    if (type_.uses(evaluation, spec.terms[term_index], operation, in_.laws, in_.domain))
                        fn(term_index); });
            }

            void scatter(const Matrix &m, const formulation::TermSpec &term, Matrix &target, Real scale = 1.0) const
            {
                const auto &rows = in_.dofs->field[static_cast<std::size_t>(term.test_field)];
                const auto &columns = in_.dofs->field[static_cast<std::size_t>(term.trial_field)];
                for (unsigned int i = 0; i < m.m(); ++i)
                    for (unsigned int j = 0; j < m.n(); ++j)
                        target(rows[i], columns[j]) += scale * m(i, j);
            }

            void scatter(const Vector &v, const formulation::TermSpec &term, Vector &target, Real sign) const
            {
                const auto &rows = in_.dofs->field[static_cast<std::size_t>(term.test_field)];
                for (unsigned int i = 0; i < v.size(); ++i)
                    target(rows[i]) += sign * v(i);
            }

            Vector trial_solution(const formulation::TermSpec &term) const
            {
                const Vector *solution = term.trial_dt == 0 ? in_.solution : term.trial_dt == 1 ? in_.solution_dt : in_.solution_ddt;
                if (!solution || !in_.dofs)
                    fail("a term needs the solution, which was not given.");
                const auto &ids = in_.dofs->field[static_cast<std::size_t>(term.trial_field)];
                Vector u(static_cast<unsigned int>(ids.size()));
                for (std::size_t j = 0; j < ids.size(); ++j)
                    u(static_cast<unsigned int>(j)) = (*solution)(ids[j]);
                return u;
            }

            /** @brief The term applies on the element: in its domain, and selected by the laws of the element. */
            bool used(const formulation::TermSpec &term) const
            {
                return ElementType::active(term, in_.domain) && ElementType::selected(term, in_.laws);
            }

            // ---- nablas, coefficients, laws -------------------------------------------------------------------

            const ValueTable &nabla(std::size_t index)
            {
                if (!table_done_[index])
                {
                    type_.nabla(index).evaluate(length_, tables_[index]);
                    table_done_[index] = true;
                }
                return tables_[index];
            }

            template <typename C>
            TypedCoefficient<C> coefficient(const formulation::CoefficientSpec &c)
            {
                TypedCoefficient<C> out;
                switch (c.kind)
                {
                case CoefficientKind::One:
                    if constexpr (std::is_same_v<C, Real>)
                        out.values.push_back(c.sign);
                    break;
                case CoefficientKind::Parameter:
                    out.values.push_back(c.sign * in_.parameters->symbol(static_cast<registry::ParameterId>(c.index)).template get<C>());
                    break;
                case CoefficientKind::LawOutput:
                    // Only in the law variants, used where the law is selected.
                    out.per_qp = true;
                    for (const ParameterValue &value : law_output(static_cast<registry::LawOutputId>(c.index)))
                        out.values.push_back(c.sign * value.template get<C>());
                    break;
                }
                return out;
            }

            const std::vector<ParameterValue> &law_output(registry::LawOutputId id)
            {
                const std::size_t l = law_position(spec, registry::info(id).law);
                const formulation::LawSpec &law = spec.laws[l];
                LawState &state = laws_[l];
                if (!state.done)
                {
                    const std::vector<ParameterValue> input = point_values(static_cast<std::size_t>(law.input));
                    state.outputs.assign(law.outputs.size(), std::vector<ParameterValue>(nq_));
                    std::vector<ParameterValue> point(law.outputs.size());
                    for (std::size_t qp = 0; qp < nq_; ++qp)
                    {
                        formulation::evaluate_law(Id, law.id, input[qp], state.curve, *in_.parameters, point.data());
                        for (std::size_t o = 0; o < point.size(); ++o)
                            state.outputs[o][qp] = point[o];
                    }
                    state.done = true;
                }
                return state.outputs[output_position(law, id)];
            }

            // ---- kernels ------------------------------------------------------------------------------------

            template <std::size_t K>
            const KernelResult &kernel()
            {
                KernelResult &result = kernels_[K];
                if (result.done)
                    return result;

                constexpr formulation::KernelSpec ks = spec.kernels[K];
                using C = cpp_t<ks.coefficient.type>;
                const TypedCoefficient<C> c = coefficient<C>(ks.coefficient);

                if constexpr (ks.test >= 0 && ks.trial >= 0)
                {
                    constexpr formulation::NablaSpec tn = spec.nablas[static_cast<std::size_t>(ks.test)];
                    constexpr formulation::NablaSpec un = spec.nablas[static_cast<std::size_t>(ks.trial)];
                    using Test = cpp_t<tn.output>;
                    using Trial = cpp_t<un.output>;
                    static_assert(kernel::has_block<Test, C, Trial>, "no product for this kernel (checked by the generator)");
                    // The mass matrices of the L2 projections are symmetric: only their upper triangle is computed.
                    constexpr bool symmetric = [] {
                        for (const auto &evaluation : spec.evaluations)
                            if (evaluation.mass_kernel == static_cast<signed short>(K))
                                return true;
                        return false;
                    }();
                    kernel::matrix<symmetric>(nabla(ks.test).template get<Test>(), c.view(),
                                              nabla(ks.trial).template get<Trial>(), weights_, result.matrix);
                }
                else if constexpr (ks.test >= 0)
                {
                    constexpr formulation::NablaSpec tn = spec.nablas[static_cast<std::size_t>(ks.test)];
                    using Test = cpp_t<tn.output>;
                    static_assert(kernel::has_source<Test, C>, "no product for this kernel (checked by the generator)");
                    kernel::vector(nabla(ks.test).template get<Test>(), c.view(), weights_, result.vector);
                }
                else if constexpr (ks.trial >= 0)
                {
                    constexpr formulation::NablaSpec un = spec.nablas[static_cast<std::size_t>(ks.trial)];
                    using Trial = cpp_t<un.output>;
                    using Out = cpp_t<ks.output>;
                    static_assert(kernel::has_value<C, Trial>, "no product for this kernel (checked by the generator)");
                    result.table.type = ks.output;
                    kernel::table(c.view(), nabla(ks.trial).template get<Trial>(), result.table.template get<Out>());
                }
                else
                {
                    result.values.clear();
                    for (std::size_t qp = 0; qp < nq_; ++qp)
                        result.values.emplace_back(c.values[c.per_qp ? qp : 0]);
                }
                result.done = true;
                return result;
            }

            // ---- field evaluations ----------------------------------------------------------------------------

            /** @brief Sum over the terms of evaluation E of C * trial(u) at each quadrature point. */
            template <std::size_t E>
            std::vector<ParameterValue> point_values()
            {
                constexpr formulation::EvaluationSpec evaluation = spec.evaluations[E];
                constexpr ValueType value = spec.fields[static_cast<std::size_t>(evaluation.field)].value;
                using Out = cpp_t<value>;
                std::vector<ParameterValue> sum(nq_, ParameterValue::zero(value));
                static_for<evaluation.n_terms>([&](auto i)
                                               {
                    constexpr formulation::TermSpec term = spec.terms[evaluation.first_term + i];
                    if (!used(term))
                        return;
                    const KernelResult &k = kernel<term.kernel>();
                    if constexpr (term.trial_field < 0)
                    {
                        for (std::size_t qp = 0; qp < nq_; ++qp)
                            sum[qp] += k.values[qp];
                    }
                    else
                    {
                        const Vector u = trial_solution(term);
                        const auto &table = k.table.template get<Out>();
                        for (std::size_t qp = 0; qp < nq_; ++qp)
                        {
                            Out &s = sum[qp].template get<Out>();
                            for (unsigned int j = 0; j < u.size(); ++j)
                                s += table[qp][j] * u(j);
                        }
                    } });
                return sum;
            }

            std::vector<ParameterValue> point_values(std::size_t evaluation)
            {
                std::vector<ParameterValue> out;
                static_for<kEvaluations>([&](auto e)
                                         {
                    if constexpr (spec.evaluations[e].kind != EvaluationKind::Block && spec.evaluations[e].kind != EvaluationKind::Rhs)
                        if (e == evaluation)
                            out = point_values<e>(); });
                return out;
            }

            template <std::size_t E>
            void average(const ElementType::FieldRange &ids, Vector &values)
            {
                constexpr formulation::FieldSpec field = spec.fields[static_cast<std::size_t>(spec.evaluations[E].field)];
                const std::vector<Real> &w = mean_weights();
                const std::vector<ParameterValue> points = point_values<E>();
                ParameterValue mean = ParameterValue::zero(field.value);
                for (std::size_t qp = 0; qp < nq_; ++qp)
                    mean += points[qp].scaled(w[qp]);
                for (std::size_t c = 0; c < ids.size(); ++c)
                {
                    if constexpr (field.value == ValueType::Vector)
                        values(ids[c]) = mean.vector()(field.components[c]);
                    else
                        values(ids[c]) = mean.scalar();
                }
            }

            template <std::size_t E>
            void project(const ElementType::FieldRange &ids, Vector &values)
            {
                constexpr formulation::EvaluationSpec evaluation = spec.evaluations[E];
                constexpr auto mass_kernel = static_cast<std::size_t>(evaluation.mass_kernel);
                Vector rhs(static_cast<unsigned int>(ids.size()));
                static_for<evaluation.n_terms>([&](auto i)
                                               {
                    constexpr formulation::TermSpec term = spec.terms[evaluation.first_term + i];
                    if (!used(term))
                        return;
                    if constexpr (term.trial_field < 0)
                        rhs += kernel<term.kernel>().vector;
                    else
                    {
                        // The kernel applied to the solution at the quadrature points, without its matrix.
                        constexpr formulation::KernelSpec ks = spec.kernels[term.kernel];
                        constexpr formulation::NablaSpec tn = spec.nablas[static_cast<std::size_t>(ks.test)];
                        constexpr formulation::NablaSpec un = spec.nablas[static_cast<std::size_t>(ks.trial)];
                        using C = cpp_t<ks.coefficient.type>;
                        using Test = cpp_t<tn.output>;
                        using Trial = cpp_t<un.output>;
                        const TypedCoefficient<C> c = coefficient<C>(ks.coefficient);
                        kernel::applied(nabla(ks.test).template get<Test>(), c.view(), nabla(ks.trial).template get<Trial>(),
                                        trial_solution(term), weights_, rhs);
                    } });

                Vector x;
                mass_factors<mass_kernel>().cholesky_solve(rhs, x);
                for (std::size_t r = 0; r < ids.size(); ++r)
                    values(ids[r]) = x(static_cast<unsigned int>(r));
            }

            /** @brief Component @p c of the field of a value. */
            template <typename T>
            static Real component_of(const T &value, unsigned int c)
            {
                if constexpr (std::is_same_v<T, RealVectorValue>)
                    return value(c);
                else
                    return value;
            }

            /** @brief Rows of an element field: the volume mean of each term (mean_weights()), by number of dt() on the trial. */
            template <std::size_t E>
            void average_operator(const ElementType::FieldRange &ids, std::array<Matrix, 3> &F, Vector &f)
            {
                constexpr formulation::EvaluationSpec evaluation = spec.evaluations[E];
                constexpr formulation::FieldSpec field = spec.fields[static_cast<std::size_t>(evaluation.field)];
                using Out = cpp_t<field.value>;
                const std::vector<Real> &w = mean_weights();
                static_for<evaluation.n_terms>([&](auto i)
                                               {
                    constexpr formulation::TermSpec term = spec.terms[evaluation.first_term + i];
                    if (!used(term))
                        return;
                    const KernelResult &k = kernel<term.kernel>();
                    for (std::size_t c = 0; c < ids.size(); ++c)
                    {
                        const unsigned int component = field.components[c];
                        if constexpr (term.trial_field < 0)
                        {
                            for (std::size_t qp = 0; qp < nq_; ++qp)
                                f(ids[c]) += w[qp] * component_of(k.values[qp].template get<Out>(), component);
                        }
                        else
                        {
                            const auto &table = k.table.template get<Out>();
                            const auto &columns = in_.dofs->field[static_cast<std::size_t>(term.trial_field)];
                            Matrix &target = F[term.trial_dt];
                            for (std::size_t qp = 0; qp < nq_; ++qp)
                                for (std::size_t j = 0; j < columns.size(); ++j)
                                    target(ids[c], columns[j]) += w[qp] * component_of(table[qp][j], component);
                        }
                    } });
            }

            /**
             * @brief Rows of a nodal field: F_k = M^-1 Y_k and f = M^-1 y. The terms are first summed by number of dt()
             * into Y_k, one column per DOF of the element (the columns of each term's trial field), and the source
             * vectors into y; then one solve per column a term touches, and one for y.
             */
            template <std::size_t E>
            void projection_operator(const ElementType::FieldRange &ids, std::array<Matrix, 3> &F, Vector &f)
            {
                constexpr formulation::EvaluationSpec evaluation = spec.evaluations[E];
                const unsigned int n_rows = static_cast<unsigned int>(ids.size());
                std::array<Matrix, 3> Y;
                std::array<std::vector<bool>, 3> touched;
                Vector y;
                static_for<evaluation.n_terms>([&](auto i)
                                               {
                    constexpr formulation::TermSpec term = spec.terms[evaluation.first_term + i];
                    if (!used(term))
                        return;
                    const KernelResult &k = kernel<term.kernel>();
                    if constexpr (term.trial_field < 0)
                    {
                        if (y.size() == 0)
                            y.resize(n_rows);
                        y += k.vector;
                    }
                    else
                    {
                        Matrix &sum = Y[term.trial_dt];
                        std::vector<bool> &columns_used = touched[term.trial_dt];
                        if (sum.m() == 0)
                        {
                            sum.resize(n_rows, in_.dofs->n_dofs);
                            columns_used.assign(in_.dofs->n_dofs, false);
                        }
                        const auto &columns = in_.dofs->field[static_cast<std::size_t>(term.trial_field)];
                        for (unsigned int j = 0; j < k.matrix.n(); ++j)
                        {
                            columns_used[columns[j]] = true;
                            for (unsigned int r = 0; r < n_rows; ++r)
                                sum(r, columns[j]) += k.matrix(r, j);
                        }
                    } });

                Matrix &mass = mass_factors<static_cast<std::size_t>(evaluation.mass_kernel)>();
                Vector column(n_rows);
                Vector x;
                for (std::size_t d = 0; d < Y.size(); ++d)
                    for (unsigned int j = 0; j < Y[d].n(); ++j)
                    {
                        if (!touched[d][j])
                            continue;
                        for (unsigned int r = 0; r < n_rows; ++r)
                            column(r) = Y[d](r, j);
                        mass.cholesky_solve(column, x);
                        for (unsigned int r = 0; r < n_rows; ++r)
                            F[d](ids[r], j) += x(r);
                    }
                if (y.size() > 0)
                {
                    mass.cholesky_solve(y, x);
                    for (unsigned int r = 0; r < n_rows; ++r)
                        f(ids[r]) += x(r);
                }
            }

            /**
             * @brief Weights of the element mean of a field: w_q / sum of w, w the integration weights (JxW, times the
             * out-of-plane length of 2D elements: the mean over the volume the element stands for).
             */
            const std::vector<Real> &mean_weights()
            {
                if (mean_weights_.size() != nq_)
                {
                    Real volume = 0.0;
                    for (const Real w : weights_)
                        volume += w;
                    mean_weights_.resize(nq_);
                    for (std::size_t qp = 0; qp < nq_; ++qp)
                        mean_weights_[qp] = weights_[qp] / volume;
                }
                return mean_weights_;
            }

            /**
             * @brief The mass matrix of kernel @p K, factorized by its first solve (DenseMatrix keeps the factors): one
             * factorization per element for all the fields projected onto the same FE. A copy, since the kernel matrix
             * may also be a block term.
             */
            template <std::size_t K>
            Matrix &mass_factors()
            {
                if (!factors_done_[K])
                {
                    factors_[K] = kernel<K>().matrix;
                    factors_done_[K] = true;
                }
                return factors_[K];
            }

            const Inputs in_;
            const ElementType &type_;
            std::size_t nq_ = 0;
            std::vector<Real> weights_;
            std::vector<Real> length_;
            std::vector<Real> mean_weights_;
            std::array<ValueTable, kNablas> tables_{};
            std::array<bool, kNablas> table_done_{};
            std::array<KernelResult, kKernels> kernels_{};
            std::array<Matrix, kKernels> factors_{};
            std::array<bool, kKernels> factors_done_{};
            std::array<LawState, kLaws> laws_{};
        };

        template <std::size_t... I>
        std::unique_ptr<ElementEvaluator::Implementation> make_implementation(const ElementEvaluator::Inputs &inputs,
                                                                             std::index_sequence<I...>)
        {
            using Maker = std::unique_ptr<ElementEvaluator::Implementation> (*)(const ElementEvaluator::Inputs &);
            static constexpr Maker makers[] = {[](const ElementEvaluator::Inputs &in) -> std::unique_ptr<ElementEvaluator::Implementation>
                                               { return std::make_unique<Evaluator<static_cast<registry::ElementTypeId>(I)>>(in); }...};
            return makers[static_cast<std::size_t>(inputs.type->id())](inputs);
        }
    } // namespace

    ElementEvaluator::ElementEvaluator(const Inputs &inputs)
        : impl_(make_implementation(inputs, std::make_index_sequence<registry::kNumElementTypeId>{}))
    {
    }

    ElementEvaluator::~ElementEvaluator() = default;
} // namespace pytec3d

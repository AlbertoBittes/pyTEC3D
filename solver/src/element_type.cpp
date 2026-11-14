/**
 * @file element_type.cpp
 * @brief ElementType: FEs, quadrature and nablas of a spec; reinit and element matrix indices.
 */
#include "element_type.h"

#include <algorithm>
#include <atomic>
#include <unordered_map>


#include "formulations.h"
#include "logging.h"

namespace pytec3d
{
    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementType: " + message);
        }

        constexpr std::size_t kNumOperations = static_cast<std::size_t>(Operation::Parameters) + 1;

        std::atomic<std::size_t> next_instance{1};

        std::size_t needed_index(registry::Domain domain, Operation operation, bool nonlinear)
        {
            return (static_cast<std::size_t>(domain) * kNumOperations + static_cast<std::size_t>(operation)) * 2 + (nonlinear ? 1 : 0);
        }
    } // namespace

    bool ElementType::dynamic(const formulation::EvaluationSpec &evaluation, registry::Domain domain, registry::LawMask laws) const
    {
        if (laws == 0)
            return false;
        for (unsigned int t = evaluation.first_term; t < evaluation.first_term + evaluation.n_terms; ++t)
        {
            const formulation::TermSpec &term = spec_->terms[t];
            if (term.law >= 0 && active(term, domain) && selected(term, laws))
                return true;
        }
        return false;
    }

    bool ElementType::has_dynamic_fields(FieldKind kind, registry::Domain domain, registry::LawMask laws) const
    {
        const auto wanted = kind == FieldKind::Nodal ? formulation::EvaluationKind::L2Projection : formulation::EvaluationKind::Average;
        for (const formulation::EvaluationSpec &evaluation : spec_->evaluations)
            if (evaluation.kind == wanted && active(spec_->fields[static_cast<std::size_t>(evaluation.field)], domain) &&
                dynamic(evaluation, domain, laws))
                return true;
        return false;
    }

    bool ElementType::uses(const formulation::EvaluationSpec &evaluation,
                           const formulation::TermSpec &term,
                           Operation operation,
                           registry::LawMask laws,
                           registry::Domain domain) const
    {
        using formulation::EvaluationKind;
        using formulation::Stage;
        if (!active(term, domain) || !selected(term, laws))
            return false;
        const bool field = evaluation.kind == EvaluationKind::L2Projection || evaluation.kind == EvaluationKind::Average;
        switch (operation)
        {
        case Operation::Matrices:
            // The operators of the static fields, stored with the element matrices: the terms with a trial only.
            if (field)
                return !dynamic(evaluation, domain, laws) && term.trial_field >= 0;
            return evaluation.kind == EvaluationKind::Block && term.stage == Stage::Linear;
        case Operation::Rhs:
            return evaluation.kind == EvaluationKind::Rhs && term.stage == Stage::Linear;
        case Operation::Nonlinear:
            return term.stage == Stage::Tangent || term.stage == Stage::Residual || term.stage == Stage::RateResidual;
        // At output: the dynamic fields, and the terms without trial of every field (e.g. Jt = Js), which read
        // the current loads.
        case Operation::NodalFields:
            return evaluation.kind == EvaluationKind::L2Projection &&
                   (dynamic(evaluation, domain, laws) || term.trial_field < 0);
        case Operation::ElementFields:
            return evaluation.kind == EvaluationKind::Average &&
                   (dynamic(evaluation, domain, laws) || term.trial_field < 0);
        case Operation::Parameters:
            // The law variants are evaluated from the solution; their parameters are read there.
            return evaluation.kind != EvaluationKind::Quadrature && term.law < 0;
        }
        return false;
    }

    ElementType::ElementType(registry::ElementTypeId id, const ElementTypeOptions &options)
        : spec_(&formulation::spec(id)), options_(options), instance_(next_instance++)
    {
        const std::string context = "element type '" + std::string(name()) + "'";
        if (options.order < 1)
            fail(context + ": order must be >= 1.");
        if (options.axisymmetric && dimension() != 2)
            fail(context + ": axisymmetric requires a 2D element type.");
        if (spec_->fes.empty())
            fail(context + " has no finite element.");

        // One Gauss rule for every FE, so that all the nablas share the quadrature points: the highest default order
        // of the FEs, 2 p + 1 (libMesh's FEType::default_quadrature_order()).
        for (const formulation::FESpec &fe : spec_->fes)
        {
            quadrature_order_ = std::max(quadrature_order_, 2 * order(fe.order) + 1);
            FEState state;
            state.family = fe.family;
            state.vdim = (fe.family == formulation::FEFamily::LAGRANGE_VEC) ? dimension() : 1;
            fes_.push_back(state);
        }

        nablas_.reserve(spec_->nablas.size());
        for (const formulation::NablaSpec &nabla : spec_->nablas)
        {
            nablas_.emplace_back(nabla, dimension(), options.axisymmetric);
            nablas_.back().request(fes_[nabla.fe]);
        }

        for (const formulation::LawSpec &law : spec_->laws)
            laws_ |= registry::bit(law.id);
        for (const formulation::TermSpec &term : spec_->terms)
        {
            if (term.stage != formulation::Stage::Linear && term.trial_dt > 0)
                law_rates_ = true;
            if (term.stage == formulation::Stage::RateResidual)
                rate_residuals_ = true;
        }

        // The operators of the static fields: the terms of the field evaluations without law.
        for (const formulation::EvaluationSpec &evaluation : spec_->evaluations)
        {
            const bool nodal = evaluation.kind == formulation::EvaluationKind::L2Projection;
            if (!nodal && evaluation.kind != formulation::EvaluationKind::Average)
                continue;
            for (std::size_t d = 0; d < registry::kNumDomain; ++d)
                for (unsigned int t = evaluation.first_term; t < evaluation.first_term + evaluation.n_terms; ++t)
                {
                    const formulation::TermSpec &term = spec_->terms[t];
                    if (active(term, static_cast<registry::Domain>(d)) && selected(term, 0))
                        field_operators_[nodal ? 0 : 1][d] |= term.trial_field < 0 ? 8u : 1u << term.trial_dt;
                }
        }

        // The FEs each operation needs: those of the nablas of the kernels of the terms it computes. With laws, those
        // of every combination of them: the terms of the laws selected and the terms of the others.
        needed_fes_.resize(registry::kNumDomain * kNumOperations * 2);
        for (std::size_t d = 0; d < registry::kNumDomain; ++d)
            for (std::size_t o = 0; o < kNumOperations; ++o)
                for (const bool nonlinear : {false, true})
                {
                    const auto domain = static_cast<registry::Domain>(d);
                    const auto operation = static_cast<Operation>(o);
                    if (operation == Operation::Parameters)
                        continue;
                    std::vector<std::size_t> &needed = needed_fes_[needed_index(domain, operation, nonlinear)];
                    bool computes = false;
                    const auto add_kernel = [&](std::size_t k)
                    {
                        for (const int nabla : {spec_->kernels[k].trial, spec_->kernels[k].test})
                        {
                            if (nabla < 0)
                                continue;
                            const std::size_t fe = spec_->nablas[static_cast<std::size_t>(nabla)].fe;
                            if (std::find(needed.begin(), needed.end(), fe) == needed.end())
                                needed.push_back(fe);
                        }
                    };
                    for (const formulation::EvaluationSpec &evaluation : spec_->evaluations)
                        for (unsigned int t = evaluation.first_term; t < evaluation.first_term + evaluation.n_terms; ++t)
                        {
                            const formulation::TermSpec &term = spec_->terms[t];
                            if (!uses(evaluation, term, operation, 0, domain) &&
                                !(nonlinear && uses(evaluation, term, operation, laws_, domain)))
                                continue;
                            computes = true;
                            add_kernel(term.kernel);
                            if (evaluation.mass_kernel >= 0)
                                add_kernel(static_cast<std::size_t>(evaluation.mass_kernel));

                            // A law output is computed from the law's input.
                            const formulation::CoefficientSpec &c = spec_->kernels[term.kernel].coefficient;
                            if (!nonlinear || c.kind != formulation::CoefficientKind::LawOutput)
                                continue;
                            const auto law_id = registry::info(static_cast<registry::LawOutputId>(c.index)).law;
                            for (const formulation::LawSpec &law : spec_->laws)
                            {
                                if (law.id != law_id)
                                    continue;
                                for (const int evaluation : {static_cast<int>(law.input), static_cast<int>(law.inputs[0]),
                                                             static_cast<int>(law.inputs[1])})
                                {
                                    if (evaluation < 0)
                                        continue;
                                    const formulation::EvaluationSpec &input = spec_->evaluations[static_cast<std::size_t>(evaluation)];
                                    for (unsigned int i = input.first_term; i < input.first_term + input.n_terms; ++i)
                                        if (active(spec_->terms[i], domain))
                                            add_kernel(spec_->terms[i].kernel);
                                }
                            }
                        }
                    // Terms without nabla (a coefficient alone) still need the quadrature points.
                    if (computes && needed.empty())
                        needed.push_back(0);
                }
    }

    SystemVariable ElementType::variable(const FieldSpec &field) const
    {
        return {field.variable, field.family, order(field.order)};
    }

    bool ElementType::active(const FieldSpec &field, registry::Domain domain)
    {
        return (field.domains & registry::bit(domain)) != 0;
    }

    bool ElementType::active(const formulation::TermSpec &term, registry::Domain domain)
    {
        return (term.domains & registry::bit(domain)) != 0;
    }

    bool ElementType::selected(const formulation::TermSpec &term, registry::LawMask laws)
    {
        return term.law >= 0 ? (laws & registry::bit(static_cast<registry::LawId>(term.law))) != 0
                             : (term.replaced_by & laws) == 0;
    }

    unsigned int ElementType::field_operators(FieldKind kind, registry::Domain domain) const
    {
        return field_operators_[kind == FieldKind::Nodal ? 0 : 1][static_cast<std::size_t>(domain)];
    }

    int ElementType::order(formulation::OrderRule rule) const
    {
        switch (rule)
        {
        case formulation::OrderRule::Option:
            return static_cast<int>(options_.order);
        case formulation::OrderRule::First:
            return 1;
        case formulation::OrderRule::Constant:
            return 0;
        }
        fail("unknown order rule.");
    }


    bool ElementType::has_stiffness(std::size_t field, registry::Domain domain) const
    {
        for (const formulation::EvaluationSpec &evaluation : spec_->evaluations)
        {
            if (evaluation.kind != formulation::EvaluationKind::Block)
                continue;
            for (unsigned int t = evaluation.first_term; t < evaluation.first_term + evaluation.n_terms; ++t)
            {
                const formulation::TermSpec &term = spec_->terms[t];
                if (term.trial_dt == 0 && term.test_field == static_cast<int>(field) && active(term, domain))
                    return true;
            }
        }
        return false;
    }

    const ElementGeometry *ElementType::reinit(const MeshElement &elem, registry::Domain domain, Operation operation,
                                               registry::LawMask laws) const
    {
        if (operation == Operation::Parameters)
            return nullptr;
        const auto &needed = needed_fes_.at(needed_index(domain, operation, laws != 0));
        if (needed.empty())
            return nullptr;

        const mfem::Geometry::Type geometry = elem.mesh->mesh().GetElementGeometry(elem.index);
        const fem::QuadratureData &quadrature = *elem.quadrature;
        const mfem::IntegrationRule &rule = quadrature.rule(elem.index);
        const int nq = rule.GetNPoints();
        const int dim = static_cast<int>(dimension());
        const double *X = quadrature.X(elem.index);
        const double *J = quadrature.J(elem.index);
        const double *det = quadrature.det(elem.index);

        Scratch &work = scratch();
        work.geometry.JxW.resize(static_cast<std::size_t>(nq));
        work.geometry.xyz.resize(static_cast<std::size_t>(nq));
        // J^-T at each point (dim x dim, in the 3 x 3 layout of QuadratureData).
        std::vector<double> inverse_transpose(static_cast<std::size_t>(9 * nq), 0.0);
        for (int q = 0; q < nq; ++q)
        {
            work.geometry.JxW[static_cast<std::size_t>(q)] = rule.IntPoint(q).weight * std::abs(det[q]);
            work.geometry.xyz[static_cast<std::size_t>(q)] = la::Point(X[3 * q], X[3 * q + 1], X[3 * q + 2]);
            const double *Jq = J + 9 * q;
            double *K = inverse_transpose.data() + 9 * q;
            const double d = det[q];
            if (dim == 3)
            {
                // (J^-1)^T(i, j) = cofactor(J)(i, j) / det.
                K[0] = (Jq[4] * Jq[8] - Jq[5] * Jq[7]) / d;
                K[1] = (Jq[5] * Jq[6] - Jq[3] * Jq[8]) / d;
                K[2] = (Jq[3] * Jq[7] - Jq[4] * Jq[6]) / d;
                K[3] = (Jq[2] * Jq[7] - Jq[1] * Jq[8]) / d;
                K[4] = (Jq[0] * Jq[8] - Jq[2] * Jq[6]) / d;
                K[5] = (Jq[1] * Jq[6] - Jq[0] * Jq[7]) / d;
                K[6] = (Jq[1] * Jq[5] - Jq[2] * Jq[4]) / d;
                K[7] = (Jq[2] * Jq[3] - Jq[0] * Jq[5]) / d;
                K[8] = (Jq[0] * Jq[4] - Jq[1] * Jq[3]) / d;
            }
            else
            {
                K[0] = Jq[4] / d;
                K[1] = -Jq[3] / d;
                K[3] = -Jq[1] / d;
                K[4] = Jq[0] / d;
            }
        }

        for (const std::size_t index : needed)
        {
            FEState &state = work.fes[index];
            const formulation::FESpec &spec = spec_->fes[index];
            const mfem::FiniteElement *fe = elem.catalog->element(spec.family, order(spec.order), geometry);
            const fem::ReferenceTable &table = elem.tables->get(*fe, rule);
            const int n = table.n_dofs;
            if (state.nedelec())
            {
                state.vphi.assign(static_cast<std::size_t>(nq), std::vector<RealVectorValue>(static_cast<std::size_t>(n)));
                if (state.need_curl)
                    state.curl.assign(static_cast<std::size_t>(nq), std::vector<RealVectorValue>(static_cast<std::size_t>(n)));
                for (int q = 0; q < nq; ++q)
                {
                    const double *K = inverse_transpose.data() + 9 * q;
                    const double *Jq = J + 9 * q;
                    for (int i = 0; i < n; ++i)
                    {
                        // Covariant: phi = J^-T phi_ref.
                        const double *b = table.B.data() + static_cast<std::size_t>((q * n + i) * dim);
                        RealVectorValue &v = state.vphi[static_cast<std::size_t>(q)][static_cast<std::size_t>(i)];
                        for (int r = 0; r < dim; ++r)
                            for (int k = 0; k < dim; ++k)
                                v(static_cast<unsigned int>(r)) += K[3 * r + k] * b[k];
                        if (!state.need_curl)
                            continue;
                        const double *c = table.G.data() + static_cast<std::size_t>((q * n + i) * table.curl_dim);
                        RealVectorValue &curl = state.curl[static_cast<std::size_t>(q)][static_cast<std::size_t>(i)];
                        if (dim == 3)
                        {
                            // curl phi = J curl_ref / det J.
                            for (int r = 0; r < 3; ++r)
                                for (int k = 0; k < 3; ++k)
                                    curl(static_cast<unsigned int>(r)) += Jq[3 * r + k] * c[k] / det[q];
                        }
                        else
                            curl(2) = c[0] / det[q];
                    }
                }
            }
            else
            {
                state.phi.assign(static_cast<std::size_t>(nq), std::vector<Real>(static_cast<std::size_t>(n)));
                if (state.need_dphi)
                    state.dphi.assign(static_cast<std::size_t>(nq), std::vector<RealVectorValue>(static_cast<std::size_t>(n)));
                for (int q = 0; q < nq; ++q)
                {
                    const double *K = inverse_transpose.data() + 9 * q;
                    for (int i = 0; i < n; ++i)
                    {
                        state.phi[static_cast<std::size_t>(q)][static_cast<std::size_t>(i)] =
                            table.B[static_cast<std::size_t>(q * n + i)];
                        if (!state.need_dphi)
                            continue;
                        // grad phi = J^-T grad_ref phi.
                        const double *g = table.G.data() + static_cast<std::size_t>((q * n + i) * dim);
                        RealVectorValue &v = state.dphi[static_cast<std::size_t>(q)][static_cast<std::size_t>(i)];
                        for (int r = 0; r < dim; ++r)
                            for (int k = 0; k < dim; ++k)
                                v(static_cast<unsigned int>(r)) += K[3 * r + k] * g[k];
                    }
                }
            }
        }
        return &work.geometry;
    }

    ElementType::Scratch &ElementType::scratch() const
    {
        thread_local std::unordered_map<std::size_t, Scratch> scratches;
        auto it = scratches.find(instance_);
        if (it == scratches.end())
            it = scratches.emplace(instance_, Scratch{fes_, ElementGeometry{}}).first;
        return it->second;
    }

    const std::vector<unsigned int> &ElementType::variable_numbers(const fem::SystemLayout &system) const
    {
        std::lock_guard<std::mutex> lock(variable_numbers_mutex_);
        auto it = variable_numbers_.find(&system);
        if (it == variable_numbers_.end())
        {
            std::vector<unsigned int> numbers(spec_->fields.size(), la::invalid_uint);
            for (std::size_t index = 0; index < spec_->fields.size(); ++index)
            {
                const std::string label(registry::info(spec_->fields[index].variable).label);
                const int number = system.variable_number(label);
                if (number >= 0)
                    numbers[index] = static_cast<unsigned int>(number);
            }
            it = variable_numbers_.emplace(&system, std::move(numbers)).first;
        }
        return it->second;
    }

    ElementType::FieldIndices ElementType::element_indices(const fem::SystemLayout &system,
                                                           const MeshElement &elem,
                                                           FieldKind kind,
                                                           registry::Domain domain) const
    {
        // The element vector of the system: the variables active on the subdomain, in variable order, each a
        // contiguous block. The buffer is reused across calls on the thread.
        const unsigned int n_vars = system.n_vars();
        thread_local std::vector<unsigned int> offset;
        offset.assign(n_vars + 1, 0);
        for (unsigned int v = 0; v < n_vars; ++v)
            offset[v + 1] = offset[v] + system.n_element_dofs(v, elem.index);

        FieldIndices out;
        out.n_dofs = offset[n_vars];
        out.field.resize(spec_->fields.size());

        const auto &numbers = variable_numbers(system);
        std::size_t covered = 0;
        for (std::size_t index = 0; index < spec_->fields.size(); ++index)
        {
            const FieldSpec &field = spec_->fields[index];
            if (field.kind != kind || !active(field, domain))
                continue;
            const unsigned int var = numbers[index];
            if (var == la::invalid_uint)
                fail("the variable of '" + std::string(registry::name(field.id)) + "' is not in the system.");

            FieldRange &range = out.field[index];
            range.start = offset[var];
            range.count = offset[var + 1] - offset[var];
            covered += range.count;
        }

        if (covered != out.n_dofs)
            fail("the fields active in the domain cover " + std::to_string(covered) + " of " + std::to_string(out.n_dofs) +
                 " element DOFs of the system.");
        return out;
    }
} // namespace pytec3d

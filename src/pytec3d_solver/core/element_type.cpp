/**
 * @file element_type.cpp
 * @brief ElementType: FEs, quadrature and nablas of a spec; reinit and element matrix indices.
 */
#include "element_type.h"

#include <algorithm>

#include <libmesh/dof_map.h>
#include <libmesh/fe_interface.h>
#include <libmesh/quadrature_gauss.h>

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
            if (field)
                return !dynamic(evaluation, domain, laws);
            return evaluation.kind == EvaluationKind::Block && term.stage == Stage::Linear;
        case Operation::Rhs:
            return evaluation.kind == EvaluationKind::Rhs && term.stage == Stage::Linear;
        case Operation::Nonlinear:
            return term.stage == Stage::Tangent || term.stage == Stage::Residual;
        case Operation::NodalFields:
            return evaluation.kind == EvaluationKind::L2Projection && dynamic(evaluation, domain, laws);
        case Operation::ElementFields:
            return evaluation.kind == EvaluationKind::Average && dynamic(evaluation, domain, laws);
        case Operation::Parameters:
            // The law variants are evaluated from the solution; their parameters are read there.
            return evaluation.kind != EvaluationKind::Quadrature && term.law < 0;
        }
        return false;
    }

    ElementType::ElementType(registry::ElementTypeId id, const ElementTypeOptions &options)
        : spec_(&formulation::spec(id)), options_(options)
    {
        const std::string context = "element type '" + std::string(name()) + "'";
        if (options.order < static_cast<unsigned int>(libMesh::FIRST))
            fail(context + ": order must be >= 1.");
        if (options.axisymmetric && dimension() != 2)
            fail(context + ": axisymmetric requires a 2D element type.");
        if (spec_->fes.empty())
            fail(context + " has no finite element.");

        // One Gauss rule for every FE, so that all the nablas share the quadrature points.
        unsigned int quadrature_order = 0;
        for (const formulation::FESpec &fe : spec_->fes)
        {
            const libMesh::FEType type(order(fe.order), fe.family);
            quadrature_order = std::max(quadrature_order, static_cast<unsigned int>(type.default_quadrature_order()));
        }
        quadrature_ = std::make_unique<libMesh::QGauss>(dimension(), static_cast<libMesh::Order>(quadrature_order));

        for (const formulation::FESpec &fe : spec_->fes)
        {
            fes_.push_back(libMesh::FEAbstract::build(dimension(), libMesh::FEType(order(fe.order), fe.family)));
            fes_.back()->attach_quadrature_rule(quadrature_.get());
        }
        // Any of them can give the geometry of an operation.
        for (auto &fe : fes_)
        {
            fe->get_JxW();
            fe->get_xyz();
        }

        nablas_.reserve(spec_->nablas.size());
        for (const formulation::NablaSpec &nabla : spec_->nablas)
            nablas_.emplace_back(nabla, *fes_[nabla.fe], dimension(), options.axisymmetric);

        for (const formulation::LawSpec &law : spec_->laws)
            laws_ |= registry::bit(law.id);
        for (const formulation::TermSpec &term : spec_->terms)
            if (term.stage != formulation::Stage::Linear && term.trial_dt > 0)
                law_rates_ = true;

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
                    std::vector<libMesh::FEAbstract *> &needed = needed_fes_[needed_index(domain, operation, nonlinear)];
                    bool computes = false;
                    const auto add_kernel = [&](std::size_t k)
                    {
                        for (const int nabla : {spec_->kernels[k].trial, spec_->kernels[k].test})
                        {
                            if (nabla < 0)
                                continue;
                            libMesh::FEAbstract *fe = fes_[spec_->nablas[static_cast<std::size_t>(nabla)].fe].get();
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
                                const formulation::EvaluationSpec &input = spec_->evaluations[static_cast<std::size_t>(law.input)];
                                for (unsigned int i = input.first_term; i < input.first_term + input.n_terms; ++i)
                                    if (active(spec_->terms[i], domain))
                                        add_kernel(spec_->terms[i].kernel);
                            }
                        }
                    // Terms without nabla (a coefficient alone) still need the quadrature points.
                    if (computes && needed.empty())
                        needed.push_back(fes_.front().get());
                }
    }

    libMesh::Order ElementType::order(formulation::OrderRule rule) const
    {
        switch (rule)
        {
        case formulation::OrderRule::Option:
            return static_cast<libMesh::Order>(options_.order);
        case formulation::OrderRule::First:
            return libMesh::FIRST;
        case formulation::OrderRule::Constant:
            return libMesh::CONSTANT;
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

    const libMesh::FEAbstract *ElementType::reinit(const libMesh::Elem &elem, registry::Domain domain, Operation operation,
                                                   registry::LawMask laws) const
    {
        if (operation == Operation::Parameters)
            return nullptr;
        const auto &needed = needed_fes_.at(needed_index(domain, operation, laws != 0));
        for (libMesh::FEAbstract *fe : needed)
            fe->reinit(&elem);
        return needed.empty() ? nullptr : needed.front();
    }

    const std::vector<unsigned int> &ElementType::variable_numbers(const libMesh::System &system) const
    {
        std::lock_guard<std::mutex> lock(variable_numbers_mutex_);
        auto it = variable_numbers_.find(system.number());
        if (it == variable_numbers_.end())
        {
            std::vector<unsigned int> numbers(spec_->fields.size(), libMesh::invalid_uint);
            for (std::size_t index = 0; index < spec_->fields.size(); ++index)
            {
                const std::string label(registry::info(spec_->fields[index].variable).label);
                if (system.has_variable(label))
                    numbers[index] = system.variable_number(label);
            }
            it = variable_numbers_.emplace(system.number(), std::move(numbers)).first;
        }
        return it->second;
    }

    ElementType::FieldIndices ElementType::element_indices(const libMesh::System &system,
                                                           const libMesh::Elem &elem,
                                                           FieldKind kind,
                                                           registry::Domain domain) const
    {
        // The element vector of the system: the variables active on the subdomain, in variable order, each a
        // contiguous block. The buffer is reused across calls on the thread.
        const unsigned int n_vars = system.n_vars();
        thread_local std::vector<unsigned int> offset;
        offset.assign(n_vars + 1, 0);
        for (unsigned int v = 0; v < n_vars; ++v)
        {
            const bool active_here = system.variable(v).active_on_subdomain(elem.subdomain_id());
            offset[v + 1] = offset[v] + (active_here ? libMesh::FEInterface::n_dofs(system.variable_type(v), &elem) : 0);
        }

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
            if (var == libMesh::invalid_uint)
                fail("the variable of '" + std::string(registry::name(field.id)) + "' is not in system '" + system.name() + "'.");

            FieldRange &range = out.field[index];
            range.start = offset[var];
            range.count = offset[var + 1] - offset[var];
            covered += range.count;
        }

        if (covered != out.n_dofs)
            fail("the fields active in the domain cover " + std::to_string(covered) + " of " + std::to_string(out.n_dofs) +
                 " element DOFs of system '" + system.name() + "'.");
        return out;
    }
} // namespace pytec3d

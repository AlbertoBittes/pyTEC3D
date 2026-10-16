/**
 * @file element_type.h
 * @brief An element type instantiated for its options: the libMesh FEs, quadrature and nablas of its spec.
 *
 * The formulation itself is the compiled spec (formulation_spec.h, generated from element_types.json). An
 * ElementType adds what depends on the options of the model data (order, axisymmetry):
 *  - the FEs, one per FESpec, at the order of the options (FIRST for NEDELEC_ONE), sharing one Gauss rule of
 *    the highest default order among them;
 *  - the nablas, bound to their FE;
 *  - the libMesh variables of each field, with family and order.
 * It reinit()s on an element only the FEs an operation needs, and tells where each field sits in the
 * element vectors of a system. Evaluating the kernels is ElementEvaluator's job.
 *
 * What an operation needs follows from the spec: the evaluations of the operation, their terms active in
 * the domain (with or without the laws selected on the element), the kernels of those terms, the nablas of the
 * kernels and, for a law variant, the terms of the law's input; the FEs are those of the nablas.
 */
#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string_view>
#include <vector>

#include <libmesh/elem.h>
#include <libmesh/fe.h>
#include <libmesh/point.h>
#include <libmesh/quadrature.h>
#include <libmesh/system.h>

#include "formulation_spec.h"
#include "nabla.h"
#include "registries.h"

namespace pytec3d
{
    using formulation::FieldKind;
    using formulation::FieldSpec;

    /** @brief What is computed on an element. */
    enum class Operation : unsigned char
    {
        Matrices,      ///< K, C, M and the operators of the static nodal and element fields (stored with them).
        Rhs,           ///< Right-hand side.
        Nonlinear,     ///< Jacobian and residual of the laws selected on an element.
        NodalFields,   ///< L2 projections of the dynamic nodal fields, at output.
        ElementFields, ///< Averages of the dynamic element fields, at output.
        Parameters,    ///< Only the coefficients of the terms (model data checks): no FE.
    };

    struct ElementTypeOptions
    {
        unsigned int order = 1;    ///< FE order of the Lagrange families.
        bool axisymmetric = false; ///< 2D elements only.
    };

    /** @brief The libMesh variable of a field. */
    struct LibMeshVariable
    {
        registry::VariableId id;
        libMesh::FEFamily family;
        libMesh::Order order;
    };

    class ElementType
    {
    public:
        ElementType(registry::ElementTypeId id, const ElementTypeOptions &options);
        ElementType(const ElementType &) = delete;
        ElementType &operator=(const ElementType &) = delete;

        const formulation::ElementTypeSpec &spec() const { return *spec_; }
        registry::ElementTypeId id() const { return spec_->id; }
        std::string_view name() const { return registry::name(spec_->id); }
        unsigned int dimension() const { return spec_->dimension; }
        const ElementTypeOptions &options() const { return options_; }

        const FieldSpec &field(std::size_t index) const { return spec_->fields[index]; }
        LibMeshVariable variable(const FieldSpec &field) const { return {field.variable, field.family, order(field.order)}; }

        static bool active(const FieldSpec &field, registry::Domain domain)
        {
            return (field.domains & registry::bit(domain)) != 0;
        }

        /** @brief True when the term applies in @p domain. */
        static bool active(const formulation::TermSpec &term, registry::Domain domain)
        {
            return (term.domains & registry::bit(domain)) != 0;
        }

        /**
         * @brief True when a term is used on an element where the laws @p laws are selected: a law variant where its
         * law is, any other term where no law replacing its coefficient is.
         */
        static bool selected(const formulation::TermSpec &term, registry::LawMask laws)
        {
            return term.law >= 0 ? (laws & registry::bit(static_cast<registry::LawId>(term.law))) != 0
                                 : (term.replaced_by & laws) == 0;
        }

        /**
         * @brief True when @p operation computes @p term of @p evaluation on an element in @p domain where the laws
         * @p laws are selected (selected()).
         *
         * Linear block and rhs terms go to Matrices and Rhs, the law terms (J, R) to Nonlinear. The terms of the
         * static nodal and element fields belong to Matrices, those of the dynamic ones to NodalFields / ElementFields.
         */
        bool uses(const formulation::EvaluationSpec &evaluation,
                  const formulation::TermSpec &term,
                  Operation operation,
                  registry::LawMask laws,
                  registry::Domain domain) const;

        /**
         * @brief True when the nodal or element field of @p evaluation depends on the solution through a law selected
         * on the element (a law variant of its terms). Static fields are linear operators on the solution, computed
         * once with the element matrices (field_operators); dynamic ones at every output.
         */
        bool dynamic(const formulation::EvaluationSpec &evaluation, registry::Domain domain, registry::LawMask laws) const;

        /** @brief True when a field of @p kind is dynamic in @p domain. */
        bool has_dynamic_fields(FieldKind kind, registry::Domain domain, registry::LawMask laws) const;

        /** @brief The laws bound to the element type (those replacing one of its parameters). */
        registry::LawMask laws() const { return laws_; }

        /** @brief True when a law term has a dt() trial: the residual needs the velocity or acceleration. */
        bool law_rates() const { return law_rates_; }

        /**
         * @brief The operators the static nodal or element fields can have in @p domain: bit k (k = 0, 1, 2) for
         * F_k, the terms with k dt() on the trial, bit 3 for the constant (terms without trial). Known from the
         * spec: the others are never allocated.
         */
        unsigned int field_operators(FieldKind kind, registry::Domain domain) const
        {
            return field_operators_[kind == FieldKind::Nodal ? 0 : 1][static_cast<std::size_t>(domain)];
        }

        const Nabla &nabla(std::size_t index) const { return nablas_[index]; }

        /**
         * @brief True when the DOF field @p field has a row in K in @p domain: it is the test of a block without dt()
         * on the trial (linear, secant or tangent). A field without one only has C or M rows, which vanish on a
         * static step.
         */
        bool has_stiffness(std::size_t field, registry::Domain domain) const;

        /**
         * @brief reinit()s on @p elem the FEs @p operation needs in @p domain.
         * @return The FE that gives JxW and the quadrature points; null when the operation computes nothing there.
         */
        const libMesh::FEAbstract *reinit(const libMesh::Elem &elem, registry::Domain domain, Operation operation,
                                          registry::LawMask laws) const;

        /**
         * @brief The positions of the values of one field in the element vector: the block of its variable, value i
         * (shape function i) at start + i.
         */
        struct FieldRange
        {
            unsigned int start = 0;
            unsigned int count = 0; ///< 0 for other kinds and inactive fields.

            unsigned int size() const { return count; }
            bool empty() const { return count == 0; }
            unsigned int operator[](std::size_t i) const { return start + static_cast<unsigned int>(i); }
        };

        /** @brief Element vector position of every shape function of every field of one kind. */
        struct FieldIndices
        {
            std::vector<FieldRange> field; ///< By field index; empty for other kinds and inactive fields.
            unsigned int n_dofs = 0;       ///< Size of the element vectors of the system.
        };

        /**
         * @brief Where, in the element vectors of @p system, each value of each field of @p kind sits.
         *
         * Follows the order of DofMap::dof_indices(): the variables active on the subdomain in system order,
         * each with FEInterface::n_dofs() DOFs (node by node, then on the element), so positions are prefix sums
         * of those counts, with no DOF lookup. The fields are those of @p kind active in @p domain; each is a
         * FieldRange, the block of its variable. The variable numbers of the fields are looked up once per system.
         */
        FieldIndices element_indices(const libMesh::System &system,
                                     const libMesh::Elem &elem,
                                     FieldKind kind,
                                     registry::Domain domain) const;

    private:
        libMesh::Order order(formulation::OrderRule rule) const;

        /** @brief Variable number in @p system of the variable of each field (by field index; libMesh::invalid_uint
         * if not in the system); looked up once. */
        const std::vector<unsigned int> &variable_numbers(const libMesh::System &system) const;

        const formulation::ElementTypeSpec *spec_;
        ElementTypeOptions options_;
        std::unique_ptr<libMesh::QBase> quadrature_;
        std::vector<std::unique_ptr<libMesh::FEAbstract>> fes_;
        std::vector<Nabla> nablas_;
        registry::LawMask laws_ = 0;
        bool law_rates_ = false;
        std::array<std::array<unsigned char, registry::kNumDomain>, 2> field_operators_{}; ///< Nodal, element; by domain.
        /// The FEs to reinit, by domain, operation and whether laws are selected; the first gives the geometry.
        std::vector<std::vector<libMesh::FEAbstract *>> needed_fes_;
        mutable std::map<unsigned int, std::vector<unsigned int>> variable_numbers_; ///< By system number.
        mutable std::mutex variable_numbers_mutex_;
    };
} // namespace pytec3d

/**
 * @file formulation_spec.h
 * @brief The compiled form of an element type of element_types.json: plain tables of enums and indices.
 *
 * generate_registries.py reads element_types.json at build time, checks it and writes one ElementTypeSpec
 * per element type (formulations.h). Nothing is parsed at run time: ElementType instantiates the libMesh
 * objects a spec asks for, ElementEvaluator runs its evaluations. Adding a formulation means editing the
 * JSON and rebuilding the module.
 *
 * An element type is made of
 *  - Fields: DOF unknowns, nodal fields (L2 projections) and element fields (quadrature averages).
 *  - FEs: one libMesh FE per distinct family, after merging the families with the same shape functions
 *    (LAGRANGE and L2_LAGRANGE, LAGRANGE_VEC and L2_LAGRANGE_VEC). Element fields have none.
 *  - Nablas: an operator (identity, grad, curl, div) on one FE; the shape functions it gives at each
 *    quadrature point are Real, RealVectorValue or RealTensorValue (output). Duplicates are merged.
 *  - Kernels: test^T * coefficient * trial over the quadrature points, with one or two nablas. Duplicates
 *    are merged, so a kernel used by several terms is computed once per element.
 *  - Terms: a kernel placed between two fields, active in some domains.
 *  - Evaluations: what the terms are combined into (EvaluationKind).
 *  - Laws (laws.json): nonlinear coefficients computed at each quadrature point from a field and a curve, in place
 *    of the parameters they replace. A term whose coefficient a law replaces has variants with the law's outputs
 *    (TermSpec::law); where the model data selects the law for an element, those are used instead of it
 *    (TermSpec::replaced_by).
 */
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

#include <libmesh/enum_fe_family.h>
#include <libmesh/libmesh_common.h>

#include "parameter_value.h"
#include "registries.h"

namespace pytec3d::formulation
{
    using ValueType = ParameterType;

    /** @brief A read-only array of a spec: pointer and size of a generated table. */
    template <typename T>
    struct Span
    {
        const T *data = nullptr;
        std::size_t count = 0;

        constexpr const T *begin() const { return data; }
        constexpr const T *end() const { return data + count; }
        constexpr std::size_t size() const { return count; }
        constexpr bool empty() const { return count == 0; }
        constexpr const T &operator[](std::size_t i) const { return data[i]; }
    };

    enum class FieldKind : unsigned char
    {
        Dof,        ///< Unknown of the DOF system.
        Nodal,      ///< Post-processed, L2-projected onto its FE.
        Element,    ///< Post-processed, averaged over the quadrature points (one value per element).
    };

    /** @brief Where an FE takes its order from. */
    enum class OrderRule : unsigned char
    {
        Option, ///< The element type option (option2): the Lagrange families.
        First,    ///< Always first order: NEDELEC_ONE.
        Constant, ///< Element fields (MONOMIAL).
    };

    enum class Operator : unsigned char
    {
        Identity,
        Grad, ///< On an out-of-plane scalar of a 2D element (active component z): -phi / L e_z, L the out-of-plane length.
        Curl, ///< On a scalar of a 2D element: the curl of the out-of-plane vector phi e_z.
        Div,
    };

    enum class CoefficientKind : unsigned char
    {
        One,       ///< No coefficient: 1.
        Parameter, ///< A registry::ParameterId, constant over the element.
        LawOutput, ///< A registry::LawOutputId, one value per quadrature point (in the law variants of the terms).
    };

    /** @brief Where a dof term goes. */
    enum class Stage : unsigned char
    {
        Linear,   ///< K, C (dt() trial), M (dt(dt()) trial) or F; or a nodal / element field term. Assembled once, cached.
        Tangent,  ///< J (law): + c_d block, d the dt() on the trial, at the current solution.
        Residual, ///< R (law): + block * d^d u / dt^d, or + the vector, at the current solution.
    };

    enum class EvaluationKind : unsigned char
    {
        Block,        ///< K and C, residual and Jacobian: terms with test and trial.
        Rhs,          ///< Right-hand side (or residual): terms with a test and no trial.
        Quadrature,   ///< A field at the quadrature points, the input of a law: terms without test.
        L2Projection, ///< A nodal field: M^-1 * terms with the field's shape functions as test.
        Average,      ///< An element field: mean over the element volume (JxW-weighted) of terms without test.
    };

    struct FieldSpec
    {
        registry::FieldId id;
        FieldKind kind;
        ValueType value;                                 ///< Of the field at a point: Scalar, or Vector (also for a single component).
        unsigned char n_components;                      ///< Values per point: 1 for a scalar or a single component.
        std::array<unsigned char, 3> components;         ///< Physical components (0=x, 1=y, 2=z) of the values, first n_components.
        signed char fe;                                  ///< Index of its FE; -1 for element fields.
        libMesh::FEFamily family;                        ///< Family of its libMesh variable.
        OrderRule order;                                 ///< Of its libMesh variable; element fields are CONSTANT.
        registry::VariableId variable;                   ///< Its libMesh variable (vector, scalar, or one component).
        registry::DomainMask domains;                    ///< Where it exists.
        registry::FieldId nodal_load;                    ///< Point array loading a DOF field on its nodes; COUNT if none.
    };

    struct FESpec
    {
        libMesh::FEFamily family;
        bool vector; ///< FEVectorBase (RealVectorValue shape functions) or FEBase (Real).
        OrderRule order;
    };

    struct NablaSpec
    {
        unsigned char fe;
        Operator op;
        signed char component; ///< A scalar FE holding component 0..2 of a vector field; -1 otherwise.
        bool out_of_plane;     ///< A scalar of a 2D element varying along z only (active component z): a fake 3D scalar.
        ValueType output;
    };

    struct CoefficientSpec
    {
        CoefficientKind kind;
        unsigned short index; ///< ParameterId or LawOutputId, by kind.
        ValueType type;
        libMesh::Real sign;
    };

    /** @brief sum over the quadrature points of test^T * coefficient * trial; times JxW only when there is a test. */
    struct KernelSpec
    {
        signed char trial; ///< Nabla index; -1 for none.
        signed char test;  ///< Nabla index; -1 for none.
        CoefficientSpec coefficient;
        ValueType output; ///< Without test: the type of coefficient * trial at a quadrature point.
    };

    struct TermSpec
    {
        unsigned short kernel;
        signed char trial_field; ///< Field index; -1 for none.
        signed char test_field;  ///< Field index; -1 for none.
        unsigned char trial_dt;  ///< Time derivatives on the trial: 0 (K), 1 (C, velocity), 2 (M, acceleration).
        Stage stage;
        signed char law;              ///< LawId of a law variant: used where the law is selected; -1 for the others.
        registry::LawMask replaced_by; ///< Laws replacing the coefficient: the term is not used where one is selected.
        registry::DomainMask domains; ///< Where it applies, already intersected with the domains of its fields.
    };

    struct EvaluationSpec
    {
        EvaluationKind kind;
        signed char field;          ///< The evaluated field; -1 for Block and Rhs.
        signed short mass_kernel;   ///< L2Projection: kernel of the field's mass matrix; -1 otherwise.
        unsigned short first_term;
        unsigned short n_terms;
    };

    /** @brief A law of laws.json bound to an element type: its input and the outputs it has there. */
    struct LawSpec
    {
        registry::LawId id;
        signed short input;                  ///< Quadrature evaluation of the field it reads.
        Span<registry::LawOutputId> outputs; ///< Secant (and tangent) of each parameter it replaces, in the order of evaluate_law().
    };

    struct ElementTypeSpec
    {
        registry::ElementTypeId id;
        unsigned int dimension;
        unsigned int time_order; ///< Highest time derivative: 0 static, 1 (K, C), 2 (K, C, M).
        signed short out_of_plane_length; ///< ParameterId of the out-of-plane length of 2D elements; -1: L = 1.
        Span<FieldSpec> fields;
        Span<FESpec> fes;
        Span<NablaSpec> nablas;
        Span<KernelSpec> kernels;
        Span<TermSpec> terms;
        Span<EvaluationSpec> evaluations;
        Span<LawSpec> laws;
    };
} // namespace pytec3d::formulation

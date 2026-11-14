/**
 * @file element_builder.h
 * @brief Per-element assembly: the matrices, vectors and fields of one element from the kernels of its element type.
 *
 * ElementBuilder collects what an element needs (material and real constants, local-to-global rotation,
 * curves, loads of the model arrays, local solution), finds where each field sits in the element vectors of the systems
 * (ElementType::element_indices) and has an ElementEvaluator run the evaluations of the element type:
 *   - Stiffness (Ke), time-derivative (Ce) and, for second-order element types, mass (Me) matrices of the DOF system;
 *   - the source (RHS) vector, and the tangent matrix / residual of the laws selected on the element;
 *   - the values of its nodal and element fields.
 * The element vectors are in the shape functions of the element (DofSpace::element_values()); the DOF system moves them
 * to the global DOFs. The formulations are those of element_types.json (formulations.h, generated).
 */
#pragma once

#include <array>
#include <vector>

#include "dof_space.h"
#include "element_evaluator.h"
#include "element_parameters.h"
#include "element_type.h"
#include "la.h"
#include "model_data.h"
#include "model_element_types.h"
#include "time_step.h"

class pyTEC3DEquationSystems;

/** @brief Per-element context bundle populated by set_element(). */
struct ElemData
{
    pytec3d::MeshElement elem;                         ///< Local element.
    const ModelElementTypes::Entry *entry = nullptr; ///< Element type and domain of the element's subdomain.
    unsigned int material_num = 0;                     ///< Material ID (MATERIAL cell array).
    unsigned int real_constant_num = 0;                ///< Real-constant ID (CONSTANT cell array).
    unsigned int coordinate_system_num = 0;            ///< Coordinate-system ID (COORDINATE_SYSTEM cell array).
    pytec3d::registry::LawMask laws = 0;               ///< Laws selected by its TYPE, MATERIAL and CONSTANT rows.
};

/** @brief The local values (DofSpace::gather) of the vectors of the DOF system an element may read. */
struct LocalVectors
{
    const std::vector<double> *solution = nullptr;
    const std::vector<double> *old_solution = nullptr;
    const std::vector<double> *velocity = nullptr;
    const std::vector<double> *acceleration = nullptr;
    const std::vector<double> *old_velocity = nullptr;
    const std::vector<double> *old_acceleration = nullptr;
};

class ElementBuilder
{
public:
    using Matrix = pytec3d::la::DenseMatrix<pytec3d::la::Number>;
    using Vector = pytec3d::la::DenseVector<pytec3d::la::Number>;

    explicit ElementBuilder(const pyTEC3DEquationSystems &equation_systems);

    /** @brief Populates the context for local element @p e. */
    void set_element(int e);
    const ElemData &get_elem_data() const { return elem_data; }
    /** @brief True when laws are selected on the active element: it has a Jacobian and residual. */
    bool nonlinear_element() const { return elem_data.laws != 0; }
    unsigned int get_n_dofs() const { return n_dofs; }

    /** @brief The element values of the DOF system vectors the builder reads. */
    void set_vectors(const LocalVectors &vectors) { vectors_ = vectors; }

    /** @brief K, C and, for second-order element types, M of the element (a zero C or M empty, 0 x 0). */
    std::vector<Matrix> &dof_matrices() { return dof_matrices_; }
    /** @brief The operators of the static nodal / element fields: {F u, F du/dt, F d2u/dt2} (empty if zero). */
    std::vector<Matrix> &field_matrices(pytec3d::FieldKind kind)
    {
        return kind == pytec3d::FieldKind::Nodal ? nodal_matrices_ : element_matrices_;
    }

    /** @brief Builds K, C, (M) and the operators of the static nodal and element fields. */
    void build_element();

    /** @brief The source vector of the element (element loads, magnetization; the nodal loads are global). */
    Vector build_element_rhs();

    /** @brief The tangent matrix and residual of the law terms at the current solution; @p scheme gives the velocity
     * and acceleration of the law terms with dt() and the scaling of their Jacobian; @p picard their secant. */
    void build_element_tangent_stiffness_residual(const TimeScheme &scheme, bool picard, Matrix &tangent, Vector &residual);

    /** @brief True when the element has law outputs under dt() (ElementType::rate_residuals and a law selected). */
    bool rate_residual_element() const;

    /** @brief The outputs q(x) of the element's law terms with @p order dt() on the trial (1 or 2), the laws evaluated
     * at the local vector @p solution: the quantity whose rates the time scheme gives. Zero on other elements. */
    Vector law_output(const std::vector<double> &solution, unsigned int order);

    /** @brief True when the nodal / element fields of the active element have terms evaluated at every output: dynamic
     * fields (a law selected) or terms without trial (a load alone, e.g. Jt = Js). */
    bool has_evaluated_fields(pytec3d::FieldKind kind) const;
    /** @brief Adds those terms to @p values, the field values of the active element in the element vector of @p layout
     * (the static operators applied to the solution by the caller). */
    void add_evaluated_fields(pytec3d::FieldKind kind, const pytec3d::fem::SystemLayout &layout, Vector &values);

private:
    pytec3d::ElementParameterValues::Inputs parameter_inputs() const;
    /** @brief The values of the local vector @p local on the active element (zeros when null). */
    Vector element_vector(const std::vector<double> *local) const;
    pytec3d::ElementEvaluator::Inputs evaluator_inputs(pytec3d::Operation operation, pytec3d::ElementParameterValues &values,
                                                       const pytec3d::ElementType::FieldIndices &dofs) const;

    const pyTEC3DEquationSystems &equation_systems_;
    const pytec3d::fem::DofSpace &dofs_;
    const ModelElementTypes &element_types_;
    const pytec3d::ModelData &model_data_;
    const pytec3d::ModelArrays &arrays_;

    ElemData elem_data;
    unsigned int n_dofs = 0;
    LocalVectors vectors_;
    std::vector<Matrix> dof_matrices_;
    std::vector<Matrix> nodal_matrices_;
    std::vector<Matrix> element_matrices_;
};

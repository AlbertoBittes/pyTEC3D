/**
 * @file element_builder.cpp
 * @brief ElementBuilder: the matrices, vectors and fields of one element (element_builder.h).
 */
#include "element_builder.h"

#include <algorithm>
#include <string>

#include "logging.h"
#include "pytec3d_equation_systems.h"
#include "pytec3d_model.h"

using pytec3d::la::Number;
using pytec3d::la::Real;

namespace
{
    [[noreturn]] void fail(const std::string &message)
    {
        pytec3d::logging::error_and_abort("ElementBuilder: " + message);
    }
} // namespace

ElementBuilder::ElementBuilder(const pyTEC3DEquationSystems &equation_systems)
    : equation_systems_(equation_systems),
      dofs_(equation_systems.dof_space()),
      element_types_(equation_systems.element_types()),
      model_data_(equation_systems.get_model_data()),
      arrays_(equation_systems.model().arrays())
{
}

void ElementBuilder::set_element(int e)
{
    elem_data.elem = pytec3d::MeshElement{&equation_systems_.mesh(), &equation_systems_.catalog(),
                                          &equation_systems_.quadrature(), &equation_systems_.reference_tables(), e};
    const int subdomain = elem_data.elem.subdomain();
    elem_data.entry = &element_types_.get(subdomain);
    const auto cell = static_cast<std::size_t>(elem_data.elem.cell());
    elem_data.material_num = static_cast<unsigned int>(arrays_.material(cell));
    elem_data.real_constant_num = static_cast<unsigned int>(arrays_.constant(cell));
    elem_data.coordinate_system_num = static_cast<unsigned int>(arrays_.coordinate_system(cell));
    elem_data.laws = model_data_.laws(static_cast<unsigned int>(subdomain), elem_data.material_num, elem_data.real_constant_num) &
                     elem_data.entry->type->laws();
    n_dofs = dofs_.n_element_dofs(e);
}

pytec3d::ElementParameterValues::Inputs ElementBuilder::parameter_inputs() const
{
    pytec3d::ElementParameterValues::Inputs inputs;
    inputs.material = model_data_.material(elem_data.material_num);
    inputs.constants = model_data_.constants(elem_data.real_constant_num);
    inputs.material_id = elem_data.material_num;
    inputs.constants_id = elem_data.real_constant_num;
    inputs.domain = elem_data.entry->domain;
    inputs.arrays = &arrays_;
    inputs.element_id = static_cast<std::size_t>(elem_data.elem.cell());
    inputs.rotation = equation_systems_.element_rotation(elem_data.elem.index);
    return inputs;
}

pytec3d::ElementEvaluator::Inputs ElementBuilder::evaluator_inputs(pytec3d::Operation operation,
                                                                   pytec3d::ElementParameterValues &values,
                                                                   const pytec3d::ElementType::FieldIndices &dofs) const
{
    pytec3d::ElementEvaluator::Inputs inputs;
    inputs.type = elem_data.entry->type;
    inputs.domain = elem_data.entry->domain;
    inputs.operation = operation;
    inputs.laws = elem_data.laws;
    inputs.elem = &elem_data.elem;
    inputs.parameters = &values;
    inputs.model_data = &model_data_;
    inputs.material_id = elem_data.material_num;
    inputs.constants_id = elem_data.real_constant_num;
    inputs.type_id = static_cast<unsigned int>(elem_data.elem.subdomain());
    inputs.dofs = &dofs;
    return inputs;
}

ElementBuilder::Vector ElementBuilder::element_vector(const std::vector<double> *local) const
{
    Vector out(n_dofs);
    if (!local)
        return out;
    std::vector<double> values;
    dofs_.element_values(elem_data.elem.index, *local, values);
    for (unsigned int i = 0; i < n_dofs; ++i)
        out(i) = values[i];
    return out;
}

void ElementBuilder::build_element()
{
    const pytec3d::ElementType &type = *elem_data.entry->type;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(dofs_, elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
    if (dofs.n_dofs != n_dofs)
        fail("the fields cover " + std::to_string(dofs.n_dofs) + " DOFs but the element has " + std::to_string(n_dofs) + ".");

    // K and C; M only for second-order element types.
    Matrix stiffness(n_dofs, n_dofs);
    Matrix damping(n_dofs, n_dofs);
    Matrix mass;
    if (type.spec().time_order >= 2)
        mass.resize(n_dofs, n_dofs);
    pytec3d::ElementEvaluator evaluator(evaluator_inputs(pytec3d::Operation::Matrices, values, dofs));
    evaluator.add_matrices(stiffness, damping, mass);

    // A zero matrix is stored empty (0 x 0): it is skipped where it is used and the element matrix files stay small.
    const auto drop_if_zero = [](Matrix &m)
    {
        const auto &entries = m.get_values();
        if (std::all_of(entries.begin(), entries.end(), [](Number x) { return x == 0.0; }))
            m.resize(0, 0);
    };
    drop_if_zero(damping);
    drop_if_zero(mass);

    dof_matrices_ = {std::move(stiffness), std::move(damping)};
    if (type.spec().time_order >= 2)
        dof_matrices_.push_back(std::move(mass));

    // The operators of the static nodal and element fields: {F u, F du/dt, F d2u/dt2} per kind. Their terms without
    // trial (a load alone, e.g. Jt = Js) are not stored: build_fields evaluates them at every output.
    for (const pytec3d::FieldKind kind : {pytec3d::FieldKind::Nodal, pytec3d::FieldKind::Element})
    {
        const pytec3d::fem::SystemLayout &layout =
            kind == pytec3d::FieldKind::Nodal ? equation_systems_.nodal_layout() : equation_systems_.element_layout();
        const auto field_indices = type.element_indices(layout, elem_data.elem, kind, elem_data.entry->domain);
        std::array<Matrix, 3> F;
        evaluator.field_operators(kind, field_indices, F);
        std::vector<Matrix> &stored = field_matrices(kind);
        stored = {std::move(F[0]), std::move(F[1]), std::move(F[2])};
        for (Matrix &m : stored)
            drop_if_zero(m);
    }
}

ElementBuilder::Vector ElementBuilder::build_element_rhs()
{
    const pytec3d::ElementType &type = *elem_data.entry->type;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(dofs_, elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
    Vector rhs(n_dofs);
    pytec3d::ElementEvaluator evaluator(evaluator_inputs(pytec3d::Operation::Rhs, values, dofs));
    evaluator.add_rhs(rhs);
    return rhs;
}

void ElementBuilder::build_element_tangent_stiffness_residual(const TimeScheme &scheme, bool picard, Matrix &tangent,
                                                              Vector &residual)
{
    tangent.resize(n_dofs, n_dofs);
    residual.resize(n_dofs);
    if (!nonlinear_element())
        return;

    const pytec3d::ElementType &type = *elem_data.entry->type;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(dofs_, elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
    const Vector u = element_vector(vectors_.solution);

    auto inputs = evaluator_inputs(pytec3d::Operation::Nonlinear, values, dofs);
    inputs.solution = &u;
    inputs.jacobian_scale = {1.0, scheme.c_v, scheme.c_a};
    inputs.picard = picard;

    // The rates of the law terms with dt(): v = c_v u + h_v, a = c_a u + h_a (history of the step).
    Vector v, a;
    if (type.law_rates() && scheme.transient())
    {
        const Vector u_n = element_vector(vectors_.old_solution);
        const Vector v_n = element_vector(vectors_.old_velocity);
        const Vector a_n = element_vector(vectors_.old_acceleration);
        v.resize(n_dofs);
        a.resize(n_dofs);
        for (unsigned int i = 0; i < n_dofs; ++i)
        {
            v(i) = scheme.c_v * u(i) + scheme.hu_v * u_n(i) + scheme.hv_v * v_n(i) + scheme.ha_v * a_n(i);
            a(i) = scheme.c_a * u(i) + scheme.hu_a * u_n(i) + scheme.hv_a * v_n(i) + scheme.ha_a * a_n(i);
        }
        inputs.solution_dt = &v;
        inputs.solution_ddt = &a;
    }
    pytec3d::ElementEvaluator evaluator(inputs);
    evaluator.add_nonlinear(tangent, residual);

    // The rate of a law output under dt() (dt(B(H))): the time scheme applied to the output q(u) as to the DOFs,
    // c_d q(u) (add_nonlinear) + hu_d q(u_n) (here) + hv_d q'_n + ha_d q''_n (the history the system stores per
    // element and adds, law_output).
    if (type.rate_residuals() && scheme.transient())
    {
        const Vector u_n = element_vector(vectors_.old_solution);
        auto old_inputs = evaluator_inputs(pytec3d::Operation::Nonlinear, values, dofs);
        old_inputs.solution = &u_n;
        pytec3d::ElementEvaluator history(old_inputs);
        history.add_law_output(residual, {0.0, scheme.hu_v, scheme.hu_a});
    }
}

bool ElementBuilder::rate_residual_element() const
{
    return nonlinear_element() && elem_data.entry->type->rate_residuals();
}

ElementBuilder::Vector ElementBuilder::law_output(const std::vector<double> &solution, unsigned int order)
{
    Vector out(n_dofs);
    if (!rate_residual_element() || order < 1 || order > 2)
        return out;
    const pytec3d::ElementType &type = *elem_data.entry->type;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(dofs_, elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
    const Vector x = element_vector(&solution);
    auto inputs = evaluator_inputs(pytec3d::Operation::Nonlinear, values, dofs);
    inputs.solution = &x;
    std::array<Real, 3> scale = {0.0, 0.0, 0.0};
    scale[order] = 1.0;
    pytec3d::ElementEvaluator evaluator(inputs);
    evaluator.add_law_output(out, scale);
    return out;
}

bool ElementBuilder::has_evaluated_fields(pytec3d::FieldKind kind) const
{
    const pytec3d::ElementType &type = *elem_data.entry->type;
    const pytec3d::registry::Domain domain = elem_data.entry->domain;
    return type.has_dynamic_fields(kind, domain, elem_data.laws) || (type.field_operators(kind, domain) & 8u) != 0;
}

void ElementBuilder::add_evaluated_fields(pytec3d::FieldKind kind, const pytec3d::fem::SystemLayout &layout, Vector &values)
{
    if (!has_evaluated_fields(kind))
        return;
    const pytec3d::ElementType &type = *elem_data.entry->type;
    const pytec3d::registry::Domain domain = elem_data.entry->domain;
    pytec3d::ElementParameterValues parameters(type, parameter_inputs());
    const auto dofs = type.element_indices(dofs_, elem_data.elem, pytec3d::FieldKind::Dof, domain);
    const auto field_indices = type.element_indices(layout, elem_data.elem, kind, domain);
    // The solution and its rates: the velocity and acceleration of the time scheme, the ones of the equations at the end
    // of the step (with theta < 1 the mean rate (u - u_n) / dt would not satisfy them).
    const Vector u = element_vector(vectors_.solution);
    const Vector u_dt = element_vector(vectors_.velocity);
    const Vector u_ddt = element_vector(vectors_.acceleration);
    auto inputs = evaluator_inputs(kind == pytec3d::FieldKind::Nodal ? pytec3d::Operation::NodalFields
                                                                      : pytec3d::Operation::ElementFields,
                                   parameters, dofs);
    inputs.solution = &u;
    inputs.solution_dt = &u_dt;
    inputs.solution_ddt = &u_ddt;
    pytec3d::ElementEvaluator evaluator(inputs);
    evaluator.field_values(kind, field_indices, values);
}

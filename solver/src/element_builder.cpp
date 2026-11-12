/**
 * @file element_builder.cpp
 * @brief ElementBuilder: per-element inputs, element vector indices and evaluator calls.
 */
#include "element_builder.h"

#include <algorithm>
#include <array>
#include <mutex>
#include <tuple>
#include <string>

#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>

#include "logging.h"
#include "pytec3d_equation_systems.h"
#include "pytec3d_model.h"
#include "pytec3d_system.h"

using libMesh::DenseMatrix;
using libMesh::DenseVector;
using libMesh::dof_id_type;
using libMesh::Number;
using libMesh::NumericVector;
using libMesh::Real;
using libMesh::System;

namespace
{
    [[noreturn]] void fail(const std::string &message)
    {
        pytec3d::logging::error_and_abort("ElementBuilder: " + message);
    }
} // namespace

ElementBuilder::ElementBuilder(pyTEC3DEquationSystems *equation_systems)
    :       matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS)
{
    equation_systems_ = equation_systems;
    dof_solution_system_ = &equation_systems_->get_system<libMesh::ImplicitSystem>(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM));

    nodal_solution_system_ = &equation_systems_->get_system(pyTEC3DSystemBase::NODAL_SOLUTION);
    element_solution_system_ = &equation_systems_->get_system(pyTEC3DSystemBase::ELEMENT_SOLUTION);
    dof_solution_dof_map_ = &dof_solution_system_->get_dof_map();

    arrays_ = &equation_systems_->model().arrays();

    element_types_ = &equation_systems_->element_types();
    if (!element_types_)
        fail("element types factory is not available.");

    model_data_ = &equation_systems_->get_model_data();
    element_local_to_global_rotations_ = &equation_systems_->get_element_local_to_global_rotations();
    node_local_to_global_rotations_ = &equation_systems_->get_node_local_to_global_rotations();
}

bool ElementBuilder::should_inject_nodal_load(const std::string &dof_variable_name, dof_id_type node_id) const
{
    if (!nodal_load_map_ || !elem_data.elem)
        return false;

    const auto var_it = nodal_load_map_->find(dof_variable_name);
    if (var_it == nodal_load_map_->end())
        return false;

    const auto &injector_elem_by_node = var_it->second;
    const auto node_it = injector_elem_by_node.find(node_id);
    return (node_it != injector_elem_by_node.end()) &&
           (node_it->second == static_cast<dof_id_type>(elem_data.elem->id()));
}

void ElementBuilder::add_nodal_loads(const pytec3d::ElementType::FieldIndices &dofs, DenseVector<Number> &rhs) const
{
    // Nodal loads are point values on the nodes, not integrals; each node is loaded by one element only
    // (nodal_load_map_). A scalar load goes to the DOF of the node. A vector load of a vector variable is given in the
    // node frame (node_angles, like the Dirichlet values _F): f = R f_node in global, and each global component the
    // field has goes to its DOF; a component the field does not have does no work on it and is dropped. A field with
    // a single component (A_z of a 2D element, a scalar variable) ignores the node frames, as its Dirichlet values and
    // couplings do: its component of the load goes to its DOF as given. A NEDELEC_ONE field has its edge DOF on the
    // midside node: a scalar load there goes to it, in the orientation of the DOF (the circulation's).
    const pytec3d::ElementType &type = *elem_data.entry->type;
    const pytec3d::registry::Domain domain = elem_data.entry->domain;
    const unsigned int sys = dof_solution_system_->number();
    for (std::size_t index = 0; index < type.spec().fields.size(); ++index)
    {
        const pytec3d::FieldSpec &field = type.field(index);
        if (field.kind != pytec3d::FieldKind::Dof || field.nodal_load == pytec3d::registry::FieldId::COUNT ||
            !pytec3d::ElementType::active(field, domain))
            continue;
        const bool edge = field.family == libMesh::NEDELEC_ONE;
        if (!edge && field.family != libMesh::LAGRANGE && field.family != libMesh::LAGRANGE_VEC)
            fail("the nodal load of '" + std::string(pytec3d::registry::name(field.id)) +
                 "' needs a Lagrange or NEDELEC_ONE field.");
        const std::string load_array(pytec3d::registry::name(field.nodal_load));
        const pytec3d::RealArray *array = arrays_->point_field(load_array);
        const std::string label(pytec3d::registry::info(field.variable).label);
        if (!array || !dof_solution_system_->has_variable(label))
            continue;
        // An edge field: a scalar array whose value on the midside point of an edge (the node of its DOF) loads it.
        const bool vector = !edge && field.value == pytec3d::formulation::ValueType::Vector;
        const bool single_component = vector && field.n_components == 1;
        if (array->n_components() != (vector ? 3u : 1u))
            fail("the nodal load '" + load_array + "' has " + std::to_string(array->n_components()) + " components, expected " +
                 (vector ? "3" : "1") + ".");

        const unsigned int var_num = dof_solution_system_->variable_number(label);
        const auto &range = dofs.field[index];
        // Position of a DOF in the element vector: searched in the block of the field (DofMap::dof_indices order),
        // so that no layout of the components of a vector variable is assumed.
        const auto position = [&](dof_id_type dof) -> unsigned int
        {
            const auto first = dof_indices_.begin() + range.start;
            const auto it = std::find(first, first + range.count, dof);
            if (it == first + range.count)
                fail("a nodal load DOF of '" + label + "' is not among the DOFs of element " +
                     std::to_string(elem_data.elem->id()) + ".");
            return static_cast<unsigned int>(it - dof_indices_.begin());
        };

        for (unsigned int i = 0; i < elem_data.elem->n_nodes(); ++i)
        {
            const libMesh::Node &node = elem_data.elem->node_ref(i);
            if (node.n_dofs(sys, var_num) == 0 || !should_inject_nodal_load(label, node.id()))
                continue;
            const auto point = static_cast<std::size_t>(node.id());
            if (!vector)
            {
                rhs(position(node.dof_number(sys, var_num, 0))) += (*array)(point);
                continue;
            }
            if (single_component)
            {
                rhs(position(node.dof_number(sys, var_num, 0))) += (*array)(point, field.components[0]);
                continue;
            }

            const auto rotation = node_local_to_global_rotations_->find(node.id());
            if (rotation == node_local_to_global_rotations_->end())
                fail("no node rotation for node " + std::to_string(node.id()) + ".");
            const libMesh::RealVectorValue global =
                rotation->second * libMesh::RealVectorValue((*array)(point, 0), (*array)(point, 1), (*array)(point, 2));
            // The values of the field on a node: its n_components physical components, DOF component j of the
            // variable (a single component is a scalar variable: DOF component 0).
            for (unsigned int j = 0; j < field.n_components; ++j)
                rhs(position(node.dof_number(sys, var_num, j))) += global(field.components[j]);
        }
    }
}

void ElementBuilder::set_element_data(const libMesh::Elem *elem)
{
    elem_data.entry = &element_types_->get(elem->subdomain_id());

    const auto elem_id = static_cast<std::size_t>(elem->id());
    elem_data.material_num = static_cast<unsigned int>(arrays_->material(elem_id));
    elem_data.real_constant_num = static_cast<unsigned int>(arrays_->constant(elem_id));
    elem_data.coordinate_system_num = static_cast<unsigned int>(arrays_->coordinate_system(elem_id));
    elem_data.elem = elem;
    elem_data.laws = model_data_->laws(elem->subdomain_id(), elem_data.material_num, elem_data.real_constant_num) &
                     elem_data.entry->type->laws();

    dof_solution_dof_map_->dof_indices(elem, dof_indices_);
    n_dofs = static_cast<unsigned int>(dof_indices_.size());
}

const ElemData &ElementBuilder::get_elem_data() const
{
    return elem_data;
}

bool ElementBuilder::nonlinear_element() const
{
    return elem_data.laws != 0;
}

pytec3d::ElementParameterValues::Inputs ElementBuilder::parameter_inputs() const
{
    pytec3d::ElementParameterValues::Inputs inputs;
    inputs.material = model_data_->material(elem_data.material_num);
    inputs.constants = model_data_->constants(elem_data.real_constant_num);
    inputs.material_id = elem_data.material_num;
    inputs.constants_id = elem_data.real_constant_num;
    inputs.domain = elem_data.entry->domain;
    inputs.arrays = arrays_;
    inputs.element_id = static_cast<std::size_t>(elem_data.elem->id());

    if (element_local_to_global_rotations_)
    {
        const auto rotation = element_local_to_global_rotations_->find(static_cast<dof_id_type>(elem_data.elem->id()));
        if (rotation != element_local_to_global_rotations_->end())
            inputs.rotation = rotation->second;
    }
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
    inputs.elem = elem_data.elem;
    inputs.parameters = &values;
    inputs.model_data = model_data_;
    inputs.material_id = elem_data.material_num;
    inputs.constants_id = elem_data.real_constant_num;
    inputs.type_id = static_cast<unsigned int>(elem_data.elem->subdomain_id());
    inputs.dofs = &dofs;
    return inputs;
}

DenseVector<Number> ElementBuilder::local_dof_vector(const NumericVector<Number> *solution) const
{
    DenseVector<Number> local(static_cast<unsigned int>(dof_indices_.size()));
    if (solution)
        for (unsigned int i = 0; i < dof_indices_.size(); ++i)
            local(i) = solution->el(dof_indices_[i]);
    return local;
}

void ElementBuilder::build_element()
{
    const pytec3d::ElementType &type = *elem_data.entry->type;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
    if (dofs.n_dofs != n_dofs)
        fail("the fields cover " + std::to_string(dofs.n_dofs) + " DOFs but the element has " + std::to_string(n_dofs) + ".");

    // K and C; M only for second-order element types.
    DenseMatrix<Number> stiffness(n_dofs, n_dofs);
    DenseMatrix<Number> damping(n_dofs, n_dofs);
    DenseMatrix<Number> mass;
    if (type.spec().time_order >= 2)
        mass.resize(n_dofs, n_dofs);
    pytec3d::ElementEvaluator evaluator(evaluator_inputs(pytec3d::Operation::Matrices, values, dofs));
    evaluator.add_matrices(stiffness, damping, mass);

    // A zero matrix is stored empty (0 x 0): it is skipped where it is used and the element matrix files stay small.
    const auto drop_if_zero = [](DenseMatrix<Number> &m)
    {
        const auto &entries = m.get_values();
        if (std::all_of(entries.begin(), entries.end(), [](Number x) { return x == 0.0; }))
            m.resize(0, 0);
    };
    drop_if_zero(damping);
    drop_if_zero(mass);

    matrices[MatricesSet::DOF_SOLUTION_MATRICES] = {std::move(stiffness), std::move(damping)};
    if (type.spec().time_order >= 2)
        matrices[MatricesSet::DOF_SOLUTION_MATRICES].push_back(std::move(mass));

    // The operators of the static nodal and element fields: {F u, F du/dt, F d2u/dt2} per kind. Their terms without
    // trial (a load alone, e.g. Jt = Js) are not stored: build_fields evaluates them at every output.
    for (const auto &[kind, set, system] : {std::make_tuple(pytec3d::FieldKind::Nodal, MatricesSet::NODAL_SOLUTION_MATRICES, nodal_solution_system_),
                                            std::make_tuple(pytec3d::FieldKind::Element, MatricesSet::ELEMENT_SOLUTION_MATRICES, element_solution_system_)})
    {
        const auto field_indices = type.element_indices(*system, *elem_data.elem, kind, elem_data.entry->domain);
        std::array<DenseMatrix<Number>, 3> F;
        evaluator.field_operators(kind, field_indices, F);

        std::vector<DenseMatrix<Number>> &stored = matrices[set];
        stored = {std::move(F[0]), std::move(F[1]), std::move(F[2])};
        for (DenseMatrix<Number> &m : stored)
            drop_if_zero(m);
    }
}

void ElementBuilder::build_element_rhs()
{
    const pytec3d::ElementType &type = *elem_data.entry->type;
    const pytec3d::registry::Domain domain = elem_data.entry->domain;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, domain);

    DenseVector<Number> rhs(n_dofs);
    pytec3d::ElementEvaluator evaluator(evaluator_inputs(pytec3d::Operation::Rhs, values, dofs));
    evaluator.add_rhs(rhs);

    add_nodal_loads(dofs, rhs);

    vectors[VectorsSet::RHS_VECTORS] = {std::move(rhs)};
}

void ElementBuilder::build_element_tangent_stiffness_residual(const TimeScheme &scheme, bool picard)
{
    DenseMatrix<Number> tangent(n_dofs, n_dofs);
    DenseVector<Number> residual(n_dofs);

    if (nonlinear_element())
    {
        const pytec3d::ElementType &type = *elem_data.entry->type;
        pytec3d::ElementParameterValues values(type, parameter_inputs());
        const auto dofs = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
        const DenseVector<Number> u = local_dof_vector(dof_solution_system_->current_local_solution.get());

        auto inputs = evaluator_inputs(pytec3d::Operation::Nonlinear, values, dofs);
        inputs.solution = &u;
        inputs.jacobian_scale = {1.0, scheme.c_v, scheme.c_a};
        inputs.picard = picard;

        // The rates of the law terms with dt(): v = c_v u + h_v, a = c_a u + h_a (history of the step).
        DenseVector<Number> v, a;
        if (type.law_rates() && scheme.transient())
        {
            const DenseVector<Number> u_n = local_dof_vector(&dof_solution_system_->get_vector("_transient_old_local_solution"));
            const DenseVector<Number> v_n = local_dof_vector(&dof_solution_system_->get_vector(pyTEC3DSystemBase::kOldVelocity));
            const DenseVector<Number> a_n = local_dof_vector(&dof_solution_system_->get_vector(pyTEC3DSystemBase::kOldAcceleration));
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
        // element and adds, ElementBuilder::law_output).
        if (type.rate_residuals() && scheme.transient())
        {
            const DenseVector<Number> u_n = local_dof_vector(&dof_solution_system_->get_vector("_transient_old_local_solution"));
            auto old_inputs = evaluator_inputs(pytec3d::Operation::Nonlinear, values, dofs);
            old_inputs.solution = &u_n;
            pytec3d::ElementEvaluator history(old_inputs);
            history.add_law_output(residual, {0.0, scheme.hu_v, scheme.hu_a});
        }
    }

    matrices[MatricesSet::TANGENT_STIFFNESS_MATRICES] = {std::move(tangent)};
    vectors[VectorsSet::RESIDUAL_VECTORS] = {std::move(residual)};
}

bool ElementBuilder::rate_residual_element() const
{
    return nonlinear_element() && elem_data.entry->type->rate_residuals();
}

DenseVector<Number> ElementBuilder::law_output(const libMesh::NumericVector<Number> &solution, unsigned int order)
{
    DenseVector<Number> out(n_dofs);
    if (!rate_residual_element() || order < 1 || order > 2)
        return out;
    const pytec3d::ElementType &type = *elem_data.entry->type;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, elem_data.entry->domain);
    const DenseVector<Number> x = local_dof_vector(&solution);
    auto inputs = evaluator_inputs(pytec3d::Operation::Nonlinear, values, dofs);
    inputs.solution = &x;
    std::array<libMesh::Real, 3> scale = {0.0, 0.0, 0.0};
    scale[order] = 1.0;
    pytec3d::ElementEvaluator evaluator(inputs);
    evaluator.add_law_output(out, scale);
    return out;
}

void ElementBuilder::build_nodal_solution()
{
    build_fields(pytec3d::FieldKind::Nodal, *nodal_solution_system_);
}

void ElementBuilder::build_element_solution()
{
    build_fields(pytec3d::FieldKind::Element, *element_solution_system_);
}

void ElementBuilder::build_fields(pytec3d::FieldKind kind, System &system)
{
    static std::mutex solution_write_mutex;

    const pytec3d::ElementType &type = *elem_data.entry->type;
    const pytec3d::registry::Domain domain = elem_data.entry->domain;
    pytec3d::ElementParameterValues values(type, parameter_inputs());
    const auto dofs = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, domain);
    const auto field_indices = type.element_indices(system, *elem_data.elem, kind, domain);

    // The static fields: their operators, stored with the element matrices (build_element), applied to the solution.
    const auto &stored = matrices[kind == pytec3d::FieldKind::Nodal ? MatricesSet::NODAL_SOLUTION_MATRICES
                                                                    : MatricesSet::ELEMENT_SOLUTION_MATRICES];
    const bool dynamic = type.has_dynamic_fields(kind, domain, elem_data.laws);
    const bool constants = (type.field_operators(kind, domain) & 8u) != 0;
    const auto needs = [&](std::size_t k) { return dynamic || (k < stored.size() && stored[k].m() > 0); };

    // The solution and, where a field reads them, its rates: the velocity and acceleration of the time scheme, the
    // ones of the equations at the end of the step (with theta < 1 the mean rate (u - u_n) / dt would not satisfy them).
    const DenseVector<Number> u = local_dof_vector(dof_solution_system_->current_local_solution.get());
    const DenseVector<Number> u_dt = needs(1) ? local_dof_vector(&dof_solution_system_->get_vector(pyTEC3DSystemBase::kVelocity))
                                              : DenseVector<Number>();
    const DenseVector<Number> u_ddt = needs(2) ? local_dof_vector(&dof_solution_system_->get_vector(pyTEC3DSystemBase::kAcceleration))
                                               : DenseVector<Number>();
    DenseVector<Number> field_values(field_indices.n_dofs);
    const DenseVector<Number> *solutions[3] = {&u, &u_dt, &u_ddt};
    for (std::size_t k = 0; k < 3 && k < stored.size(); ++k)
        if (stored[k].m() > 0)
            stored[k].vector_mult_add(field_values, 1.0, *solutions[k]);

    // Evaluated now: the dynamic fields (a law selected on the element), which depend on the solution, and the terms
    // without trial of the static ones, which read the current loads.
    if (dynamic || constants)
    {
        auto inputs = evaluator_inputs(kind == pytec3d::FieldKind::Nodal ? pytec3d::Operation::NodalFields : pytec3d::Operation::ElementFields,
                                       values, dofs);
        inputs.solution = &u;
        inputs.solution_dt = &u_dt;
        inputs.solution_ddt = &u_ddt;
        pytec3d::ElementEvaluator evaluator(inputs);
        evaluator.field_values(kind, field_indices, field_values);
    }

    std::vector<dof_id_type> system_dofs;
    system.get_dof_map().dof_indices(elem_data.elem, system_dofs);
    if (system_dofs.size() != field_values.size())
        fail("the fields have " + std::to_string(field_values.size()) + " values but system '" + system.name() + "' has " +
             std::to_string(system_dofs.size()) + " DOFs on the element.");

    std::lock_guard<std::mutex> lock(solution_write_mutex);
    for (unsigned int i = 0; i < system_dofs.size(); ++i)
        if ((system.solution->first_local_index() <= system_dofs[i]) && (system_dofs[i] < system.solution->last_local_index()))
            system.solution->set(system_dofs[i], field_values(i));
}

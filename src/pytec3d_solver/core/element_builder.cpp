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
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

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
    : FEParentObjectBase<pyTEC3DEquationSystems>(equation_systems),
      matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS)
{
    equation_systems_ = equation_systems;
    dof_solution_system_ = &equation_systems_->get_system<libMesh::ImplicitSystem>(
        pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM));

    nodal_solution_system_ = &equation_systems_->get_system(pyTEC3DSystemBase::NODAL_SOLUTION);
    element_solution_system_ = &equation_systems_->get_system(pyTEC3DSystemBase::ELEMENT_SOLUTION);
    dof_solution_dof_map_ = &dof_solution_system_->get_dof_map();

    vtk_grid_ = equation_systems_->get_parent_object()->get_active_vtk_grid();
    if (!vtk_grid_)
        fail("active VTK grid is not available.");
    material_cell_array_ = vtk_grid_->GetCellData()->GetArray("MATERIAL");
    if (!material_cell_array_)
        fail("MATERIAL array is not available in the active VTK grid.");
    real_constant_cell_array_ = vtk_grid_->GetCellData()->GetArray("CONSTANT");
    if (!real_constant_cell_array_)
        fail("CONSTANT array is not available in the active VTK grid.");
    coordinate_system_cell_array_ = vtk_grid_->GetCellData()->GetArray("COORDINATE_SYSTEM");
    if (!coordinate_system_cell_array_)
        fail("COORDINATE_SYSTEM array is not available in the active VTK grid.");

    element_types_factory_ = equation_systems_->get_element_types_factory().get();
    if (!element_types_factory_)
        fail("element types factory is not available.");

    model_data_ = &equation_systems_->get_model_data();
    element_local_to_global_rotations_ = &equation_systems_->get_element_local_to_global_rotations();
}

bool ElementBuilder::should_inject_nodal_load_dof(const std::string &dof_variable_name,
                                                  dof_id_type dof_id) const
{
    if (!nodal_load_map_ || !elem_data.elem)
        return false;

    const auto var_it = nodal_load_map_->find(dof_variable_name);
    if (var_it == nodal_load_map_->end())
        return false;

    const auto &injector_elem_by_dof = var_it->second;
    const auto dof_it = injector_elem_by_dof.find(dof_id);
    return (dof_it != injector_elem_by_dof.end()) &&
           (dof_it->second == static_cast<dof_id_type>(elem_data.elem->id()));
}

void ElementBuilder::set_element_data(const libMesh::Elem *elem)
{
    elem_data.entry = &element_types_factory_->get(elem->subdomain_id());

    const vtkIdType elem_id = static_cast<vtkIdType>(elem->id());
    elem_data.material_num = static_cast<unsigned int>(material_cell_array_->GetComponent(elem_id, 0));
    elem_data.real_constant_num = static_cast<unsigned int>(real_constant_cell_array_->GetComponent(elem_id, 0));
    elem_data.coordinate_system_num = static_cast<unsigned int>(coordinate_system_cell_array_->GetComponent(elem_id, 0));
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
    inputs.cell_data = vtk_grid_->GetCellData();
    inputs.element_id = static_cast<vtkIdType>(elem_data.elem->id());

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

    // The operators of the static nodal and element fields: {F u, F du/dt, F d2u/dt2, constant} per kind.
    for (const auto &[kind, set, system] : {std::make_tuple(pytec3d::FieldKind::Nodal, MatricesSet::NODAL_SOLUTION_MATRICES, nodal_solution_system_),
                                            std::make_tuple(pytec3d::FieldKind::Element, MatricesSet::ELEMENT_SOLUTION_MATRICES, element_solution_system_)})
    {
        const auto field_indices = type.element_indices(*system, *elem_data.elem, kind, elem_data.entry->domain);
        std::array<DenseMatrix<Number>, 3> F;
        DenseVector<Number> f;
        evaluator.field_operators(kind, field_indices, F, f);

        DenseMatrix<Number> constant(f.size(), f.size() ? 1 : 0);
        for (unsigned int i = 0; i < f.size(); ++i)
            constant(i, 0) = f(i);
        std::vector<DenseMatrix<Number>> &stored = matrices[set];
        stored = {std::move(F[0]), std::move(F[1]), std::move(F[2]), std::move(constant)};
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

    // Nodal loads are point values on the nodes, not integrals.
    for (std::size_t index = 0; index < type.spec().fields.size(); ++index)
    {
        const pytec3d::FieldSpec &field = type.field(index);
        if (field.kind != pytec3d::FieldKind::Dof || field.nodal_load == pytec3d::registry::FieldId::COUNT ||
            !pytec3d::ElementType::active(field, domain))
            continue;
        const std::string load_array(pytec3d::registry::name(field.nodal_load));
        vtkDataArray *array = vtk_grid_->GetPointData()->GetArray(load_array.c_str());
        const std::string label(pytec3d::registry::info(field.variable).label);
        if (!array || !dof_solution_system_->has_variable(label))
            continue;

        const unsigned int var_num = dof_solution_system_->variable_number(label);
        const auto &ids = dofs.field[index];
        for (unsigned int i = 0; i < ids.size(); ++i)
        {
            const libMesh::Node &node = elem_data.elem->node_ref(i);
            const dof_id_type dof_id = node.dof_number(dof_solution_system_->number(), var_num, 0);
            if (should_inject_nodal_load_dof(label, dof_id))
                rhs(ids[i]) += array->GetComponent(static_cast<vtkIdType>(node.id()), 0);
        }
    }

    vectors[VectorsSet::RHS_VECTORS] = {std::move(rhs)};
}

void ElementBuilder::build_element_tangent_stiffness_residual(const TimeScheme &scheme)
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
    }

    matrices[MatricesSet::TANGENT_STIFFNESS_MATRICES] = {std::move(tangent)};
    vectors[VectorsSet::RESIDUAL_VECTORS] = {std::move(residual)};
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
    const auto needs = [&](std::size_t k) { return dynamic || (k < stored.size() && stored[k].m() > 0); };

    // The solution and, where a field reads them, its rates: first-order models report the mean rate over the
    // step, second-order ones the scheme velocity.
    const DenseVector<Number> u = local_dof_vector(dof_solution_system_->current_local_solution.get());
    const char *velocity = equation_systems_->time_order() >= 2 ? pyTEC3DSystemBase::kVelocity : pyTEC3DSystemBase::kMeanVelocity;
    const DenseVector<Number> u_dt = needs(1) ? local_dof_vector(&dof_solution_system_->get_vector(velocity)) : DenseVector<Number>();
    const DenseVector<Number> u_ddt = needs(2) ? local_dof_vector(&dof_solution_system_->get_vector(pyTEC3DSystemBase::kAcceleration))
                                               : DenseVector<Number>();
    DenseVector<Number> field_values(field_indices.n_dofs);
    const DenseVector<Number> *solutions[3] = {&u, &u_dt, &u_ddt};
    for (std::size_t k = 0; k < 3 && k < stored.size(); ++k)
        if (stored[k].m() > 0)
            stored[k].vector_mult_add(field_values, 1.0, *solutions[k]);
    if (stored.size() > 3 && stored[3].m() > 0)
        for (unsigned int i = 0; i < field_values.size(); ++i)
            field_values(i) += stored[3](i, 0);

    // The dynamic fields (a law selected on the element) depend on the solution: evaluated now.
    if (dynamic)
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

/**
 * @file element_builder.cpp
 * @brief ElementBuilder: per-element reinit, matrix indices and kernel calls.
 */
#include "element_builder.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>

#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include "bh_curve.h"
#include "kernels.h"
#include "logging.h"
#include "nabla.h"
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
        pytec3d::logging::error_and_abort("ElementBuilder: active VTK grid is not available.");
    material_cell_array_ = vtk_grid_->GetCellData()->GetArray("MATERIAL");
    if (!material_cell_array_)
        pytec3d::logging::error_and_abort("ElementBuilder: MATERIAL array is not available in the active VTK grid.");
    real_constant_cell_array_ = vtk_grid_->GetCellData()->GetArray("CONSTANT");
    if (!real_constant_cell_array_)
        pytec3d::logging::error_and_abort("ElementBuilder: CONSTANT array is not available in the active VTK grid.");
    coordinate_system_cell_array_ = vtk_grid_->GetCellData()->GetArray("COORDINATE_SYSTEM");
    if (!coordinate_system_cell_array_)
        pytec3d::logging::error_and_abort("ElementBuilder: COORDINATE_SYSTEM array is not available in the active VTK grid.");

    element_types_factory_ = equation_systems_->get_element_types_factory().get();
    if (!element_types_factory_)
        pytec3d::logging::error_and_abort("ElementBuilder: element types factory is not available.");

    model_data_ = &equation_systems_->get_model_data();
    bh_curves_ = &model_data_->bh_curves();
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

Real ElementBuilder::thickness() const
{
    const pytec3d::ConstantRow *constants = model_data_->constants(elem_data.real_constant_num);
    const pytec3d::ParameterValue *value = constants ? constants->find(pytec3d::registry::RealConstant::thickness) : nullptr;
    return value ? value->scalar : 1.0;
}

void ElementBuilder::set_element_data(const libMesh::Elem *elem)
{
    elem_data.entry = &element_types_factory_->get(elem->subdomain_id());

    const vtkIdType elem_id = static_cast<vtkIdType>(elem->id());
    elem_data.material_num = static_cast<unsigned int>(material_cell_array_->GetComponent(elem_id, 0));
    elem_data.real_constant_num = static_cast<unsigned int>(real_constant_cell_array_->GetComponent(elem_id, 0));
    elem_data.coordinate_system_num = static_cast<unsigned int>(coordinate_system_cell_array_->GetComponent(elem_id, 0));
    elem_data.elem = elem;

    std::vector<dof_id_type> dof_indices;
    dof_solution_dof_map_->dof_indices(elem, dof_indices);
    n_dofs = static_cast<unsigned int>(dof_indices.size());
}

const ElemData &ElementBuilder::get_elem_data() const
{
    return elem_data;
}

bool ElementBuilder::nonlinear_element() const
{
    return bh_curves_ && bh_curves_->find(elem_data.material_num) != bh_curves_->end();
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

    if (bh_curves_)
    {
        const auto curve = bh_curves_->find(elem_data.material_num);
        if (curve != bh_curves_->end())
            inputs.bh_curve = &curve->second;
    }
    return inputs;
}

namespace
{
    using pytec3d::ElementParameterValues;
    using pytec3d::ParameterType;

    [[noreturn]] void fail(const std::string &message)
    {
        pytec3d::logging::error_and_abort("ElementBuilder: " + message);
    }

    // Calls fn with the scalar or tensor value of a parameter used between a test and a trial.
    template <typename Fn>
    void with_matrix_parameter(ElementParameterValues &values, const pytec3d::ParameterRef &reference, Fn &&fn)
    {
        const pytec3d::ParameterValue value = values.value(reference);
        if (value.type == ParameterType::Scalar)
            fn(value.scalar);
        else if (value.type == ParameterType::Tensor)
            fn(value.tensor);
        else
            fail("parameter '" + std::string(reference.name()) + "' is a " + pytec3d::type_name(value.type) + " and cannot sit between a test and a trial.");
    }

    // Calls fn with the scalar or vector value of a source parameter.
    template <typename Fn>
    void with_source_parameter(ElementParameterValues &values, const pytec3d::ParameterRef &reference, Fn &&fn)
    {
        const pytec3d::ParameterValue value = values.value(reference);
        if (value.type == ParameterType::Scalar)
            fn(value.scalar);
        else if (value.type == ParameterType::Vector)
            fn(value.vector);
        else
            fail("parameter '" + std::string(reference.name()) + "' is a " + pytec3d::type_name(value.type) + " and cannot be a source.");
    }
} // namespace

/** @brief One element evaluation: the reinit()ed FEs, the integration weights and the operators evaluated on them. */
class ElementBuilder::Context
{
public:
    Context(const ElementTypesFactory::Entry &entry, const libMesh::Elem &elem, Real thickness)
        : type_(*entry.type)
    {
        type_.reinit(elem, entry.domain);

        weights_ = type_.jxw();
        length_.assign(weights_.size(), 1.0);
        if (type_.dimension() == 2)
        {
            const auto &xyz = type_.xyz();
            for (std::size_t qp = 0; qp < weights_.size(); ++qp)
            {
                // Out-of-plane length: the thickness, or the circumference of the axisymmetric element.
                length_[qp] = type_.options().axisymmetric ? thickness * 2.0 * libMesh::pi * std::abs(xyz[qp](0)) : thickness;
                weights_[qp] *= length_[qp];
            }
        }
    }

    const std::vector<Real> &weights() const { return weights_; }

    const pytec3d::Nabla::Operand &operand(const pytec3d::OperandRef &reference) { return evaluate(*reference.nabla); }

    /** @brief The shape functions of @p field. */
    const pytec3d::Nabla::Operand &basis(const pytec3d::Field &field)
    {
        return evaluate(type_.nabla(field, pytec3d::NablaOperator::Identity));
    }

private:
    const pytec3d::Nabla::Operand &evaluate(const pytec3d::Nabla &nabla)
    {
        auto it = cache_.find(&nabla);
        if (it == cache_.end())
        {
            const pytec3d::Field &field = nabla.field();
            pytec3d::Nabla::Operand value;
            if (nabla.op() == pytec3d::NablaOperator::Grad && type_.dimension() == 2 && !field.physical_vector)
            {
                // A lumped 2D potential varies along z only: its gradient is -phi / L, out of plane.
                value = type_.nabla(field, pytec3d::NablaOperator::Identity).evaluate();
                for (std::size_t qp = 0; qp < value.values.size(); ++qp)
                    value.values[qp].scale(-1.0 / length_[qp]);
                value.components = {2};
            }
            else
                value = nabla.evaluate();
            it = cache_.emplace(&nabla, std::move(value)).first;
        }
        return it->second;
    }

    const pytec3d::ElementType &type_;
    std::vector<Real> weights_;
    std::vector<Real> length_;
    std::map<const pytec3d::Nabla *, pytec3d::Nabla::Operand> cache_;
};

DenseVector<Number> ElementBuilder::local_dof_vector(const NumericVector<Number> *solution) const
{
    std::vector<dof_id_type> dof_indices;
    dof_solution_dof_map_->dof_indices(elem_data.elem, dof_indices);

    DenseVector<Number> local(static_cast<unsigned int>(dof_indices.size()));
    if (solution)
        for (unsigned int i = 0; i < dof_indices.size(); ++i)
            local(i) = solution->el(dof_indices[i]);
    return local;
}

void ElementBuilder::build_element()
{
    const ElementTypesFactory::Entry &entry = *elem_data.entry;
    const pytec3d::ElementType &type = *entry.type;
    const auto inputs = parameter_inputs();
    ElementParameterValues values(type, inputs);

    Context context(entry, *elem_data.elem, thickness());
    const auto indices = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, entry.domain);
    if (indices.n_dofs != n_dofs)
        fail("the fields cover " + std::to_string(indices.n_dofs) + " DOFs but the element has " + std::to_string(n_dofs) + ".");

    DenseMatrix<Number> stiffness(n_dofs, n_dofs);
    DenseMatrix<Number> time_derivative(n_dofs, n_dofs);

    for (const pytec3d::DofKernel &kernel : type.dof_kernels())
    {
        // Nonlinear laws enter through the tangent and residual.
        if (!kernel.has_trial || !type.active(kernel, entry.domain) || values.nonlinear(kernel.parameter))
            continue;

        const pytec3d::Nabla::Operand &test = context.operand(kernel.test);
        const pytec3d::Nabla::Operand &trial = context.operand(kernel.trial);
        DenseMatrix<Number> &target = kernel.trial.time_derivative ? time_derivative : stiffness;
        const auto &test_ids = indices.field[kernel.test.field().index];
        const auto &trial_ids = indices.field[kernel.trial.field().index];

        with_matrix_parameter(values, kernel.parameter, [&](const auto &value)
                              { pytec3d::kernel::Block(value, trial, test, context.weights(), trial_ids, test_ids, target); });
    }

    matrices[MatricesSet::DOF_SOLUTION_MATRICES] = {std::move(stiffness), std::move(time_derivative)};
}

void ElementBuilder::build_element_rhs()
{
    const ElementTypesFactory::Entry &entry = *elem_data.entry;
    const pytec3d::ElementType &type = *entry.type;
    const auto inputs = parameter_inputs();
    ElementParameterValues values(type, inputs);

    Context context(entry, *elem_data.elem, thickness());
    const auto indices = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, entry.domain);
    DenseVector<Number> rhs(n_dofs);

    for (const pytec3d::DofKernel &kernel : type.dof_kernels())
    {
        if (kernel.has_trial || !type.active(kernel, entry.domain))
            continue;

        const pytec3d::Nabla::Operand &test = context.operand(kernel.test);
        const auto &test_ids = indices.field[kernel.test.field().index];
        with_source_parameter(values, kernel.parameter, [&](const auto &value)
                              { pytec3d::kernel::Source(value, test, context.weights(), test_ids, rhs); });
    }

    // Nodal loads are point values on the nodes, not integrals.
    for (const pytec3d::Field &field : type.fields(pytec3d::FieldKind::Dof))
    {
        if (!field.nodal_loads || !field.active_in(entry.domain))
            continue;
        const std::string load_array(pytec3d::registry::name(*field.nodal_loads));
        vtkDataArray *array = vtk_grid_->GetPointData()->GetArray(load_array.c_str());
        const std::string &label = field.libmesh_variables.front().label;
        if (!array || !dof_solution_system_->has_variable(label))
            continue;

        const unsigned int var_num = dof_solution_system_->variable_number(label);
        const auto &ids = indices.field[field.index];
        for (unsigned int i = 0; i < ids.size(); ++i)
        {
            const libMesh::Node &node = elem_data.elem->node_ref(i);
            const dof_id_type dof_id = node.dof_number(dof_solution_system_->number(), var_num, 0);
            if (should_inject_nodal_load_dof(std::string(field.name()), dof_id))
                rhs(ids[i]) += array->GetComponent(static_cast<vtkIdType>(node.id()), 0);
        }
    }

    vectors[VectorsSet::RHS_VECTORS] = {std::move(rhs)};
}

void ElementBuilder::build_element_tangent_stiffness_residual()
{
    const ElementTypesFactory::Entry &entry = *elem_data.entry;
    const pytec3d::ElementType &type = *entry.type;
    const auto inputs = parameter_inputs();
    ElementParameterValues values(type, inputs);

    DenseMatrix<Number> tangent(n_dofs, n_dofs);
    DenseVector<Number> residual(n_dofs);

    if (inputs.bh_curve)
    {
        Context context(entry, *elem_data.elem, thickness());
        const auto indices = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, entry.domain);
        const DenseVector<Number> u = local_dof_vector(dof_solution_system_->current_local_solution.get());

        for (const pytec3d::DofKernel &kernel : type.dof_kernels())
        {
            if (!kernel.has_trial || !type.active(kernel, entry.domain) || !values.nonlinear(kernel.parameter))
                continue;

            const auto law = [&](const Real *B, Real *H, Real(*J)[3])
            { values.apply_law(kernel.parameter, B, H, J); };
            pytec3d::kernel::NonlinearBlock(law,
                                            context.operand(kernel.trial),
                                            context.operand(kernel.test),
                                            context.weights(),
                                            indices.field[kernel.trial.field().index],
                                            indices.field[kernel.test.field().index],
                                            u,
                                            tangent,
                                            residual);
        }
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

    const ElementTypesFactory::Entry &entry = *elem_data.entry;
    const pytec3d::ElementType &type = *entry.type;
    const auto inputs = parameter_inputs();
    ElementParameterValues values(type, inputs);

    Context context(entry, *elem_data.elem, thickness());
    const auto dof_indices = type.element_indices(*dof_solution_system_, *elem_data.elem, pytec3d::FieldKind::Dof, entry.domain);
    const auto field_indices = type.element_indices(system, *elem_data.elem, kind, entry.domain);

    const DenseVector<Number> u = local_dof_vector(dof_solution_system_->current_local_solution.get());
    const DenseVector<Number> u_dt = local_dof_vector(dof_solution_system_->request_vector("solution_time_derivative"));

    // field = F * u + F_dt * u_dt + (what the nonlinear laws give), F and F_dt acting on the local DOF solution and its time derivative.
    DenseMatrix<Number> F(field_indices.n_dofs, dof_indices.n_dofs);
    DenseMatrix<Number> F_dt(field_indices.n_dofs, dof_indices.n_dofs);
    DenseVector<Number> field_values(field_indices.n_dofs);

    const auto &kernels = kind == pytec3d::FieldKind::Nodal ? type.nodal_kernels() : type.element_kernels();
    for (const pytec3d::FieldKernel &kernel : kernels)
    {
        if (!type.active(kernel, entry.domain))
            continue;

        const pytec3d::Nabla::Operand &trial = context.operand(kernel.trial);
        DenseMatrix<Number> &target = kernel.trial.time_derivative ? F_dt : F;
        const auto &trial_ids = dof_indices.field[kernel.trial.field().index];
        const auto &field_ids = field_indices.field[kernel.field->index];
        const pytec3d::Nabla::Operand *basis = kind == pytec3d::FieldKind::Nodal ? &context.basis(*kernel.field) : nullptr;

        if (values.nonlinear(kernel.parameter))
        {
            const auto law = [&](const Real *B, Real *H, Real(*J)[3])
            { values.apply_law(kernel.parameter, B, H, J); };
            if (basis)
                pytec3d::kernel::NonlinearNodalField(law, trial, *basis, context.weights(), trial_ids, field_ids, u, field_values);
            else
                pytec3d::kernel::NonlinearElementField(law, trial, kernel.field->active_components, trial_ids, field_ids, u, field_values);
        }
        else if (basis)
            with_matrix_parameter(values, kernel.parameter, [&](const auto &value)
                                  { pytec3d::kernel::NodalField(value, trial, *basis, context.weights(), trial_ids, field_ids, target); });
        else
            with_matrix_parameter(values, kernel.parameter, [&](const auto &value)
                                  { pytec3d::kernel::ElementField(value, trial, kernel.field->active_components, trial_ids, field_ids, target); });
    }

    for (unsigned int i = 0; i < F.m(); ++i)
        for (unsigned int j = 0; j < F.n(); ++j)
            field_values(i) += F(i, j) * u(j) + F_dt(i, j) * u_dt(j);

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

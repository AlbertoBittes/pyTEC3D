/**
 * @file element_builder.cpp
 * @brief ElementBuilder: gathers per-element data and delegates matrices, vectors and fields to ElementKernels.
 */
#include "element_builder.h"

#include <algorithm>
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

    bh_curves_ = &equation_systems_->get_bh_curves();
    materials_rows_ = &equation_systems_->get_material_rows();
    constants_rows_ = &equation_systems_->get_constant_rows();
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

const std::vector<Real> &ElementBuilder::material_row(dof_id_type id) const
{
    const auto it = materials_rows_->find(id);
    if (it == materials_rows_->end())
        pytec3d::logging::error_and_abort("ElementBuilder: material " + std::to_string(id) + " is not defined.");
    return it->second;
}

const std::vector<Real> &ElementBuilder::constant_row(dof_id_type id) const
{
    const auto it = constants_rows_->find(id);
    if (it == constants_rows_->end())
        pytec3d::logging::error_and_abort("ElementBuilder: real constants " + std::to_string(id) + " are not defined.");
    return it->second;
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

pytec3d::ElementParameters ElementBuilder::element_parameters() const
{
    pytec3d::ElementParameters parameters;

    const auto &material = material_row(elem_data.material_num);
    parameters.material.relative_permeability = {material[MAT_MUR_X], material[MAT_MUR_Y], material[MAT_MUR_Z]};
    parameters.material.resistivity = {material[MAT_RHO_X], material[MAT_RHO_Y], material[MAT_RHO_Z]};
    parameters.material.magnetization = {material[MAT_M_X], material[MAT_M_Y], material[MAT_M_Z]};

    const auto &constants = constant_row(elem_data.real_constant_num);
    parameters.constants.thickness = constants[CONST_THK];
    parameters.constants.cross_section = constants[CONST_CS];
    parameters.constants.turns = constants[CONST_TURNS];
    parameters.constants.volume = constants[CONST_VOL];
    parameters.constants.radius = constants[CONST_RAD];
    parameters.constants.resistance = constants[CONST_R];
    parameters.constants.symmetry = constants[CONST_S];
    parameters.constants.current_direction = {constants[CONST_T_X], constants[CONST_T_Y], constants[CONST_T_Z]};

    if (element_local_to_global_rotations_)
    {
        const auto rotation = element_local_to_global_rotations_->find(static_cast<dof_id_type>(elem_data.elem->id()));
        if (rotation != element_local_to_global_rotations_->end())
            parameters.local_to_global = rotation->second;
    }

    if (bh_curves_)
    {
        const auto curve = bh_curves_->find(elem_data.material_num);
        if (curve != bh_curves_->end())
            parameters.bh_curve = &curve->second;
    }

    parameters.dofs = element_dofs();
    return parameters;
}

const std::vector<VariablesMap::SystemVariable> &ElementBuilder::system_variables(pyTEC3DSystemBase::SystemType system_type) const
{
    const auto key = std::make_pair(elem_data.elem->subdomain_id(), static_cast<int>(system_type));
    auto it = system_variables_cache_.find(key);
    if (it == system_variables_cache_.end())
        it = system_variables_cache_.emplace(key, VariablesMap::get_system_variables_for_element(*elem_data.entry->type, system_type)).first;
    return it->second;
}

std::vector<dof_id_type> ElementBuilder::column_dofs(const libMesh::System &system,
                                                     pyTEC3DSystemBase::SystemType system_type,
                                                     const pytec3d::VariableSpec &spec) const
{
    std::vector<std::vector<dof_id_type>> per_variable;
    for (const auto &variable : system_variables(system_type))
    {
        if (variable.variable != spec.name)
            continue;
        std::vector<dof_id_type> variable_dofs;
        system.get_dof_map().dof_indices(elem_data.elem, variable_dofs,
                                         system.variable_number(VariablesMap::get_info_by_key(variable.key).global_label));
        per_variable.push_back(std::move(variable_dofs));
    }
    if (per_variable.empty())
        pytec3d::logging::error_and_abort("ElementBuilder: variable '" + spec.name + "' is not in system '" + system.name() + "'.");
    if (per_variable.size() == 1)
        return per_variable.front();

    // One scalar libMesh variable per component; shape functions run node by node over the components.
    const std::size_t n_components = per_variable.size();
    const std::size_t n_nodes = per_variable.front().size();
    std::vector<dof_id_type> columns(n_components * n_nodes);
    for (std::size_t component = 0; component < n_components; ++component)
    {
        if (per_variable[component].size() != n_nodes)
            pytec3d::logging::error_and_abort("ElementBuilder: components of '" + spec.name + "' have different DOF counts.");
        for (std::size_t node = 0; node < n_nodes; ++node)
            columns[node * n_components + component] = per_variable[component][node];
    }
    return columns;
}

pytec3d::ElementDofs ElementBuilder::element_dofs() const
{
    std::vector<dof_id_type> element_dof_ids;
    dof_solution_dof_map_->dof_indices(elem_data.elem, element_dof_ids);

    std::unordered_map<dof_id_type, unsigned int> position;
    for (unsigned int i = 0; i < element_dof_ids.size(); ++i)
        position[element_dof_ids[i]] = i;

    pytec3d::ElementDofs dofs;
    dofs.n_dofs = static_cast<unsigned int>(element_dof_ids.size());
    std::size_t covered = 0;
    for (const pytec3d::VariableSpec &spec : elem_data.entry->type->variables())
    {
        if (spec.kind != pytec3d::FieldKind::Dof)
            continue;
        std::vector<unsigned int> &columns = dofs.columns[spec.name];
        for (const dof_id_type dof : column_dofs(*dof_solution_system_, pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM, spec))
        {
            const auto it = position.find(dof);
            if (it == position.end())
                pytec3d::logging::error_and_abort("ElementBuilder: DOF of variable '" + spec.name + "' is not a DOF of the element.");
            columns.push_back(it->second);
        }
        covered += columns.size();
    }
    if (covered != element_dof_ids.size())
        pytec3d::logging::error_and_abort("ElementBuilder: element variables cover " + std::to_string(covered) + " of " +
                                          std::to_string(element_dof_ids.size()) + " element DOFs.");
    return dofs;
}

pytec3d::ElementLoads ElementBuilder::element_loads() const
{
    pytec3d::ElementLoads loads;

    for (const std::string &name : elem_data.entry->kernels->element_load_names())
    {
        vtkDataArray *array = vtk_grid_->GetCellData()->GetArray(name.c_str());
        if (!array)
            continue;
        const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
        std::array<Real, 3> value{0.0, 0.0, 0.0};
        for (int c = 0; c < std::min(array->GetNumberOfComponents(), 3); ++c)
            value[c] = array->GetComponent(elem_id, c);
        loads.element[name] = value;
    }

    for (const pytec3d::VariableSpec &variable : elem_data.entry->type->variables())
    {
        if (variable.kind != pytec3d::FieldKind::Dof)
            continue;

        if (!variable.nodal_loads.empty())
        {
            vtkDataArray *array = vtk_grid_->GetPointData()->GetArray(variable.nodal_loads.c_str());
            if (!array || !dof_solution_system_->has_variable(variable.name))
                continue;

            const unsigned int var_num = dof_solution_system_->variable_number(variable.name);
            std::vector<dof_id_type> variable_dofs;
            dof_solution_dof_map_->dof_indices(elem_data.elem, variable_dofs, var_num);

            pytec3d::ElementLoads::Nodal nodal;
            for (unsigned int i = 0; i < variable_dofs.size(); ++i)
            {
                const libMesh::Node &node = elem_data.elem->node_ref(i);
                const dof_id_type dof_id = node.dof_number(dof_solution_system_->number(), var_num, 0);
                nodal.values.push_back(array->GetComponent(static_cast<vtkIdType>(node.id()), 0));
                nodal.inject.push_back(should_inject_nodal_load_dof(variable.name, dof_id));
            }
            loads.nodal[variable.name] = std::move(nodal);
        }
    }
    return loads;
}

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
    DenseMatrix<Number> stiffness;
    DenseMatrix<Number> time_derivative;
    elem_data.entry->kernels->build_matrices(*elem_data.elem, element_parameters(), stiffness, time_derivative);

    if (stiffness.m() != n_dofs)
        pytec3d::logging::error_and_abort("ElementBuilder: kernel DOF layout (" + std::to_string(stiffness.m()) +
                                          ") does not match the element DOFs (" + std::to_string(n_dofs) + ").");

    matrices[MatricesSet::DOF_SOLUTION_MATRICES] = {std::move(stiffness), std::move(time_derivative)};
}

void ElementBuilder::build_element_rhs()
{
    DenseVector<Number> rhs;
    elem_data.entry->kernels->build_rhs(*elem_data.elem, element_parameters(), element_loads(), rhs);
    vectors[VectorsSet::RHS_VECTORS] = {std::move(rhs)};
}

void ElementBuilder::build_element_tangent_stiffness_residual()
{
    DenseMatrix<Number> tangent;
    DenseVector<Number> residual;
    elem_data.entry->kernels->build_tangent_residual(*elem_data.elem,
                                                     element_parameters(),
                                                     local_dof_vector(dof_solution_system_->current_local_solution.get()),
                                                     tangent,
                                                     residual);

    matrices[MatricesSet::TANGENT_STIFFNESS_MATRICES] = {std::move(tangent)};
    vectors[VectorsSet::RESIDUAL_VECTORS] = {std::move(residual)};
}

void ElementBuilder::build_nodal_solution()
{
    const auto fields = elem_data.entry->kernels->nodal_fields(
        *elem_data.elem,
        element_parameters(),
        local_dof_vector(dof_solution_system_->current_local_solution.get()),
        local_dof_vector(dof_solution_system_->request_vector("solution_time_derivative")));
    write_fields(*nodal_solution_system_, pyTEC3DSystemBase::NODAL_SOLUTION, pytec3d::FieldKind::Nodal, fields);
}

void ElementBuilder::build_element_solution()
{
    const auto fields = elem_data.entry->kernels->element_fields(
        *elem_data.elem,
        element_parameters(),
        local_dof_vector(dof_solution_system_->current_local_solution.get()),
        local_dof_vector(dof_solution_system_->request_vector("solution_time_derivative")));
    write_fields(*element_solution_system_, pyTEC3DSystemBase::ELEMENT_SOLUTION, pytec3d::FieldKind::Element, fields);
}

void ElementBuilder::write_fields(System &system,
                                  pyTEC3DSystemBase::SystemType system_type,
                                  pytec3d::FieldKind kind,
                                  const std::vector<DenseVector<Number>> &fields)
{
    static std::mutex solution_write_mutex;

    std::size_t index = 0;
    for (const pytec3d::VariableSpec &spec : elem_data.entry->type->variables())
    {
        if (spec.kind != kind)
            continue;
        if (index >= fields.size())
            pytec3d::logging::error_and_abort("ElementBuilder: missing values for field '" + spec.name + "'.");
        const DenseVector<Number> &values = fields[index++];

        const std::vector<dof_id_type> dofs = column_dofs(system, system_type, spec);
        if (dofs.size() != values.size())
            pytec3d::logging::error_and_abort("ElementBuilder: field '" + spec.name + "' has " + std::to_string(values.size()) +
                                              " values but " + std::to_string(dofs.size()) + " DOFs in system '" + system.name() + "'.");

        std::lock_guard<std::mutex> lock(solution_write_mutex);
        for (unsigned int i = 0; i < dofs.size(); ++i)
            if ((system.solution->first_local_index() <= dofs[i]) && (dofs[i] < system.solution->last_local_index()))
                system.solution->set(dofs[i], values(i));
    }
}

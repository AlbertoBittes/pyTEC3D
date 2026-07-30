/**
 * @file element_builder.cpp
 * @brief ElementBuilder implementation: construction, element-data setup, matrix/vector sizing,
 *        dispatch to formulation-specific builders, and solution post-processing.
 */
#include "element_builder.h"
#include "binary_reader.h"
#include "coordinate_system.h"
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <string>
#include <pybind11/stl.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>

#include "element_types.h"
#include "apec_m_system.h"
#include "apec_m_solver.h"
#include "element_types_factory.h"

#include "variables_map.h"
#include "apec_m_logging.h"

#include <stdexcept>

namespace py = pybind11;

using libMesh::DenseVector;
using libMesh::dof_id_type;
using libMesh::DofMap;
using libMesh::Elem;
using libMesh::Number;
using libMesh::NumericVector;
using libMesh::Real;
using libMesh::System;

// ---------------------------------------------------------------------------
// Construction and initialisation
// ---------------------------------------------------------------------------

ElementBuilder::ElementBuilder(ApecMEquationSystems *equation_systems)
    : FEParentObjectBase<ApecMEquationSystems>(equation_systems),
      matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS)
{
    equation_systems_ = equation_systems;
    dof_solution_system_ = &equation_systems_->get_system<libMesh::ImplicitSystem>(
        ApecMSystemBase::system_name_from_type(ApecMSystemBase::DOF_SOLUTION_SYSTEM));

    nodal_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::NODAL_SOLUTION);
    element_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::ELEMENT_SOLUTION);
    dof_solution_dof_map_ = &dof_solution_system_->get_dof_map();

    vtk_grid_ = equation_systems_->get_parent_object()->get_active_vtk_grid();
    material_cell_array_ = vtk_grid_->GetCellData()->GetArray("MATERIAL");
    real_constant_cell_array_ = vtk_grid_->GetCellData()->GetArray("CONSTANT");
    coordinate_system_cell_array_ = vtk_grid_->GetCellData()->GetArray("COORDINATE_SYSTEM");

    const auto *solver = equation_systems_->get_parent_object();
    if (solver && solver->has_bh_curves())
        bh_curves_ = &solver->get_bh_curves();

    py::object apec_object = equation_systems_->get_parent_object()->get_apec_object();
    py::object materials_attribute = apec_object.attr("materials");
    py::object constants_attribute = apec_object.attr("constants");

    for (const py::handle row_obj : materials_attribute.attr("to_dict")("records"))
    {
        const py::dict row = py::reinterpret_borrow<py::dict>(row_obj);
        const dof_id_type id = py::cast<dof_id_type>(row["MATERIAL_ID"]);
        auto &values = materials_rows_[id];
        values.resize(MAT_NUM_COLS);
        values[MAT_MUR_X] = py::cast<Real>(row["MUR_X"]);
        values[MAT_MUR_Y] = py::cast<Real>(row["MUR_Y"]);
        values[MAT_MUR_Z] = py::cast<Real>(row["MUR_Z"]);
        values[MAT_RHO_X] = py::cast<Real>(row["RHO_X"]);
        values[MAT_RHO_Y] = py::cast<Real>(row["RHO_Y"]);
        values[MAT_RHO_Z] = py::cast<Real>(row["RHO_Z"]);
        values[MAT_M_X] = py::cast<Real>(row["M_X"]);
        values[MAT_M_Y] = py::cast<Real>(row["M_Y"]);
        values[MAT_M_Z] = py::cast<Real>(row["M_Z"]);
    }

    for (const py::handle row_obj : constants_attribute.attr("to_dict")("records"))
    {
        const py::dict row = py::reinterpret_borrow<py::dict>(row_obj);
        const dof_id_type id = py::cast<dof_id_type>(row["CONSTANTS_ID"]);
        auto &values = constants_rows_[id];
        values.resize(CONST_NUM_COLS);
        values[CONST_THK] = py::cast<Real>(row["THK"]);
        values[CONST_CS] = py::cast<Real>(row["CS"]);
        values[CONST_TURNS] = py::cast<Real>(row["TURNS"]);
        values[CONST_VOL] = py::cast<Real>(row["VOL"]);
        values[CONST_RAD] = py::cast<Real>(row["RAD"]);
        values[CONST_T_X] = py::cast<Real>(row["T_X"]);
        values[CONST_T_Y] = py::cast<Real>(row["T_Y"]);
        values[CONST_T_Z] = py::cast<Real>(row["T_Z"]);
        values[CONST_R] = py::cast<Real>(row["R"]);
        values[CONST_S] = py::cast<Real>(row["S"]);
    }

    const py::object coordinate_systems_object = apec_object.attr("coordinate_systems");
    if (!coordinate_systems_object.is_none())
    {
        const py::dict coordinate_systems_dict = py::cast<py::dict>(coordinate_systems_object);
        for (const auto &item : coordinate_systems_dict)
        {
            const dof_id_type id = py::cast<dof_id_type>(item.first);
            coordinate_systems_.emplace(
                id,
                py::cast<std::shared_ptr<apec_m::CoordinateSystem>>(item.second));
        }
    }
}

//
//
// Methods to size element matrices and vectors
//
//

void ElementBuilder::build_vectors(std::vector<std::vector<unsigned int>> *sizes_map, unsigned int vectors_set)
{
    this->vectors[vectors_set].resize(sizes_map->size());
    for (unsigned int i = 0; i < sizes_map->size(); ++i)
    {
        const auto &size_pair = (*sizes_map)[i];
        unsigned int rows = size_pair[0];
        this->vectors[vectors_set][i].resize(rows);
    }
}

void ElementBuilder::build_matrices(std::vector<std::vector<unsigned int>> *sizes_map, unsigned int matrices_set)
{

    unsigned int num_matrices = static_cast<unsigned int>(sizes_map->size());
    this->matrices[matrices_set].resize(num_matrices);
    for (unsigned int i = 0; i < num_matrices; ++i)
    {
        const auto &size_pair = (*sizes_map)[i];
        unsigned int rows = size_pair[0];
        unsigned int cols = size_pair[1];
        this->matrices[matrices_set][i].resize(rows, cols);
    }
}

void ElementBuilder::prepare_element_matrix_building()
{
    // Two square DOF matrices (Ke stiffness, Ce capacitance) and one RHS vector.
    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs},
        {n_dofs, n_dofs}};

    this->build_matrices(&matrix_sizes, MatricesSet::DOF_SOLUTION_MATRICES);
    this->build_matrices(elem_data.element_type->get_nodal_solution_matrix_idx_size_map().get(), MatricesSet::NODAL_SOLUTION_MATRICES);
    this->build_matrices(elem_data.element_type->get_element_solution_matrix_idx_size_map().get(), MatricesSet::ELEMENT_SOLUTION_MATRICES);
}

// ---------------------------------------------------------------------------
// Element-data setup
// ---------------------------------------------------------------------------

void ElementBuilder::set_element_data(Elem *elem)
{
    auto element_types = equation_systems_->get_parent_object()->get_element_types_factory();

    elem_data.element_type = element_types->get_element_type(elem->subdomain_id());

    const vtkIdType elem_id = static_cast<vtkIdType>(elem->id());
    elem_data.material_num = static_cast<unsigned int>(material_cell_array_->GetComponent(elem_id, 0));
    elem_data.real_constant_num = static_cast<unsigned int>(real_constant_cell_array_->GetComponent(elem_id, 0));
    elem_data.coordinate_system_num = static_cast<unsigned int>(coordinate_system_cell_array_->GetComponent(elem_id, 0));
    elem_data.elem = elem;
    elem_data.coordinate_system = nullptr;

    if (elem_data.coordinate_system_num != 0)
    {
        const auto it = coordinate_systems_.find(static_cast<dof_id_type>(elem_data.coordinate_system_num));
        elem_data.coordinate_system = it->second.get();
    }

    std::vector<dof_id_type> dof_indices;
    dof_solution_dof_map_->dof_indices(elem, dof_indices);
    this->n_dofs = static_cast<unsigned int>(dof_indices.size());
}

ElemData &ElementBuilder::get_elem_data()
{
    return elem_data;
}

const ElemData &ElementBuilder::get_elem_data() const
{
    return elem_data;
}

bool ElementBuilder::nonlinear_element() const
{
    if (!bh_curves_)
        return false;

    const dof_id_type material_id = static_cast<dof_id_type>(elem_data.material_num);
    return bh_curves_->find(material_id) != bh_curves_->end();
}

void ElementBuilder::rotate_vector_for_elem(const Real *local_components,
                                            Real *global_components) const
{
    if (!elem_data.coordinate_system)
    {
        global_components[0] = local_components[0];
        global_components[1] = local_components[1];
        global_components[2] = local_components[2];
        return;
    }

    const libMesh::Point center = elem_data.elem->vertex_average();
    const libMesh::Real position[3] = {
        static_cast<libMesh::Real>(center(0)),
        static_cast<libMesh::Real>(center(1)),
        static_cast<libMesh::Real>(center(2))};
    elem_data.coordinate_system->vector_components_to_global_raw(local_components, position, 1, global_components);
}

void ElementBuilder::rotate_diag_tensor_for_elem(const Real diag[3],
                                                 Real global_tensor[3][3]) const
{
    if (!elem_data.coordinate_system)
    {
        global_tensor[0][0] = diag[0];
        global_tensor[0][1] = 0.0;
        global_tensor[0][2] = 0.0;
        global_tensor[1][0] = 0.0;
        global_tensor[1][1] = diag[1];
        global_tensor[1][2] = 0.0;
        global_tensor[2][0] = 0.0;
        global_tensor[2][1] = 0.0;
        global_tensor[2][2] = diag[2];
        return;
    }

    const libMesh::Point center = elem_data.elem->vertex_average();
    const libMesh::Real position[3] = {
        static_cast<libMesh::Real>(center(0)),
        static_cast<libMesh::Real>(center(1)),
        static_cast<libMesh::Real>(center(2))};
    const libMesh::Real local_raw[9] = {
        diag[0], 0.0, 0.0,
        0.0, diag[1], 0.0,
        0.0, 0.0, diag[2]};
    libMesh::Real global_raw[9] = {0.0};
    elem_data.coordinate_system->rotate_tensor_to_global_raw(local_raw, position, global_raw);
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            global_tensor[i][j] = global_raw[i * 3 + j];
}

// ---------------------------------------------------------------------------
// Element matrix and RHS build dispatch
// ---------------------------------------------------------------------------

void ElementBuilder::build_element()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_OPTIONs()[0]);

    switch (etype)
    {
    case ElementTypeEnum::AV3DNODAL:
        this->build_av3dnodal();
        break;
    case ElementTypeEnum::AV3DEDGE:
        this->build_av3dedge();
        break;
    case ElementTypeEnum::AV2DNODAL:
        this->build_av2dnodal();
        break;
    default:
        apec::logging::error("Unsupported element type enum={} while building element", static_cast<unsigned int>(etype));
        apec::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

// ---------------------------------------------------------------------------
// Tangent stiffness matrix and RHS build dispatch
// ---------------------------------------------------------------------------

void ElementBuilder::build_element_tangent_stiffness_residual()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_OPTIONs()[0]);

    switch (etype)
    {
    case ElementTypeEnum::AV3DNODAL:
        this->build_av3dnodal_tangent_stiffness_residual();
        break;
    case ElementTypeEnum::AV3DEDGE:
        this->build_av3dedge_tangent_stiffness_residual();
        break;
    case ElementTypeEnum::AV2DNODAL:
        this->build_av2dnodal_tangent_stiffness_residual();
        break;
    default:
        apec::logging::error("Unsupported element type enum={} while building element", static_cast<unsigned int>(etype));
        apec::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_element_rhs()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_OPTIONs()[0]);

    switch (etype)
    {
    case ElementTypeEnum::AV3DNODAL:
        this->build_av3dnodal_rhs();
        break;
    case ElementTypeEnum::AV3DEDGE:
        this->build_av3dedge_rhs();
        break;
    case ElementTypeEnum::AV2DNODAL:
        this->build_av2dnodal_rhs();
        break;
    default:
        apec::logging::error("Unsupported element type enum={} while building element RHS", static_cast<unsigned int>(etype));
        apec::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

// ---------------------------------------------------------------------------
// Post-processing: nodal and element solution fields
// ---------------------------------------------------------------------------

void ElementBuilder::build_nodal_solution()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_OPTIONs()[0]);

    switch (etype)
    {
    case ElementTypeEnum::AV3DNODAL:
        this->build_av3dnodal_nodal_solution();
        break;
    case ElementTypeEnum::AV3DEDGE:
        this->build_av3dedge_nodal_solution();
        break;
    case ElementTypeEnum::AV2DNODAL:
        this->build_av2dnodal_nodal_solution();
        break;
    default:
        apec::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_element_solution()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_OPTIONs()[0]);

    switch (etype)
    {
    case ElementTypeEnum::AV3DNODAL:
        this->build_av3dnodal_element_solution();
        break;
    case ElementTypeEnum::AV3DEDGE:
        this->build_av3dedge_element_solution();
        break;
    case ElementTypeEnum::AV2DNODAL:
        this->build_av2dnodal_element_solution();
        break;
    default:
        apec::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

std::vector<DenseVector<Number>> ElementBuilder::get_variables_solution(const NumericVector<Number> *solution,
                                                                        int dof_map_index)
{
    System *sys = &equation_systems_->get_system(dof_map_index);
    DofMap *dof_map = (dof_map_index == ApecMSystemBase::DOF_SOLUTION_SYSTEM) ? dof_solution_dof_map_ : &sys->get_dof_map();
    // Cast dof_map_index to SystemType for SYSTEM_INFO_MAP lookup
    ApecMSystemBase::SystemType sys_type = static_cast<ApecMSystemBase::SystemType>(dof_map_index);

    std::vector<DenseVector<Number>> sol_vectors;
    if (!solution || !dof_map)
        return sol_vectors;

    const auto vars = VariablesMap::get_active_variables_for_element(*elem_data.element_type, sys_type);
    sol_vectors.reserve(vars.size());

    for (const auto &var_enum : vars)
    {
        unsigned int system_var_num = sys->variable_number(VariablesMap::get_info_by_key(var_enum).global_label);
        std::vector<dof_id_type> dof_indices_variable;
        dof_map->dof_indices(elem_data.elem, dof_indices_variable, system_var_num);
        DenseVector<Number> local_sol(dof_indices_variable.size());
        for (unsigned int i = 0; i < dof_indices_variable.size(); ++i)
            local_sol(i) = solution->el(dof_indices_variable[i]);
        sol_vectors.emplace_back(std::move(local_sol));
    }

    return sol_vectors;
}

void ElementBuilder::append_and_write_solution(System &system,
                                               const std::vector<DenseVector<Number>> &ordered_vectors)
{
    DenseVector<Number> values;
    for (const auto &vec : ordered_vectors)
    {
        values.append(vec);
    }

    if (values.size() == 0)
        return;

    std::vector<dof_id_type> dof_indices;
    DofMap *dof_map = &system.get_dof_map();
    dof_map->dof_indices(elem_data.elem, dof_indices);

    for (unsigned int i = 0; i < dof_indices.size() && i < values.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((system.solution->first_local_index() <= dof_index) &&
            (dof_index < system.solution->last_local_index()))
            system.solution->set(dof_index, values(i));
    }
}

#include "element_builder.h"
#include "binary_reader.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <string>
#include <pybind11/stl.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>

// Define Gauss quadrature rules.
#include "element_types.h" // Include the header defining ElementType
#include "apec_m_system.h"
#include "apec_m_solver.h"
#include "element_types_factory.h"

#include "variables_map.h"

namespace py = pybind11;

//
//
// constructor and init
//
//

ElementBuilder::ElementBuilder(EquationSystems *equation_systems)
    : FEParentObjectBase<EquationSystems>(equation_systems),
      matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS)
{
    equation_systems_ = dynamic_cast<ApecMEquationSystems *>(equation_systems);
    dof_solution_system_ = dynamic_cast<ApecMLinearSystem *>(&equation_systems_->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM));

    nodal_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::NODAL_SOLUTION);
    element_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::ELEMENT_SOLUTION);
    dof_solution_dof_map_ = &dof_solution_system_->get_dof_map();

    vtk_grid_ = equation_systems_->get_parent_object()->active_vtk_grid;
    material_cell_array_ = vtk_grid_->GetCellData()->GetArray("MATERIAL");
    real_constant_cell_array_ = vtk_grid_->GetCellData()->GetArray("CONSTANT");

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
}

const std::vector<Real> &ElementBuilder::material_row(dof_id_type id) const
{
    const auto it = materials_rows_.find(id);
    if (it == materials_rows_.end())
        libmesh_error_msg("Material ID not found in cached materials map: " + std::to_string(id));
    return it->second;
}

const std::vector<Real> &ElementBuilder::constant_row(dof_id_type id) const
{
    const auto it = constants_rows_.find(id);
    if (it == constants_rows_.end())
        libmesh_error_msg("Constant ID not found in cached constants map: " + std::to_string(id));
    return it->second;
}

Real ElementBuilder::material_value(dof_id_type id, unsigned int col_idx) const
{
    const auto &row = material_row(id);
    if (col_idx >= row.size())
        libmesh_error_msg("Material cache column index out of range: " + std::to_string(col_idx));
    return row[col_idx];
}

Real ElementBuilder::constant_value(dof_id_type id, unsigned int col_idx) const
{
    const auto &row = constant_row(id);
    if (col_idx >= row.size())
        libmesh_error_msg("Constants cache column index out of range: " + std::to_string(col_idx));
    return row[col_idx];
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
    // reinit to be done inside each element
    // as well as the reception of shape functions
    // this->elem_data.element_type->reinit_fe(this->elem);
    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs},
        {n_dofs, n_dofs}};

    // One vector of size n_dofs
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};

    // Usage:
    this->build_matrices(&matrix_sizes, MatricesSet::DOF_SOLUTION_MATRICES);
    this->build_vectors(&vector_sizes, VectorsSet::RHS_VECTORS);
    this->build_matrices(elem_data.element_type->get_solution_matrix_idx_size_map(true).get(), MatricesSet::NODAL_SOLUTION_MATRICES);
    this->build_matrices(elem_data.element_type->get_solution_matrix_idx_size_map(false).get(), MatricesSet::ELEMENT_SOLUTION_MATRICES);
    // this->build_vectors(vector_idx_size_map, 1);
}

//
//
// Setters
//
//

void ElementBuilder::set_element_data(Elem *elem)
{
    auto element_types = equation_systems_->get_parent_object()->element_types_factory;

    elem_data.element_type = element_types->get_element_type(elem->subdomain_id());

    if (!elem->infinite())
    {
        const vtkIdType elem_id = static_cast<vtkIdType>(elem->id());
        elem_data.material_num = static_cast<unsigned int>(material_cell_array_->GetComponent(elem_id, 0));
        elem_data.real_constant_num = static_cast<unsigned int>(real_constant_cell_array_->GetComponent(elem_id, 0));
    }
    elem_data.elem = elem;

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

//
//
// Methods to build element matrices and rhs
//
//

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
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

//
//
// solution calculation methods
//
//

void ElementBuilder::build_nodal_solution()
{
    if (!this->compute_nodal_solution_vectors)
        return;

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
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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

    // This is now maybe too slow, but it works for now
    // TODO: Check if there is room for optimization later
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

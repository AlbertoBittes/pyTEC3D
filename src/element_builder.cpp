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
#include <array>
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

namespace
{
    const Real kMu0 = 4.0 * libMesh::pi * 1e-7;
}

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

    apec_object_ = equation_systems_->get_parent_object()->get_apec_object();
    vtk_grid_ = equation_systems_->get_parent_object()->active_vtk_grid;
    materials_attribute_ = apec_object_.attr("materials");
    real_constants_attribute_ = apec_object_.attr("real_constants");
    material_cell_array_ = vtk_grid_->GetCellData()->GetArray("MATERIAL");
    real_constant_cell_array_ = vtk_grid_->GetCellData()->GetArray("CONSTANT");
}

py::object ElementBuilder::get_row_by_id(const py::object &table, const char *id_column, unsigned int id) const
{
    const py::object mask = table.attr("__getitem__")(id_column).attr("__eq__")(id);
    return table.attr("loc").attr("__getitem__")(mask).attr("iloc").attr("__getitem__")(0);
}

std::array<Real, 3> ElementBuilder::get_mu_matrix(unsigned int material_num) const
{
    const py::object row = get_row_by_id(materials_attribute_, "MATERIAL", material_num);
    const Real murx = py::cast<Real>(row.attr("__getitem__")("MURX"));
    const Real mury = py::cast<Real>(row.attr("__getitem__")("MURY"));
    const Real murz = py::cast<Real>(row.attr("__getitem__")("MURZ"));
    return {kMu0 * murx, kMu0 * mury, kMu0 * murz};
}

std::array<Real, 3> ElementBuilder::get_sigma_matrix(unsigned int material_num) const
{
    const py::object row = get_row_by_id(materials_attribute_, "MATERIAL", material_num);
    const Real rsvx = py::cast<Real>(row.attr("__getitem__")("RSVX"));
    const Real rsvy = py::cast<Real>(row.attr("__getitem__")("RSVY"));
    const Real rsvz = py::cast<Real>(row.attr("__getitem__")("RSVZ"));
    return {1.0 / rsvx, 1.0 / rsvy, 1.0 / rsvz};
}

std::vector<Real> ElementBuilder::get_real_constants(unsigned int constant_num) const
{
    const py::object row = get_row_by_id(real_constants_attribute_, "CONSTANT", constant_num);
    static const std::array<const char *, 20> columns = {
        "THK", "SC", "NC", "VC", "RAD",
        "TX", "TY", "TZ", "R", "SYM",
        "VELOX", "VELOY", "VELOZ",
        "OMEGAX", "OMEGAY", "OMEGAZ",
        "XLOC", "YLOC", "ZLOC", "CPSIGN"};

    std::vector<Real> values;
    values.reserve(columns.size());
    for (const char *col : columns)
        values.push_back(py::cast<Real>(row.attr("__getitem__")(col)));
    return values;
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
    case ElementTypeEnum::SOLID97:
        this->build_solid97();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233();
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
    case ElementTypeEnum::SOLID97:
        this->build_solid97_rhs();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_rhs();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_rhs();
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
    case ElementTypeEnum::SOLID97:
        this->build_solid97_nodal_solution();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_nodal_solution();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_nodal_solution();
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
    case ElementTypeEnum::SOLID97:
        this->build_solid97_element_solution();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_element_solution();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_element_solution();
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

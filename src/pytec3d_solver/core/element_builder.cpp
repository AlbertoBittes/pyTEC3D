/**
 * @file element_builder.cpp
 * @brief ElementBuilder implementation: construction, element-data setup, matrix/vector sizing,
 *        dispatch to formulation-specific builders, and solution post-processing.
 */
#include "element_builder.h"
#include <libmesh/dof_map.h>
#include <libmesh/elem_range.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/vector_value.h>
#include <string>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include "element_types.h"
#include "pytec3d_system.h"
#include "pytec3d_model.h"
#include "element_types_factory.h"

#include "variables_map.h"
#include "logging.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>

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

void ElementBuilder::set_element_data(const Elem *elem)
{
    elem_data.element_type = element_types_factory_->get_element_type(elem->subdomain_id());

    const vtkIdType elem_id = static_cast<vtkIdType>(elem->id());
    elem_data.material_num = static_cast<unsigned int>(material_cell_array_->GetComponent(elem_id, 0));
    elem_data.real_constant_num = static_cast<unsigned int>(real_constant_cell_array_->GetComponent(elem_id, 0));
    elem_data.coordinate_system_num = static_cast<unsigned int>(coordinate_system_cell_array_->GetComponent(elem_id, 0));
    elem_data.elem = elem;

    std::vector<dof_id_type> dof_indices;
    dof_solution_dof_map_->dof_indices(elem, dof_indices);
    this->n_dofs = static_cast<unsigned int>(dof_indices.size());
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

void ElementBuilder::rotate_vector_for_elem(libMesh::RealVectorValue &components) const
{
    if (!elem_data.elem || !element_local_to_global_rotations_)
        return;

    const auto rot_it = element_local_to_global_rotations_->find(static_cast<dof_id_type>(elem_data.elem->id()));
    if (rot_it == element_local_to_global_rotations_->end())
        return;

    const auto &r = rot_it->second;
    components = r * components;
}

void ElementBuilder::rotate_diag_tensor_for_elem(libMesh::TensorValue<Real> &tensor) const
{
    if (!elem_data.elem || !element_local_to_global_rotations_)
        return;

    const auto rot_it = element_local_to_global_rotations_->find(static_cast<dof_id_type>(elem_data.elem->id()));
    if (rot_it == element_local_to_global_rotations_->end())
        return;

    const auto &r = rot_it->second;
    tensor = r * tensor * r.transpose();
}

void ElementBuilder::evaluate_isotropic_bh(const BHCurve &curve, const Real B[3], unsigned int dim,
                                           Real H[3], Real jac[3][3], Real *nu_out) const
{
    Real b_mag_sq = 0.0;
    for (unsigned int k = 0; k < dim; ++k)
        b_mag_sq += B[k] * B[k];
    const Real b_mag = std::sqrt(b_mag_sq);

    const Real nu = curve.evaluate_reluctivity_from_b(b_mag);
    const Real dhdb = curve.slope_h_vs_b(b_mag);
    if (nu_out)
        *nu_out = nu;

    for (unsigned int k = 0; k < 3; ++k)
    {
        H[k] = 0.0;
        for (unsigned int l = 0; l < 3; ++l)
            jac[k][l] = 0.0;
    }

    // Isotropic differential reluctivity tensor: nu(B)*I + (dH/dB - nu(B)) * (B tensor B)/B^2.
    constexpr Real kZeroB = 1e-14;
    for (unsigned int k = 0; k < dim; ++k)
    {
        H[k] = nu * B[k];
        jac[k][k] += nu;
        if (b_mag > kZeroB)
            for (unsigned int l = 0; l < dim; ++l)
                jac[k][l] += (dhdb - nu) * B[k] * B[l] / b_mag_sq;
    }
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
        pytec3d::logging::error("Unsupported element type enum={} while building element", static_cast<unsigned int>(etype));
        pytec3d::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
        pytec3d::logging::error("Unsupported element type enum={} while building element", static_cast<unsigned int>(etype));
        pytec3d::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
        pytec3d::logging::error("Unsupported element type enum={} while building element RHS", static_cast<unsigned int>(etype));
        pytec3d::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
        pytec3d::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
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
        pytec3d::logging::error_and_abort("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

std::vector<DenseVector<Number>> ElementBuilder::get_variables_solution(const NumericVector<Number> *solution,
                                                                        int dof_map_index)
{
    System *sys = &equation_systems_->get_system(dof_map_index);
    DofMap *dof_map = (dof_map_index == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM) ? dof_solution_dof_map_ : &sys->get_dof_map();
    // Cast dof_map_index to SystemType for SYSTEM_INFO_MAP lookup
    pyTEC3DSystemBase::SystemType sys_type = static_cast<pyTEC3DSystemBase::SystemType>(dof_map_index);

    std::vector<DenseVector<Number>> sol_vectors;
    if (!solution || !dof_map)
        return sol_vectors;

    std::vector<unsigned int> vars_numbers;
    sol_vectors.reserve(sys->n_vars());

    for (unsigned int var_num = 0; var_num < sys->n_vars(); ++var_num)
    {
        auto variable = dof_map->variable(var_num);
        if (!variable.active_on_subdomain(elem_data.elem->subdomain_id()))
        {
            sol_vectors.emplace_back(DenseVector<Number>(0));
            continue;
        }
        std::vector<dof_id_type> dof_indices_variable;
        dof_map->dof_indices(elem_data.elem, dof_indices_variable, var_num);
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
    static std::mutex solution_write_mutex;

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

    std::lock_guard<std::mutex> lock(solution_write_mutex);

    for (unsigned int i = 0; i < dof_indices.size() && i < values.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((system.solution->first_local_index() <= dof_index) &&
            (dof_index < system.solution->last_local_index()))
            system.solution->set(dof_index, values(i));
    }
}

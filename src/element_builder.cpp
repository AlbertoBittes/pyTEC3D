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
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include "element_types.h"
#include "apec_m_system.h"
#include "apec_m_model.h"
#include "element_types_factory.h"

#include "variables_map.h"
#include "apec_m_logging.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
    std::filesystem::path get_emat_folderpath_from_solver()
    {
        return std::filesystem::current_path() / "emat";
    }
}

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

    element_types_factory_ = equation_systems_->get_element_types_factory().get();
    bh_curves_ = &equation_systems_->get_bh_curves();
    materials_rows_ = &equation_systems_->get_material_rows();
    constants_rows_ = &equation_systems_->get_constant_rows();
    coordinate_systems_ = &equation_systems_->get_coordinate_systems();
}

void ElementBuilder::write_element_data()
{
    this->write_element_matrices();
    this->write_element_rhs();
}

void ElementBuilder::write_element_matrices()
{
    const libMesh::MeshBase &mesh = equation_systems_->get_mesh();
    const auto &comm = mesh.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;

    const std::vector<std::string> filenames = {
        "mtrx_" + rank_str + ".emat",
        "nodal_sol_mtrx_" + rank_str + ".emat",
        "element_sol_mtrx_" + rank_str + ".emat"};

    const bool emat_files_exist = std::all_of(filenames.begin(), filenames.end(), [&](const std::string &filename)
                                              { return std::filesystem::exists(emat_folderpath / filename); });

    if (!emat_files_exist)
    {
        apec::logging::info("Computing element matrices for current time-step");

        auto streams = reader.open_output_streams(emat_folderpath, filenames);
        auto &mtrx_output_file = *streams[0];
        auto &nodal_sol_mtrx_output_file = *streams[1];
        auto &element_sol_mtrx_output_file = *streams[2];

        this->run_for_elements(mesh, [&](const libMesh::Elem *)
                               {
                                   this->build_element();

                                   BinaryReader::write_element_matrices(*this, mtrx_output_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
                                   BinaryReader::write_element_matrices(*this, nodal_sol_mtrx_output_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES);
                                   BinaryReader::write_element_matrices(*this, element_sol_mtrx_output_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES); });

        reader.close_file_streams(streams);
    }
    else
        apec::logging::info("Reusing element matrices from previous time step");
}

void ElementBuilder::build_nodal_load_injector_map(const libMesh::MeshBase &mesh,
                                                   const std::string &dof_variable_name,
                                                   const std::string &nodal_load_array_name)
{
    if (!dof_solution_system_->has_variable(dof_variable_name))
        return;

    vtkDataArray *load_array = vtk_grid_->GetPointData()->GetArray(nodal_load_array_name.c_str());
    if (!load_array)
        return;

    const unsigned int dof_var_num = dof_solution_system_->variable_number(dof_variable_name);
    const unsigned int sys_num = dof_solution_system_->number();
    const auto &dof_var = dof_solution_system_->variable(dof_var_num);

    std::map<dof_id_type, dof_id_type> local_injector_elem_by_dof;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        if (!dof_var.active_on_subdomain(elem->subdomain_id()))
            continue;

        for (unsigned int i = 0; i < elem->n_nodes(); ++i)
        {
            const auto &node = elem->node_ref(i);
            const auto node_id = static_cast<vtkIdType>(node.id());
            const Real load_value = load_array->GetComponent(node_id, 0);
            if (load_value == 0.0)
                continue;

            const dof_id_type dof_id = node.dof_number(sys_num, dof_var_num, 0);
            local_injector_elem_by_dof.emplace(dof_id, static_cast<dof_id_type>(elem->id()));
        }
    }

    auto global_injector_elem_by_dof = local_injector_elem_by_dof;
    mesh.comm().set_union(global_injector_elem_by_dof);
    auto &injector_elem_by_dof = nodal_load_injector_elem_by_dof_by_variable[dof_variable_name];
    for (const auto &entry : global_injector_elem_by_dof)
        injector_elem_by_dof[entry.first] = entry.second;
}

bool ElementBuilder::should_inject_nodal_load_dof(const std::string &dof_variable_name,
                                                  dof_id_type dof_id) const
{
    const auto var_it = nodal_load_injector_elem_by_dof_by_variable.find(dof_variable_name);
    if (var_it == nodal_load_injector_elem_by_dof_by_variable.end())
        return false;

    const auto &injector_elem_by_dof = var_it->second;
    const auto dof_it = injector_elem_by_dof.find(dof_id);
    return (dof_it != injector_elem_by_dof.end()) &&
           (dof_it->second == static_cast<dof_id_type>(elem_data.elem->id()));
}

void ElementBuilder::write_element_rhs()
{
    const libMesh::MeshBase &mesh = equation_systems_->get_mesh();
    const auto &comm = mesh.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;

    const std::vector<std::string> out_filenames = {"rhs_" + rank_str + ".emat"};
    auto streams_out = reader.open_output_streams(emat_folderpath, out_filenames);
    auto &rhs_output_file = *streams_out[0];

    apec::logging::info("Calculating elements RHS");

    this->clear_nodal_load_injector_map();
    const auto dof_nodal_load_map =
        VariablesMap::get_dof_variable_nodal_load_map_for_system(*dof_solution_system_);
    for (const auto &entry : dof_nodal_load_map)
        this->build_nodal_load_injector_map(mesh, entry.first, entry.second);

    this->run_for_elements(mesh, [&](const libMesh::Elem *)
                           {
                               this->build_element_rhs();
                               BinaryReader::write_element_vectors(*this, rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS); });

    reader.close_file_streams(streams_out);
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
    elem_data.element_type = element_types_factory_->get_element_type(elem->subdomain_id());

    const vtkIdType elem_id = static_cast<vtkIdType>(elem->id());
    elem_data.material_num = static_cast<unsigned int>(material_cell_array_->GetComponent(elem_id, 0));
    elem_data.real_constant_num = static_cast<unsigned int>(real_constant_cell_array_->GetComponent(elem_id, 0));
    elem_data.coordinate_system_num = static_cast<unsigned int>(coordinate_system_cell_array_->GetComponent(elem_id, 0));
    elem_data.elem = elem;
    elem_data.coordinate_system = nullptr;

    if (elem_data.coordinate_system_num != 0)
    {
        const auto it = coordinate_systems_->find(static_cast<dof_id_type>(elem_data.coordinate_system_num));
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

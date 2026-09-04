/**
 * @file apec_m_system.cpp
 * @brief ApecMSystemBase and concrete system implementations.
 */
#include "apec_m_solver/apec_m_system.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "apec_m_solver/apec_m_equation_systems.h"
#include "apec_m_solver/apec_m_model.h"
#include "apec_m_solver/apec_m_logging.h"
#include "apec_m_solver/binary_reader.h"
#include "apec_m_solver/element_builder.h"
#include "apec_m_solver/element_types_factory.h"
#include "apec_m_solver/global_constraints_object.h"
#include "apec_m_solver/system_matrix_reduction.h"
#include "apec_m_solver/variables_map.h"

#include "libmesh/libmesh.h"
#include "libmesh/elem_range.h"
#include "libmesh/numeric_vector.h"

#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>

namespace
{
    std::filesystem::path get_emat_folderpath_from_solver()
    {
        return std::filesystem::current_path() / "emat";
    }

    using ElemRecordIndex = std::unordered_map<unsigned int, std::streampos>;

    NodalLoadInjectorElemByVarMap build_nodal_load_injector_map_for_system(
        const libMesh::MeshBase &mesh,
        const libMesh::ImplicitSystem &dof_system,
        vtkPointData *point_data)
    {
        NodalLoadInjectorElemByVarMap nodal_load_map;
        if (!point_data)
            return nodal_load_map;

        const auto dof_nodal_load_map =
            VariablesMap::get_dof_variable_nodal_load_map_for_system(dof_system);

        for (const auto &entry : dof_nodal_load_map)
        {
            const std::string &dof_variable_name = entry.first;
            const std::string &nodal_load_array_name = entry.second;

            if (!dof_system.has_variable(dof_variable_name))
                continue;

            vtkDataArray *load_array = point_data->GetArray(nodal_load_array_name.c_str());
            if (!load_array)
                continue;

            const unsigned int dof_var_num = dof_system.variable_number(dof_variable_name);
            const unsigned int sys_num = dof_system.number();
            const auto &dof_var = dof_system.variable(dof_var_num);

            std::map<libMesh::dof_id_type, libMesh::dof_id_type> local_injector_elem_by_dof;
            for (const auto &elem : mesh.active_local_element_ptr_range())
            {
                if (!dof_var.active_on_subdomain(elem->subdomain_id()))
                    continue;

                for (unsigned int i = 0; i < elem->n_nodes(); ++i)
                {
                    const auto &node = elem->node_ref(i);
                    const auto node_id = static_cast<vtkIdType>(node.id());
                    const libMesh::Real load_value = load_array->GetComponent(node_id, 0);
                    if (load_value == 0.0)
                        continue;

                    const libMesh::dof_id_type dof_id = node.dof_number(sys_num, dof_var_num, 0);
                    local_injector_elem_by_dof.emplace(dof_id, static_cast<libMesh::dof_id_type>(elem->id()));
                }
            }

            auto global_injector_elem_by_dof = local_injector_elem_by_dof;
            mesh.comm().set_union(global_injector_elem_by_dof);
            auto &injector_elem_by_dof = nodal_load_map[dof_variable_name];
            for (const auto &injector_entry : global_injector_elem_by_dof)
                injector_elem_by_dof[injector_entry.first] = injector_entry.second;
        }

        return nodal_load_map;
    }

}

const std::map<ApecMSystemBase::SystemType, std::string> ApecMSystemBase::SYSTEM_INFO_MAP = {
    {ApecMSystemBase::DOF_SOLUTION_SYSTEM, "DOF Solution"},
    {ApecMSystemBase::NODAL_SOLUTION, "Nodal Solution"},
    {ApecMSystemBase::ELEMENT_SOLUTION, "Element Solution"},
};

const std::string &ApecMSystemBase::system_name_from_type(ApecMSystemBase::SystemType system_type)
{
    auto it = ApecMSystemBase::SYSTEM_INFO_MAP.find(system_type);
    if (it == ApecMSystemBase::SYSTEM_INFO_MAP.end())
        apec::logging::error_and_abort("System type not found.");
    return it->second;
}

ApecMSystemBase::SystemType ApecMSystemBase::system_type_from_name(const std::string &system_name)
{
    for (const auto &pair : ApecMSystemBase::SYSTEM_INFO_MAP)
        if (pair.second == system_name)
            return pair.first;

    apec::logging::error_and_abort("System type not found for system name: " + system_name);
    return ApecMSystemBase::NUM_SYSTEMS;
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::prepare_time_step_solve(const pybind11::dict &time_step_dict)
{
    this->typed_equation_systems().prepare_time_step_solve(time_step_dict);
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::init_time_step(const pybind11::dict &time_step_dict)
{
    this->prepare_time_step_solve(time_step_dict);

    const auto *time_step = this->current_time_step();
    if (!time_step)
        apec::logging::error_and_abort("ApecMDOFSystem::init_time_step: time_step is null after prepare_time_step_solve.");

    if (time_step->updated)
    {
        *this->old_local_solution = *this->current_local_solution;
        this->old_local_solution->close();
        apec::logging::info("Advanced to new time step");
        apec::logging::table_root(apec::logging::format_time_step_table(*time_step));
    }

    this->get_dof_map().swap_dof_constraints();
    this->reinit_constraints();
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::compute_solution_time_derivative()
{
    const auto *time_step = this->current_time_step();
    if (!time_step)
        apec::logging::error_and_abort("ApecMDOFSystem::compute_solution_time_derivative: time_step is null.");

    const bool transient_effect = time_step->transient_effects;
    const libMesh::Real dt = time_step->dt;
    const libMesh::Real theta = time_step->theta;
    libMesh::NumericVector<libMesh::Number> *time_derivative_vector = this->request_vector("solution_time_derivative");
    if (!transient_effect)
    {
        this->older_local_solution->zero();
        time_derivative_vector->zero();
    }
    else
    {
        libMesh::Real theta_x_dt = theta * dt;
        libMesh::Real theta_minus_1_dt = (theta - 1) * dt;

        this->older_local_solution->scale(theta_minus_1_dt);
        *this->older_local_solution += *this->current_local_solution;
        *this->older_local_solution -= *this->old_local_solution;
        this->older_local_solution->scale(1 / theta_x_dt);

        time_derivative_vector->zero();
        *time_derivative_vector += *this->current_local_solution;
        *time_derivative_vector -= *this->old_local_solution;
        time_derivative_vector->scale(1 / dt);
    }
    this->older_local_solution->close();
    time_derivative_vector->close();
}

template <typename ImplicitSys>
bool ApecMDOFSystem<ImplicitSys>::check_update_system_matrix()
{
    auto bc_object = this->get_global_constraints_object();

    if (bc_object && bc_object->check_update_boundary_conditions())
    {
        apec::logging::info("New constraints: updating system matrix");
        return true;
    }

    bool update_matrix = false;
    const auto *current_step = this->current_time_step();
    const auto *old_step = this->previous_time_step();

    if (!current_step)
        update_matrix = true;
    else if (current_step->updated)
    {
        if (!old_step)
            update_matrix = true;
        else
        {
            const bool dt_differs = std::abs(current_step->dt - old_step->dt) > std::numeric_limits<libMesh::Real>::epsilon();
            const bool theta_differs = std::abs(current_step->theta - old_step->theta) > std::numeric_limits<libMesh::Real>::epsilon();
            const bool transient_effect_differs = current_step->transient_effects != old_step->transient_effects;
            if (dt_differs || theta_differs || transient_effect_differs)
                update_matrix = true;
        }
    }

    if (update_matrix)
        apec::logging::info("Updating system matrix for current time-step");
    else
        apec::logging::info("Reusing system matrix from previous time step");

    return update_matrix;
}

bool ApecMNonLinearSystem::check_update_system_matrix()
{
    if (this->diff_solver->total_inner_iterations() > 0)
    {
        apec::logging::info("Reusing system matrix from previous time step");
        return false;
    }
    else
        return ApecMDOFSystem<libMesh::NonlinearImplicitSystem>::check_update_system_matrix();
}

template <typename ImplicitSys>
bool ApecMDOFSystem<ImplicitSys>::check_update_element_matrices()
{
    return false;
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::write_element_data()
{
    this->write_element_matrices();
    this->write_element_rhs();
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::write_element_matrices()
{
    auto &equation_systems = this->typed_equation_systems();
    const libMesh::MeshBase &mesh = equation_systems.get_mesh();
    const auto &comm = mesh.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;

    const std::vector<std::string> filenames = {
        "mtrx_" + rank_str + ".emat",
        "nodal_sol_mtrx_" + rank_str + ".emat",
        "element_sol_mtrx_" + rank_str + ".emat"};

    const bool emat_files_exist = std::all_of(
        filenames.begin(),
        filenames.end(),
        [&](const std::string &filename)
        {
            return std::filesystem::exists(emat_folderpath / filename);
        });

    if (emat_files_exist)
    {
        apec::logging::info("Reusing element matrices from previous time step");
        return;
    }

    apec::logging::info("Computing element matrices for current time-step");

    auto streams = reader.open_output_streams(emat_folderpath, filenames);
    auto &mtrx_output_file = *streams[0];
    auto &nodal_sol_mtrx_output_file = *streams[1];
    auto &element_sol_mtrx_output_file = *streams[2];

    ElementBuilder builder(&equation_systems);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        builder.build_element();

        const unsigned int elem_id = static_cast<unsigned int>(elem->id());
        BinaryReader::write_element_matrices(
            mtrx_output_file,
            elem_id,
            builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
        BinaryReader::write_element_matrices(
            nodal_sol_mtrx_output_file,
            elem_id,
            builder.get_matrices()[ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES]);
        BinaryReader::write_element_matrices(
            element_sol_mtrx_output_file,
            elem_id,
            builder.get_matrices()[ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES]);
    }

    reader.close_file_streams(streams);
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::write_element_rhs()
{
    auto &equation_systems = this->typed_equation_systems();
    const libMesh::MeshBase &mesh = equation_systems.get_mesh();
    const auto &comm = mesh.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;

    const std::vector<std::string> out_filenames = {"rhs_" + rank_str + ".emat"};
    auto streams_out = reader.open_output_streams(emat_folderpath, out_filenames);
    auto &rhs_output_file = *streams_out[0];

    apec::logging::info("Calculating elements RHS");

    const auto *nodal_load_map = this->ensure_nodal_load_injector_map();

    ElementBuilder builder(&equation_systems);
    builder.set_nodal_load_map(nodal_load_map);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        builder.build_element_rhs();

        const unsigned int elem_id = static_cast<unsigned int>(elem->id());
        BinaryReader::write_element_vectors(
            rhs_output_file,
            elem_id,
            builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS]);
    }

    reader.close_file_streams(streams_out);
}

template <typename ImplicitSys>
const NodalLoadInjectorElemByVarMap *ApecMDOFSystem<ImplicitSys>::ensure_nodal_load_injector_map()
{
    if (!nodal_load_injector_map_)
    {
        auto &equation_systems = this->typed_equation_systems();
        auto *solver_parent = equation_systems.get_parent_object();
        if (!solver_parent)
            apec::logging::error_and_abort("ApecMDOFSystem::ensure_nodal_load_injector_map: solver parent not available.");

        auto *vtk_grid = solver_parent->get_active_vtk_grid();
        if (!vtk_grid)
            apec::logging::error_and_abort("ApecMDOFSystem::ensure_nodal_load_injector_map: active VTK grid not available.");

        nodal_load_injector_map_ = std::make_shared<NodalLoadInjectorElemByVarMap>(
            build_nodal_load_injector_map_for_system(
                equation_systems.get_mesh(),
                *this,
                vtk_grid->GetPointData()));
    }

    return nodal_load_injector_map_.get();
}

template <typename ImplicitSys>
void ApecMDOFSystem<ImplicitSys>::post_init()
{
    apec::logging::info("Initializing system matrix");
    this->zero_out_matrix_and_rhs = false;
    this->get_system_matrix().zero();
    apec::logging::info("Setting linear solver configurations");
}

void ApecMLinearSystem::release_subset_solve()
{
    if (this->_subset != nullptr)
    {
        subset_solve_holder_.reset();
        this->_subset = nullptr;
    }
}

void ApecMLinearSystem::delete_rows_cols_constrained_dofs()
{
    apec::logging::info("Removing constrained DOFs from system matrix and RHS");
    this->release_subset_solve();
    subset_solve_holder_ = std::make_shared<SystemMatrixReduction>(*this);
    this->restrict_solve_to(subset_solve_holder_.get(), libMesh::SUBSET_ZERO);
}

void ApecMNonLinearSystem::delete_rows_cols_constrained_dofs()
{
    apec::logging::info("Removing constrained DOFs from system matrix and RHS");
    subset_solve_holder_ = std::make_shared<SystemMatrixReduction>(*this);
    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    auto *newton_petsc_solver =
        dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(&newton_solver->get_linear_solver());
    newton_petsc_solver->restrict_solve_to(&subset_solve_holder_.get()->dof_ids(), libMesh::SUBSET_ZERO);
}

void ApecMLinearSystem::solve()
{
    libMesh::TransientSystem<libMesh::LinearImplicitSystem>::solve();
    this->get_dof_map().enforce_constraints_exactly(*this);
    this->update();
    this->compute_solution_time_derivative();
    apec::logging::info("Finished linear system solve");
}

void ApecMLinearSystem::post_init()
{
    ApecMDOFSystem<libMesh::LinearImplicitSystem>::post_init();

    auto &equation_systems = this->get_apec_equation_systems();
    auto *solver_parent = equation_systems.get_parent_object();
    auto *petsc_solver = this->get_petsc_solver();
    if (solver_parent && petsc_solver)
    {
        this->solver_config_copy_ = std::make_unique<PetscSolverConfiguration>(
            solver_parent->get_apec_object().attr("solver_settings"),
            petsc_solver);
        this->solver_config_copy_->configure_solver();
    }
    else
        apec::logging::error_and_abort("Solver parent or PETSc linear solver not available.");
}

void ApecMLinearSystem::assemble_system()
{
    libMesh::PerfLog perf_log("Matrix Assembly");

    auto &system = *this;
    libMesh::SparseMatrix<libMesh::Number> &matrix = system.get_system_matrix();
    system.rhs->zero();

    const libMesh::DofMap &dof_map = system.get_dof_map();
    const libMesh::MeshBase &mesh = this->typed_equation_systems().get_mesh();
    auto *solver = this->get_petsc_solver();
    const auto &comm = system.comm();

    const bool update_matrix = this->check_update_system_matrix();
    if (update_matrix)
    {
        matrix.zero();
        solver->reuse_preconditioner(false);
    }
    else
        solver->reuse_preconditioner(true);

    const auto *time_step = this->current_time_step();
    const libMesh::Real dt = time_step->dt;
    const libMesh::Real theta = time_step->theta;
    const bool transient_effect = time_step->transient_effects;
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();
    const std::string rank_str = std::to_string(comm.rank());

    BinaryReader reader;
    const std::filesystem::path mtrx_filepath = emat_folderpath / ("mtrx_" + rank_str + ".emat");
    const std::filesystem::path rhs_filepath = emat_folderpath / ("rhs_" + rank_str + ".emat");

    auto streams = reader.open_input_streams(emat_folderpath,
                                             {"mtrx_" + rank_str + ".emat",
                                              "rhs_" + rank_str + ".emat"});
    auto &mtrx_input_file = *streams[0];
    auto &rhs_input_file = *streams[1];

    apec::logging::info("Building global matrix and RHS");

    ElementBuilder builder(&this->typed_equation_systems());
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::Real> old_solution_values;
    std::vector<libMesh::Real> old_time_derivative_values;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        dof_map.dof_indices(elem, dof_indices);

        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(mtrx_input_file,
                                            elem_id,
                                            builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
        BinaryReader::read_element_vectors(rhs_input_file,
                                           elem_id,
                                           builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS]);

        auto &Ke = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][0];
        auto &Ce = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][1];

        if (transient_effect)
        {
            Ce.scale(1 / (dt * theta));
            Ke += Ce;
        }

        auto &rhs = builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];

        if (transient_effect)
        {
            system.old_local_solution->get(dof_indices, old_solution_values);
            Ce.vector_mult_add(rhs, 1.0, old_solution_values);

            if (theta != 1.0)
            {
                system.older_local_solution->get(dof_indices, old_time_derivative_values);
                Ce.vector_mult_add(rhs, dt * (1.0 - theta), old_time_derivative_values);
            }
        }

        if (update_matrix)
        {
            dof_map.heterogeneously_constrain_element_matrix_and_vector(Ke, rhs, dof_indices);
            matrix.add_matrix(Ke, dof_indices);
        }
        else
            dof_map.heterogeneously_constrain_element_vector(Ke, rhs, dof_indices);

        system.rhs->add_vector(rhs, dof_indices);
    }

    reader.close_file_streams(streams);
    comm.barrier();
    apec::logging::info("Linear solver output info");
}

void ApecMNonLinearSystem::assemble_system_impl(bool assemble_jacobian, bool assemble_residual)
{
    libMesh::PerfLog perf_log("Matrix Assembly");

    auto &system = *this;
    libMesh::SparseMatrix<libMesh::Number> &matrix = system.get_system_matrix();
    libMesh::SparseMatrix<libMesh::Number> &linear_matrix = system.get_matrix("linear_stiffness");

    const libMesh::DofMap &dof_map = system.get_dof_map();
    const libMesh::MeshBase &mesh = this->typed_equation_systems().get_mesh();
    const auto &comm = system.comm();

    const bool update_matrix = this->check_update_system_matrix();
    const auto *time_step = this->current_time_step();
    const libMesh::Real dt = time_step->dt;
    const libMesh::Real theta = time_step->theta;
    const bool transient_effect = time_step->transient_effects;

    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();
    const std::string rank_str = std::to_string(comm.rank());

    if (assemble_jacobian)
    {
        if (update_matrix)
            linear_matrix.zero();
        matrix.zero();
    }
    if (assemble_residual)
        system.rhs->zero();

    BinaryReader reader;
    const std::filesystem::path mtrx_filepath = emat_folderpath / ("mtrx_" + rank_str + ".emat");
    const std::filesystem::path rhs_filepath = emat_folderpath / ("rhs_" + rank_str + ".emat");

    auto streams = reader.open_input_streams(emat_folderpath,
                                             {"mtrx_" + rank_str + ".emat",
                                              "rhs_" + rank_str + ".emat"});
    auto &mtrx_input_file = *streams[0];
    auto &rhs_input_file = *streams[1];

    if (assemble_jacobian && assemble_residual)
        apec::logging::info("Building global Jacobian and residual");
    else if (assemble_jacobian)
        apec::logging::info("Building global Jacobian");
    else
        apec::logging::info("Building global residual");

    ElementBuilder builder(&this->typed_equation_systems());
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::Real> current_solution_values;
    std::vector<libMesh::Real> old_solution_values;
    std::vector<libMesh::Real> old_time_derivative_values;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        dof_map.dof_indices(elem, dof_indices);

        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(mtrx_input_file,
                                            elem_id,
                                            builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
        BinaryReader::read_element_vectors(rhs_input_file,
                                           elem_id,
                                           builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS]);

        auto &linear_ke_base = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][0];
        auto &ce = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][1];
        auto &linear_rhs = builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];

        auto linear_ke = linear_ke_base;
        if (transient_effect)
        {
            ce.scale(1 / (dt * theta));
            linear_ke += ce;
        }

        system.current_local_solution->get(dof_indices, current_solution_values);
        linear_ke.vector_mult_add(linear_rhs, -1.0, current_solution_values);

        if (transient_effect)
        {
            system.old_local_solution->get(dof_indices, old_solution_values);
            ce.vector_mult_add(linear_rhs, 1.0, old_solution_values);

            if (theta != 1.0)
            {
                system.older_local_solution->get(dof_indices, old_time_derivative_values);
                ce.vector_mult_add(linear_rhs, dt * (1.0 - theta), old_time_derivative_values);
            }
        }

        // NewtonSolver solves J * du = r and then updates u <- u - du, so r
        // must be the true residual F(u), not -F(u).
        linear_rhs.scale(-1.0);

        const bool nonlinear_element = builder.nonlinear_element();
        if (nonlinear_element)
        {
            builder.build_element_tangent_stiffness_residual();
            auto &tangent_ke = builder.get_matrices()[ElementBuilder::MatricesSet::TANGENT_STIFFNESS_MATRICES][0];
            auto &residual = builder.get_vectors()[ElementBuilder::VectorsSet::RESIDUAL_VECTORS][0];

            if (assemble_jacobian)
                linear_ke += tangent_ke;
            if (assemble_residual)
                linear_rhs += residual;

            if (assemble_jacobian)
            {
                if (assemble_residual)
                    dof_map.heterogeneously_constrain_element_jacobian_and_residual(
                        linear_ke,
                        linear_rhs,
                        dof_indices,
                        *system.current_local_solution);
                else
                    dof_map.heterogeneously_constrain_element_matrix_and_vector(
                        linear_ke,
                        linear_rhs,
                        dof_indices);

                matrix.add_matrix(linear_ke, dof_indices);
            }

            if (assemble_residual)
            {
                if (!assemble_jacobian)
                    dof_map.heterogeneously_constrain_element_residual(
                        linear_rhs,
                        dof_indices,
                        *system.current_local_solution);

                system.rhs->add_vector(linear_rhs, dof_indices);
            }
        }
        else
        {
            if (assemble_jacobian && update_matrix)
            {
                if (assemble_residual)
                    dof_map.heterogeneously_constrain_element_jacobian_and_residual(
                        linear_ke,
                        linear_rhs,
                        dof_indices,
                        *system.current_local_solution);
                else
                    dof_map.heterogeneously_constrain_element_matrix_and_vector(
                        linear_ke,
                        linear_rhs,
                        dof_indices);
                linear_matrix.add_matrix(linear_ke, dof_indices);
            }
            else if (assemble_residual)
                dof_map.heterogeneously_constrain_element_residual(
                    linear_rhs,
                    dof_indices,
                    *system.current_local_solution);

            if (assemble_residual)
                system.rhs->add_vector(linear_rhs, dof_indices);
        }
    }

    if (assemble_jacobian)
    {
        // MatAXPY (add()) requires both operands assembled; add_matrix() above leaves them open.
        matrix.close();
        linear_matrix.close();
        matrix.add(1.0, linear_matrix);
    }

    reader.close_file_streams(streams);
    comm.barrier();
    if (assemble_jacobian && assemble_residual)
        apec::logging::info("Newton-Rhapson solver output info");
}

void ApecMNonLinearSystem::assemble_system()
{
    this->assemble_system_impl(true, true);
}

void ApecMNonLinearSystem::assemble_jacobian_only()
{
    this->assemble_system_impl(true, false);
}

void ApecMNonLinearSystem::assemble_residual_only()
{
    this->assemble_system_impl(false, true);
}

void ApecMNonLinearSystem::solve()
{
    libMesh::NonlinearImplicitSystem::solve();
    this->compute_solution_time_derivative();
    apec::logging::info("Finished nonlinear system solve");
}

void ApecMNonLinearSystem::post_init()
{
    ApecMDOFSystem<libMesh::NonlinearImplicitSystem>::post_init();

    auto &equation_systems = this->get_apec_equation_systems();
    auto *solver_parent = equation_systems.get_parent_object();
    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    if (!solver_parent || !newton_solver)
        apec::logging::error_and_abort("ApecMNonLinearSystem::post_init: solver parent or Newton solver not available.");

    auto *newton_petsc_solver =
        dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(&newton_solver->get_linear_solver());
    if (!newton_petsc_solver)
        apec::logging::error_and_abort("ApecMNonLinearSystem::post_init: Newton linear solver is not PETSc.");

    this->nonlinear_solver_config_copy_ = std::make_unique<PetscSolverConfiguration>(
        solver_parent->get_apec_object().attr("solver_settings"),
        newton_petsc_solver);
    this->nonlinear_solver_config_copy_->configure_solver();
}

void ApecMNodalSolutionSystem::assemble_system()
{
    apec::logging::info("Computing nodal solution");

    auto &eq_sys = this->typed_equation_systems();

    const libMesh::MeshBase &mesh = eq_sys.get_mesh();
    const auto &comm = mesh.comm();
    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"nodal_sol_mtrx_" + rank_str + ".emat"};
    const std::filesystem::path nodal_mtrx_filepath = emat_folderpath / inp_filenames[0];

    auto inp_streams = reader.open_input_streams(emat_folderpath, inp_filenames);
    auto &nodal_solution_mtrx_input_file = *inp_streams[0];

    ElementBuilder builder(&eq_sys);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(
            nodal_solution_mtrx_input_file,
            elem_id,
            builder.get_matrices()[ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES]);

        builder.set_element_data(elem);
        builder.build_nodal_solution();
    }

    reader.close_file_streams(inp_streams);

    this->solution->close();
    apec::logging::info("Nodal solution updated");
}

void ApecMElementSolutionSystem::assemble_system()
{
    apec::logging::info("Computing elemental solution");

    auto &eq_sys = this->typed_equation_systems();

    const libMesh::MeshBase &mesh = eq_sys.get_mesh();
    const auto &comm = mesh.comm();
    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"element_sol_mtrx_" + rank_str + ".emat"};
    const std::filesystem::path element_mtrx_filepath = emat_folderpath / inp_filenames[0];

    auto inp_streams = reader.open_input_streams(emat_folderpath, inp_filenames);
    auto &element_solution_mtrx_input_file = *inp_streams[0];

    ElementBuilder builder(&eq_sys);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(
            element_solution_mtrx_input_file,
            elem_id,
            builder.get_matrices()[ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES]);

        builder.set_element_data(elem);
        builder.build_element_solution();
    }

    reader.close_file_streams(inp_streams);

    this->solution->close();
    apec::logging::info("Element solution updated");
}

template <typename SysType>
void ApecMSystem<SysType>::init_variables()
{
    auto &equation_systems = this->typed_equation_systems();
    auto element_types_collector = equation_systems.get_element_types_factory();
    const auto type = ApecMSystemBase::system_type_from_name(this->name());

    if (!element_types_collector || !element_types_collector->get_element_types_map())
        apec::logging::error_and_abort("Element types map not available for system variables.");

    auto vars_to_subdomain_ids = VariablesMap::associate_variables_with_element_IDs(
        *element_types_collector->get_element_types_map(),
        type);

    std::set<libMesh::subdomain_id_type> mesh_subdomain_ids;
    this->get_mesh().subdomain_ids(mesh_subdomain_ids);

    for (auto &pair : vars_to_subdomain_ids)
    {
        const auto &variable_string = VariablesMap::get_info_by_key(pair.first).global_label;
        auto &subdomain_ids = pair.second;
        for (auto it = subdomain_ids.begin(); it != subdomain_ids.end();)
        {
            if (!mesh_subdomain_ids.count(*it))
                it = subdomain_ids.erase(it);
            else
                ++it;
        }
        if (subdomain_ids.empty())
            continue;

        auto elem_type = element_types_collector->get_element_type(*subdomain_ids.begin());
        auto fe_family = VariablesMap::get_variable_fe_family(*elem_type, pair.first);
        auto order = VariablesMap::get_variable_order(*elem_type, pair.first);

        this->add_variable(
            variable_string,
            order,
            fe_family,
            &subdomain_ids);
    }
}

template class ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class ApecMSystem<libMesh::TransientSystem<libMesh::NonlinearImplicitSystem>>;
template class ApecMSystem<libMesh::ExplicitSystem>;
template class ApecMDOFSystem<libMesh::LinearImplicitSystem>;
template class ApecMDOFSystem<libMesh::NonlinearImplicitSystem>;

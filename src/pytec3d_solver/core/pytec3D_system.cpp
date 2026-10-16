/**
 * @file pytec3d_system.cpp
 * @brief pyTEC3DSystemBase and concrete system implementations.
 */
#include "pytec3d_system.h"

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

#include "pytec3d_equation_systems.h"
#include "pytec3d_model.h"
#include "logging.h"
#include "binary_reader.h"
#include "element_builder.h"
#include "element_types_factory.h"
#include "constraints.h"
#include "system_matrix_reduction.h"
#include "variables_map.h"

#include <libmesh/libmesh.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/elem_range.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/petsc_matrix.h> // for libMesh::PetscMatrix
#include <libmesh/petsc_macro.h>
#include <libmesh/petsc_vector.h>  // for libMesh::PetscVector
#include <libmesh/wrapped_petsc.h> // for WrappedPetsc

#include <petscmat.h> // for MatResetHash

#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkUnstructuredGrid.h>

using libMesh::libmesh_ignore;

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

namespace
{
    /**
     * @brief The time scheme on one element: Keff = K + c_v C + c_a M and rhs -= C h_v + M h_a.
     *
     * @p matrices are the element K, C and, for second-order element types, M; h_v and h_a come from the
     * state at the start of the step (u_n, v_n, a_n) on the element DOFs.
     */
    libMesh::DenseMatrix<libMesh::Number> apply_time_scheme(const std::vector<libMesh::DenseMatrix<libMesh::Number>> &matrices,
                                                            const TimeScheme &scheme,
                                                            const libMesh::System &system,
                                                            const libMesh::NumericVector<libMesh::Number> &old_solution,
                                                            const std::vector<libMesh::dof_id_type> &dof_indices,
                                                            libMesh::DenseVector<libMesh::Number> &rhs)
    {
        libMesh::DenseMatrix<libMesh::Number> keff = matrices[0];
        if (!scheme.transient())
            return keff;

        std::vector<libMesh::Number> u_n, v_n, a_n;
        old_solution.get(dof_indices, u_n);
        system.get_vector(pyTEC3DSystemBase::kOldVelocity).get(dof_indices, v_n);
        system.get_vector(pyTEC3DSystemBase::kOldAcceleration).get(dof_indices, a_n);

        const auto add_term = [&](const libMesh::DenseMatrix<libMesh::Number> &matrix, libMesh::Real c,
                                  libMesh::Real hu, libMesh::Real hv, libMesh::Real ha)
        {
            if (matrix.m() == 0 || (c == 0.0 && hu == 0.0 && hv == 0.0 && ha == 0.0))
                return;
            keff.add(c, matrix);
            libMesh::DenseVector<libMesh::Number> history(static_cast<unsigned int>(dof_indices.size()));
            for (unsigned int i = 0; i < history.size(); ++i)
                history(i) = hu * u_n[i] + hv * v_n[i] + ha * a_n[i];
            matrix.vector_mult_add(rhs, -1.0, history);
        };
        add_term(matrices[1], scheme.c_v, scheme.hu_v, scheme.hv_v, scheme.ha_v);
        if (matrices.size() > 2)
            add_term(matrices[2], scheme.c_a, scheme.hu_a, scheme.hv_a, scheme.ha_a);
        return keff;
    }
} // namespace

const std::map<pyTEC3DSystemBase::SystemType, std::string> pyTEC3DSystemBase::SYSTEM_INFO_MAP = {
    {pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM, "DOF Solution"},
    {pyTEC3DSystemBase::NODAL_SOLUTION, "Nodal Solution"},
    {pyTEC3DSystemBase::ELEMENT_SOLUTION, "Element Solution"},
};

const std::string &pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::SystemType system_type)
{
    auto it = pyTEC3DSystemBase::SYSTEM_INFO_MAP.find(system_type);
    if (it == pyTEC3DSystemBase::SYSTEM_INFO_MAP.end())
        pytec3d::logging::error_and_abort("System type not found.");
    return it->second;
}

pyTEC3DSystemBase::SystemType pyTEC3DSystemBase::system_type_from_name(const std::string &system_name)
{
    for (const auto &pair : pyTEC3DSystemBase::SYSTEM_INFO_MAP)
        if (pair.second == system_name)
            return pair.first;

    pytec3d::logging::error_and_abort("System type not found for system name: " + system_name);
    return pyTEC3DSystemBase::NUM_SYSTEMS;
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::prepare_time_step_solve(const pybind11::dict &time_step_dict)
{
    this->typed_equation_systems().prepare_time_step_solve(time_step_dict);
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::init_time_step(const pybind11::dict &time_step_dict)
{
    this->prepare_time_step_solve(time_step_dict);

    const auto *time_step = this->current_time_step();
    if (!time_step)
        pytec3d::logging::error_and_abort("pyTEC3DDOFSystem::init_time_step: time_step is null after prepare_time_step_solve.");

    if (time_step->updated)
    {
        // The end of the previous step is the start of this one.
        *this->old_local_solution = *this->current_local_solution;
        this->old_local_solution->close();
        this->get_vector(pyTEC3DSystemBase::kOldVelocity) = this->get_vector(pyTEC3DSystemBase::kVelocity);
        this->get_vector(pyTEC3DSystemBase::kOldVelocity).close();
        this->get_vector(pyTEC3DSystemBase::kOldAcceleration) = this->get_vector(pyTEC3DSystemBase::kAcceleration);
        this->get_vector(pyTEC3DSystemBase::kOldAcceleration).close();
        pytec3d::logging::info("Advanced to new time step");
        pytec3d::logging::table_root(pytec3d::logging::format_time_step_table(*time_step));
    }

    this->get_dof_map().swap_dof_constraints();
    this->reinit_constraints();
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::update_time_derivatives()
{
    // v = c_v u + hu_v u_n + hv_v v_n + ha_v a_n, a likewise (time_step.h); zero on a static step.
    const TimeScheme scheme = this->typed_equation_systems().time_scheme();
    const auto &v_n = this->get_vector(pyTEC3DSystemBase::kOldVelocity);
    const auto &a_n = this->get_vector(pyTEC3DSystemBase::kOldAcceleration);
    const auto update = [&](libMesh::NumericVector<libMesh::Number> &x, libMesh::Real c, libMesh::Real hu, libMesh::Real hv,
                            libMesh::Real ha)
    {
        x.zero();
        if (scheme.transient())
        {
            x.add(c, *this->current_local_solution);
            x.add(hu, *this->old_local_solution);
            x.add(hv, v_n);
            x.add(ha, a_n);
        }
        x.close();
    };
    update(this->get_vector(pyTEC3DSystemBase::kVelocity), scheme.c_v, scheme.hu_v, scheme.hv_v, scheme.ha_v);
    update(this->get_vector(pyTEC3DSystemBase::kAcceleration), scheme.c_a, scheme.hu_a, scheme.hv_a, scheme.ha_a);
    const libMesh::Real rate = scheme.transient() ? 1.0 / this->current_time_step()->dt : 0.0;
    update(this->get_vector(pyTEC3DSystemBase::kMeanVelocity), rate, -rate, 0.0, 0.0);
}

template <typename ImplicitSys>
bool pyTEC3DDOFSystem<ImplicitSys>::check_update_system_matrix()
{
    auto bc_object = this->get_global_constraints_object();

    if (bc_object && bc_object->check_update_boundary_conditions())
    {
        pytec3d::logging::info("New constraints: updating system matrix");
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
            const auto differs = [](libMesh::Real a, libMesh::Real b)
            { return std::abs(a - b) > std::numeric_limits<libMesh::Real>::epsilon(); };
            const bool scheme_differs = differs(current_step->theta, old_step->theta) ||
                                        differs(current_step->beta, old_step->beta) ||
                                        differs(current_step->gamma, old_step->gamma);
            const bool transient_effect_differs = current_step->transient_effects != old_step->transient_effects;
            if (dt_differs || scheme_differs || transient_effect_differs)
                update_matrix = true;
        }
    }

    if (update_matrix)
        pytec3d::logging::info("Updating system matrix for current time-step");
    else
        pytec3d::logging::info("Reusing system matrix from previous time step");

    return update_matrix;
}

bool pyTEC3DNonLinearSystem::check_update_system_matrix()
{
    if (this->diff_solver->total_inner_iterations() > 0)
    {
        pytec3d::logging::info("Reusing system matrix from previous time step");
        return false;
    }
    else
        return pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>::check_update_system_matrix();
}

template <typename ImplicitSys>
bool pyTEC3DDOFSystem<ImplicitSys>::check_update_element_matrices()
{
    return false;
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::write_element_matrices()
{
    auto &equation_systems = this->typed_equation_systems();
    const libMesh::MeshBase &mesh = equation_systems.get_mesh();
    const auto &comm = mesh.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;

    // K, C (M), and the operators of the static nodal and element fields, element by element in the order of the
    // local elements, read back in the same order by the DOF, nodal and element systems.
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
        pytec3d::logging::info("Reusing element matrices from previous time step");
        return;
    }

    pytec3d::logging::info("Computing element matrices for current time-step");

    auto streams = reader.open_output_streams(emat_folderpath, filenames);

    ElementBuilder builder(&equation_systems);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        builder.build_element();

        const auto elem_id = static_cast<unsigned int>(elem->id());
        BinaryReader::write_element_matrices(*streams[0], elem_id,
                                             builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
        BinaryReader::write_element_matrices(*streams[1], elem_id,
                                             builder.get_matrices()[ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES]);
        BinaryReader::write_element_matrices(*streams[2], elem_id,
                                             builder.get_matrices()[ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES]);
    }

    reader.close_file_streams(streams);
}

template <typename ImplicitSys>
const NodalLoadInjectorElemByVarMap *pyTEC3DDOFSystem<ImplicitSys>::ensure_nodal_load_injector_map()
{
    if (!nodal_load_injector_map_)
    {
        auto &equation_systems = this->typed_equation_systems();
        auto *solver_parent = equation_systems.get_parent_object();
        if (!solver_parent)
            pytec3d::logging::error_and_abort("pyTEC3DDOFSystem::ensure_nodal_load_injector_map: solver parent not available.");

        auto *vtk_grid = solver_parent->get_active_vtk_grid();
        if (!vtk_grid)
            pytec3d::logging::error_and_abort("pyTEC3DDOFSystem::ensure_nodal_load_injector_map: active VTK grid not available.");

        nodal_load_injector_map_ = std::make_shared<NodalLoadInjectorElemByVarMap>(
            build_nodal_load_injector_map_for_system(
                equation_systems.get_mesh(),
                *this,
                vtk_grid->GetPointData()));
    }

    return nodal_load_injector_map_.get();
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::post_init()
{
    pytec3d::logging::info("Initializing system matrix");
    this->zero_out_matrix_and_rhs = false;
    this->get_system_matrix().zero();
    pytec3d::logging::info("Setting linear solver configurations");
}

void pyTEC3DLinearSystem::release_subset_solve()
{
    if (this->_subset != nullptr)
    {
        subset_solve_holder_.reset();
        this->_subset = nullptr;
    }

    clear_reduced_solve_cache();
}

void pyTEC3DLinearSystem::clear_reduced_solve_cache()
{
    reduced_subset_dofs_cache_.clear();
    reduced_subset_is_.reset_to_zero();
    reduced_scatter_.reset_to_zero();
    reduced_submat_.reset_to_zero();
    reduced_subprecond_.reset_to_zero();
    reduced_subrhs_.reset_to_zero();
    reduced_subsolution_.reset_to_zero();
    reduced_parent_mat_cache_ = nullptr;
    reduced_parent_precond_mat_cache_ = nullptr;
    reduced_parent_rhs_vec_cache_ = nullptr;
    reduced_cache_ready_ = false;
}

void pyTEC3DLinearSystem::delete_rows_cols_constrained_dofs()
{
    pytec3d::logging::info("Removing constrained DOFs from system matrix and RHS");
    this->release_subset_solve();
    subset_solve_holder_ = std::make_shared<SystemMatrixReduction>(*this);
    this->restrict_solve_to(subset_solve_holder_.get(), libMesh::SUBSET_ZERO);
}

void pyTEC3DNonLinearSystem::delete_rows_cols_constrained_dofs()
{
    pytec3d::logging::info("Removing constrained DOFs from system matrix and RHS");
    subset_solve_holder_ = std::make_shared<SystemMatrixReduction>(*this);
    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    auto *newton_petsc_solver =
        dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(&newton_solver->get_linear_solver());
    newton_petsc_solver->restrict_solve_to(&subset_solve_holder_.get()->dof_ids(), libMesh::SUBSET_ZERO);
}

void pyTEC3DLinearSystem::solve()
{
    if (this->assemble_before_solve)
        // Assemble the linear system
        this->assemble();

    // Get the user-specified linear solver tolerance
    const auto [maxits, tol] = this->get_linear_solve_parameters();

    // Solve the linear system.  Several cases:
    std::pair<unsigned int, libMesh::Real> rval = std::make_pair(0, 0.0);

    auto *matrix = static_cast<libMesh::PetscMatrixBase<libMesh::Number> *>(this->matrix);
    auto *precond = static_cast<libMesh::PetscMatrixBase<libMesh::Number> *>(this->request_matrix("Preconditioner"));
    auto *solution = static_cast<libMesh::PetscVector<libMesh::Number> *>(this->solution.get());
    auto *rhs = static_cast<libMesh::PetscVector<libMesh::Number> *>(this->rhs);

    if (!precond)
        precond = matrix;

    matrix->close();

    if (precond != matrix)
        precond->close();
    solution->close();
    rhs->close();

    const auto &subset_dofs = subset_solve_holder_->dof_ids();
    if (subset_dofs.empty())
        pytec3d::logging::error_and_abort("pyTEC3DLinearSystem::solve: subset solve requested with empty subset dof list.");

    const Mat parent_mat = matrix->mat();
    const Mat parent_precond_mat = precond->mat();
    const bool allow_reuse = linear_solver->get_same_preconditioner();

    if (!allow_reuse)
    {
        clear_reduced_solve_cache();

        linear_solver->clear();
        if (this->prefix_with_name())
            linear_solver->init(this->prefix().c_str());
        else
            linear_solver->init();
        linear_solver->init_systems(*this);

        PetscInt *petsc_dofs = nullptr;
        LibmeshPetscCall(PetscMalloc(subset_dofs.size() * sizeof(PetscInt), &petsc_dofs));
        for (std::size_t i = 0; i < subset_dofs.size(); ++i)
            petsc_dofs[i] = static_cast<PetscInt>(subset_dofs[i]);

        LibmeshPetscCall(ISCreateGeneral(this->comm().get(),
                                         static_cast<PetscInt>(subset_dofs.size()),
                                         petsc_dofs,
                                         PETSC_OWN_POINTER,
                                         reduced_subset_is_.get()));

        PetscInt subset_local_size = 0;
        LibmeshPetscCall(ISGetLocalSize(reduced_subset_is_, &subset_local_size));

        LibmeshPetscCall(VecCreate(this->comm().get(), reduced_subrhs_.get()));
        LibmeshPetscCall(VecSetSizes(reduced_subrhs_, subset_local_size, PETSC_DECIDE));
        LibmeshPetscCall(VecSetFromOptions(reduced_subrhs_));

        LibmeshPetscCall(VecCreate(this->comm().get(), reduced_subsolution_.get()));
        LibmeshPetscCall(VecSetSizes(reduced_subsolution_, subset_local_size, PETSC_DECIDE));
        LibmeshPetscCall(VecSetFromOptions(reduced_subsolution_));

        LibmeshPetscCall(VecScatterCreate(rhs->vec(),
                                          reduced_subset_is_,
                                          reduced_subrhs_,
                                          nullptr,
                                          reduced_scatter_.get()));

        LibmeshPetscCall(LibMeshCreateSubMatrix(parent_mat,
                                                reduced_subset_is_,
                                                reduced_subset_is_,
                                                MAT_INITIAL_MATRIX,
                                                reduced_submat_.get()));

        LibmeshPetscCall(LibMeshCreateSubMatrix(parent_mat,
                                                reduced_subset_is_,
                                                reduced_subset_is_,
                                                MAT_INITIAL_MATRIX,
                                                reduced_subprecond_.get()));
        // reduced_subprecond_ = reduced_submat_;

        reduced_subset_dofs_cache_ = subset_dofs;
        reduced_parent_mat_cache_ = parent_mat;
        reduced_parent_precond_mat_cache_ = parent_precond_mat;
        reduced_parent_rhs_vec_cache_ = rhs->vec();
        reduced_cache_ready_ = true;
    }

    VecScatterBeginEnd(this->comm(), reduced_scatter_, rhs->vec(), reduced_subrhs_, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterBeginEnd(this->comm(), reduced_scatter_, solution->vec(), reduced_subsolution_, INSERT_VALUES, SCATTER_FORWARD);

    libMesh::PetscMatrix<libMesh::Number> submat_wrapper(reduced_submat_, this->comm());
    libMesh::PetscVector<libMesh::Number> subrhs_wrapper(reduced_subrhs_, this->comm());
    libMesh::PetscVector<libMesh::Number> subsolution_wrapper(reduced_subsolution_, this->comm());

    libMesh::PetscMatrix<libMesh::Number> subprecond_wrapper(reduced_subprecond_, this->comm());
    rval = linear_solver->solve(submat_wrapper,
                                subprecond_wrapper,
                                subsolution_wrapper,
                                subrhs_wrapper,
                                tol,
                                maxits);

    LibmeshPetscCall(VecZeroEntries(solution->vec()));
    VecScatterBeginEnd(this->comm(), reduced_scatter_, reduced_subsolution_, solution->vec(), INSERT_VALUES, SCATTER_REVERSE);
    solution->close();

    // Store the number of linear iterations required to
    // solve and the final residual.
    _n_linear_iterations = rval.first;
    _final_linear_residual = rval.second;

    // Update the system after the solve
    this->update();

    this->get_dof_map().enforce_constraints_exactly(*this);
    this->update();
    this->update_time_derivatives();
    pytec3d::logging::info("Finished linear system solve");
}

void pyTEC3DLinearSystem::post_init()
{
    pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>::post_init();

    auto &equation_systems = this->get_pytec3d_equation_systems();
    auto *solver_parent = equation_systems.get_parent_object();
    auto *linear_solver = this->get_linear_solver();
    if (solver_parent && linear_solver)
    {
        this->solver_config_copy_ = std::make_unique<PetscSolverConfiguration>(
            equation_systems.get_solver_settings(),
            linear_solver);
        this->solver_config_copy_->configure_solver();
    }
    else
        pytec3d::logging::error_and_abort("Solver parent or PETSc linear solver not available.");
}

void pyTEC3DLinearSystem::assemble_system()
{
    libMesh::PerfLog perf_log("Matrix Assembly");

    auto &system = *this;
    libMesh::SparseMatrix<libMesh::Number> &matrix = system.get_system_matrix();
    system.rhs->zero();

    const libMesh::DofMap &dof_map = system.get_dof_map();
    const libMesh::MeshBase &mesh = this->typed_equation_systems().get_mesh();
    auto *solver = this->get_linear_solver();
    const auto &comm = system.comm();

    const bool update_matrix = this->check_update_system_matrix();
    if (update_matrix)
    {
        matrix.zero();
        solver->reuse_preconditioner(false);

        auto &petsc_matrix = static_cast<libMesh::PetscMatrix<libMesh::Number> &>(matrix);
        PetscErrorCode ierr = MatResetHash(petsc_matrix.mat());
        CHKERRABORT(this->comm().get(), ierr);
    }
    else
        solver->reuse_preconditioner(true);

    const TimeScheme scheme = this->typed_equation_systems().time_scheme();
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();
    const std::string rank_str = std::to_string(comm.rank());

    BinaryReader reader;
    const std::filesystem::path mtrx_filepath = emat_folderpath / ("mtrx_" + rank_str + ".emat");

    auto streams = reader.open_input_streams(emat_folderpath,
                                             {"mtrx_" + rank_str + ".emat"});
    auto &mtrx_input_file = *streams[0];

    pytec3d::logging::info("Building global matrix and RHS");

    const auto *nodal_load_map = this->ensure_nodal_load_injector_map();

    ElementBuilder builder(&this->typed_equation_systems());
    builder.set_nodal_load_map(nodal_load_map);
    std::vector<libMesh::dof_id_type> dof_indices;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        dof_map.dof_indices(elem, dof_indices);

        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(mtrx_input_file,
                                            elem_id,
                                            builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
        builder.build_element_rhs();

        auto &rhs = builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];
        auto Ke = apply_time_scheme(builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES],
                                    scheme, system, *system.old_local_solution, dof_indices, rhs);

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
    pytec3d::logging::info("Linear solver output info");
}

void pyTEC3DNonLinearSystem::assemble_system_impl(bool assemble_jacobian, bool assemble_residual)
{
    libMesh::PerfLog perf_log("Matrix Assembly");

    auto &system = *this;
    libMesh::SparseMatrix<libMesh::Number> &matrix = system.get_system_matrix();
    libMesh::SparseMatrix<libMesh::Number> &linear_matrix = system.get_matrix("linear_stiffness");

    const libMesh::DofMap &dof_map = system.get_dof_map();
    const libMesh::MeshBase &mesh = this->typed_equation_systems().get_mesh();
    const auto &comm = system.comm();

    const bool update_matrix = this->check_update_system_matrix();
    const TimeScheme scheme = this->typed_equation_systems().time_scheme();

    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();
    const std::string rank_str = std::to_string(comm.rank());

    if (assemble_jacobian)
    {
        if (update_matrix)
        {
            linear_matrix.zero();

            auto &petsc_linear_matrix = static_cast<libMesh::PetscMatrix<libMesh::Number> &>(linear_matrix);
            PetscErrorCode ierr = MatResetHash(petsc_linear_matrix.mat());
            CHKERRABORT(this->comm().get(), ierr);
        }
        matrix.zero();

        auto &petsc_matrix = static_cast<libMesh::PetscMatrix<libMesh::Number> &>(matrix);
        PetscErrorCode ierr = MatResetHash(petsc_matrix.mat());
        CHKERRABORT(this->comm().get(), ierr);
    }
    if (assemble_residual)
        system.rhs->zero();

    BinaryReader reader;
    const std::filesystem::path mtrx_filepath = emat_folderpath / ("mtrx_" + rank_str + ".emat");

    auto streams = reader.open_input_streams(emat_folderpath,
                                             {"mtrx_" + rank_str + ".emat"});
    auto &mtrx_input_file = *streams[0];

    if (assemble_jacobian && assemble_residual)
        pytec3d::logging::info("Building global Jacobian and residual");
    else if (assemble_jacobian)
        pytec3d::logging::info("Building global Jacobian");
    else
        pytec3d::logging::info("Building global residual");

    const auto *nodal_load_map = this->ensure_nodal_load_injector_map();

    ElementBuilder builder(&this->typed_equation_systems());
    builder.set_nodal_load_map(nodal_load_map);
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::Real> current_solution_values;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        dof_map.dof_indices(elem, dof_indices);

        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(mtrx_input_file,
                                            elem_id,
                                            builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
        builder.build_element_rhs();

        auto &linear_rhs = builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];
        auto linear_ke = apply_time_scheme(builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES],
                                           scheme, system, *system.old_local_solution, dof_indices, linear_rhs);

        system.current_local_solution->get(dof_indices, current_solution_values);
        linear_ke.vector_mult_add(linear_rhs, -1.0, current_solution_values);

        // NewtonSolver solves J * du = r and then updates u <- u - du, so r
        // must be the true residual F(u), not -F(u).
        linear_rhs.scale(-1.0);

        const bool nonlinear_element = builder.nonlinear_element();
        if (nonlinear_element)
        {
            builder.build_element_tangent_stiffness_residual(scheme);
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
        pytec3d::logging::info("Newton-Rhapson solver output info");
}

void pyTEC3DNonLinearSystem::assemble_system()
{
    this->assemble_system_impl(true, true);
}

void pyTEC3DNonLinearSystem::assemble_jacobian_only()
{
    this->assemble_system_impl(true, false);
}

void pyTEC3DNonLinearSystem::assemble_residual_only()
{
    this->assemble_system_impl(false, true);
}

void pyTEC3DNonLinearSystem::solve()
{
    libMesh::NonlinearImplicitSystem::solve();
    this->update_time_derivatives();
    pytec3d::logging::info("Finished nonlinear system solve");
}

void pyTEC3DNonLinearSystem::post_init()
{
    pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>::post_init();

    auto &equation_systems = this->get_pytec3d_equation_systems();
    auto *solver_parent = equation_systems.get_parent_object();
    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    if (!solver_parent || !newton_solver)
        pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::post_init: solver parent or Newton solver not available.");
    newton_solver->require_residual_reduction = false;

    this->nonlinear_solver_config_copy_ = std::make_unique<PetscSolverConfiguration>(
        equation_systems.get_solver_settings(),
        &newton_solver->get_linear_solver());
    this->nonlinear_solver_config_copy_->configure_solver();
}

void pyTEC3DNodalSolutionSystem::assemble_system()
{
    pytec3d::logging::info("Computing nodal solution");

    auto &eq_sys = this->typed_equation_systems();

    const libMesh::MeshBase &mesh = eq_sys.get_mesh();

    BinaryReader reader;
    const std::string rank_str = std::to_string(mesh.comm().rank());
    auto streams = reader.open_input_streams(get_emat_folderpath_from_solver(), {"nodal_sol_mtrx_" + rank_str + ".emat"});

    ElementBuilder builder(&eq_sys);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(*streams[0], elem_id,
                                            builder.get_matrices()[ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES]);
        builder.build_nodal_solution();
    }
    reader.close_file_streams(streams);

    this->solution->close();
    pytec3d::logging::info("Nodal solution updated");
}

void pyTEC3DElementSolutionSystem::assemble_system()
{
    pytec3d::logging::info("Computing elemental solution");

    auto &eq_sys = this->typed_equation_systems();

    const libMesh::MeshBase &mesh = eq_sys.get_mesh();

    BinaryReader reader;
    const std::string rank_str = std::to_string(mesh.comm().rank());
    auto streams = reader.open_input_streams(get_emat_folderpath_from_solver(), {"element_sol_mtrx_" + rank_str + ".emat"});

    ElementBuilder builder(&eq_sys);
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        unsigned int elem_id = 0;
        BinaryReader::read_element_matrices(*streams[0], elem_id,
                                            builder.get_matrices()[ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES]);
        builder.build_element_solution();
    }
    reader.close_file_streams(streams);

    this->solution->close();
    pytec3d::logging::info("Element solution updated");
}

template <typename SysType>
void pyTEC3DSystem<SysType>::init_variables()
{
    auto &equation_systems = this->typed_equation_systems();
    auto element_types_collector = equation_systems.get_element_types_factory();
    const auto type = pyTEC3DSystemBase::system_type_from_name(this->name());

    if (!element_types_collector)
        pytec3d::logging::error_and_abort("Element types are not available for system variables.");

    auto vars_to_subdomain_ids = VariablesMap::associate_variables_with_element_IDs(
        *element_types_collector,
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

        const auto &entry = element_types_collector->get(*subdomain_ids.begin());
        const auto system_variables = VariablesMap::get_system_variables_for_element(*entry.type, entry.domain, type);
        const auto variable = std::find_if(system_variables.begin(), system_variables.end(),
                                           [&pair](const VariablesMap::SystemVariable &v)
                                           { return v.key == pair.first; });

        this->add_variable(
            variable_string,
            variable->order,
            variable->family,
            &subdomain_ids);
    }
}

template class pyTEC3DSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class pyTEC3DSystem<libMesh::TransientSystem<libMesh::NonlinearImplicitSystem>>;
template class pyTEC3DSystem<libMesh::ExplicitSystem>;
template class pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>;
template class pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>;

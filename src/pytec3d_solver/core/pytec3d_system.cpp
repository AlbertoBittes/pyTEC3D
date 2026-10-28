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
#include "report_tables.h"
#include "solver_diagnostics.h"
#include "binary_reader.h"
#include "element_builder.h"
#include "model_element_types.h"
#include "constraints.h"
#include "system_matrix_reduction.h"
#include "variables_map.h"

#include <libmesh/libmesh.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/elem_range.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/diff_solver.h>
#include <libmesh/newton_solver.h>
#include <libmesh/nonlinear_solver.h>
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
    /** @brief Keeps the norm of each full Newton step (before the line search), for the Newton history. */
    class FullStepMonitor final : public libMesh::LinearSolutionMonitor
    {
    public:
        explicit FullStepMonitor(libMesh::Real &norm) : norm_(norm) {}
        void operator()(const libMesh::NumericVector<libMesh::Number> & /*delta_u*/, const libMesh::Real &norm_delta_u,
                        const libMesh::NumericVector<libMesh::Number> & /*u*/, const libMesh::Real & /*norm_u*/,
                        const libMesh::NumericVector<libMesh::Number> & /*residual*/,
                        const libMesh::Real & /*norm_residual*/, const unsigned int /*iteration*/) override
        {
            norm_ = norm_delta_u;
        }

    private:
        libMesh::Real &norm_;
    };

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

// ---- pyTEC3DSystem ----------------------------------------------------------------------------------------

template <typename SysLike>
pyTEC3DSystem<SysLike>::pyTEC3DSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number)
    : SysLike(es, name, number),
      system_type_(pyTEC3DSystemBase::system_type_from_name(name))
{
    this->attach_assemble_function(&pyTEC3DSystem<SysLike>::assemble_dispatch);
}

template <typename SysLike>
void pyTEC3DSystem<SysLike>::solve()
{
    SysLike::solve();
}

template <typename SysLike>
bool pyTEC3DSystem<SysLike>::is_nodal() const
{
    return system_type_ != pyTEC3DSystemBase::ELEMENT_SOLUTION;
}

template <typename SysLike>
void pyTEC3DSystem<SysLike>::post_init()
{
}

template <typename SysLike>
pyTEC3DEquationSystems &pyTEC3DSystem<SysLike>::get_pytec3d_equation_systems()
{
    return this->typed_equation_systems();
}

template <typename SysLike>
const pyTEC3DEquationSystems &pyTEC3DSystem<SysLike>::get_pytec3d_equation_systems() const
{
    return this->typed_equation_systems();
}

template <typename SysLike>
const TimeStep *pyTEC3DSystem<SysLike>::current_time_step() const
{
    return this->typed_equation_systems().current_time_step();
}

template <typename SysLike>
const TimeStep *pyTEC3DSystem<SysLike>::previous_time_step() const
{
    return this->typed_equation_systems().previous_time_step();
}

template <typename SysLike>
pyTEC3DEquationSystems &pyTEC3DSystem<SysLike>::typed_equation_systems()
{
    return static_cast<pyTEC3DEquationSystems &>(this->get_equation_systems());
}

template <typename SysLike>
const pyTEC3DEquationSystems &pyTEC3DSystem<SysLike>::typed_equation_systems() const
{
    return static_cast<const pyTEC3DEquationSystems &>(this->get_equation_systems());
}

template <typename SysLike>
void pyTEC3DSystem<SysLike>::assemble_dispatch(libMesh::EquationSystems &es, const std::string &name)
{
    auto *system = dynamic_cast<pyTEC3DSystemBase *>(&es.get_system(name));
    if (!system)
        pytec3d::logging::error_and_abort("pyTEC3DSystem::assemble_dispatch: system is not a pyTEC3DSystemBase.");
    system->assemble_system();
}

// ---- pyTEC3DDOFSystem -------------------------------------------------------------------------------------

template <typename ImplicitSys>
pyTEC3DDOFSystem<ImplicitSys>::pyTEC3DDOFSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number)
    : pyTEC3DSystem<libMesh::TransientSystem<ImplicitSys>>(es, name, number)
{
    global_constraints_object_ = std::make_unique<GlobalConstraintsObject<ImplicitSys>>(*this);
    this->attach_constraint_object(*global_constraints_object_);
    for (const char *vector_name : {pyTEC3DSystemBase::kVelocity, pyTEC3DSystemBase::kAcceleration,
                                    pyTEC3DSystemBase::kOldVelocity, pyTEC3DSystemBase::kOldAcceleration})
        this->add_vector(vector_name, true, libMesh::GHOSTED);
}

template <typename ImplicitSys>
pyTEC3DDOFSystem<ImplicitSys>::~pyTEC3DDOFSystem() noexcept = default;

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::solve()
{
    libMesh::TransientSystem<ImplicitSys>::solve();
}

template <typename ImplicitSys>
GlobalConstraintsObject<ImplicitSys> *pyTEC3DDOFSystem<ImplicitSys>::get_global_constraints_object() const
{
    return global_constraints_object_.get();
}

template <typename ImplicitSys>
libMesh::LinearSolver<libMesh::Number> *pyTEC3DDOFSystem<ImplicitSys>::get_linear_solver() const
{
    return this->linear_solver.get();
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::prepare_time_step_solve(const TimeStepInput &time_step)
{
    this->typed_equation_systems().prepare_time_step_solve(time_step);
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::init_time_step(const TimeStepInput &time_step_input)
{
    this->prepare_time_step_solve(time_step_input);

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
        pytec3d::logging::table_root(pytec3d::report::format_time_step_table(*time_step, this->typed_equation_systems().time_order()));
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
}

template <typename ImplicitSys>
bool pyTEC3DDOFSystem<ImplicitSys>::check_update_system_matrix()
{
    auto *bc_object = this->get_global_constraints_object();

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

    StepReport &report = equation_systems.step_report();
    if (element_matrices_written_)
    {
        pytec3d::logging::info("Reusing element matrices from previous time step");
        return;
    }

    report.element_matrices_computed = true;
    pytec3d::logging::ScopedTimer timer("Element matrices computed (" +
                                            std::to_string(mesh.n_active_local_elem()) + " local elements)",
                                        &report.seconds_element_matrices);

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
    element_matrices_written_ = true;
}

template <typename ImplicitSys>
const NodalLoadInjectorElemByVarMap *pyTEC3DDOFSystem<ImplicitSys>::ensure_nodal_load_injector_map()
{
    if (!nodal_load_injector_map_)
    {
        auto &equation_systems = this->typed_equation_systems();
        auto *solver_parent = &equation_systems.model();
        if (!solver_parent)
            pytec3d::logging::error_and_abort("pyTEC3DDOFSystem::ensure_nodal_load_injector_map: solver parent not available.");

        auto *vtk_grid = solver_parent->get_active_vtk_grid();
        if (!vtk_grid)
            pytec3d::logging::error_and_abort("pyTEC3DDOFSystem::ensure_nodal_load_injector_map: active VTK grid not available.");

        nodal_load_injector_map_ = std::make_unique<NodalLoadInjectorElemByVarMap>(
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

// ---- pyTEC3DLinearSystem ----------------------------------------------------------------------------------

pyTEC3DLinearSystem::pyTEC3DLinearSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number)
    : pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>(es, name, number)
{
    this->linear_solver = std::make_unique<ReducedPetscLinearSolver>(es.comm());
}

pyTEC3DLinearSystem::~pyTEC3DLinearSystem() noexcept = default;

void pyTEC3DLinearSystem::add_matrices()
{
    libMesh::ImplicitSystem::add_matrices();
    this->get_system_matrix().use_hash_table(true);
}

ReducedPetscLinearSolver &pyTEC3DLinearSystem::reduced_linear_solver()
{
    auto *solver = dynamic_cast<ReducedPetscLinearSolver *>(this->linear_solver.get());
    if (!solver)
        pytec3d::logging::error_and_abort("pyTEC3DLinearSystem: the linear solver is not a ReducedPetscLinearSolver.");
    return *solver;
}

void pyTEC3DLinearSystem::delete_rows_cols_constrained_dofs()
{
    pytec3d::logging::info("Removing constrained DOFs from system matrix and RHS");
    subset_solve_holder_ = std::make_unique<SystemMatrixReduction>(*this);
    this->reduced_linear_solver().restrict_solve_to(&subset_solve_holder_->dof_ids(), libMesh::SUBSET_ZERO);
}

void pyTEC3DNonLinearSystem::delete_rows_cols_constrained_dofs()
{
    pytec3d::logging::info("Removing constrained DOFs from system matrix and RHS");
    subset_solve_holder_ = std::make_unique<SystemMatrixReduction>(*this);
    this->reduced_linear_solver().restrict_solve_to(&subset_solve_holder_->dof_ids(), libMesh::SUBSET_ZERO);
}

void pyTEC3DLinearSystem::solve()
{
    StepReport &report = this->typed_equation_systems().step_report();
    pytec3d::logging::ScopedTimer timer("Linear system assembled and solved", &report.seconds_solve);
    if (this->assemble_before_solve)
        // Assemble the linear system
        this->assemble();

    // Get the user-specified linear solver tolerance
    const auto [maxits, tol] = this->get_linear_solve_parameters();

    if (subset_solve_holder_ && subset_solve_holder_->dof_ids().empty())
        pytec3d::logging::error_and_abort("pyTEC3DLinearSystem::solve: subset solve requested with empty subset dof list.");

    // On the unconstrained DOFs; the factorization reused, or only its numeric part, when possible.
    ReducedPetscLinearSolver &solver = this->reduced_linear_solver();
    const ReducedPetscLinearSolver::Counts before = solver.counts();
    const std::pair<unsigned int, libMesh::Real> rval = solver.solve(*this->matrix, *this->solution, *this->rhs, tol, maxits);
    report.full_factorizations += solver.counts().full - before.full;
    report.numeric_factorizations += solver.counts().numeric - before.numeric;
    report.true_residual = pytec3d::diagnostics::relative_true_residual(solver.reduced_matrix(), solver.reduced_solution(),
                                                                        solver.reduced_rhs());
    report.factorization = pytec3d::diagnostics::factorization_report(solver.ksp());
    this->solution->close();

    // Store the number of linear iterations required to
    // solve and the final residual.
    _n_linear_iterations = rval.first;
    _final_linear_residual = rval.second;

    // Update the system after the solve
    this->update();

    this->get_dof_map().enforce_constraints_exactly(*this);
    this->update();
    this->update_time_derivatives();
    const auto reason = solver.get_converged_reason();
    PetscInt local_unconstrained = 0;
    LibmeshPetscCall(VecGetLocalSize(solver.reduced_rhs(), &local_unconstrained));
    std::size_t unconstrained = static_cast<std::size_t>(local_unconstrained);
    this->comm().sum(unconstrained);
    report.linear_iterations = rval.first;
    report.linear_residual = rval.second;
    report.linear_reason = libMesh::Utility::enum_to_string(reason);
    report.unconstrained_dofs = unconstrained;
    pytec3d::logging::info("Linear solve: {} iteration(s), residual {:.3e}, true residual {:.3e}, {} ({} unconstrained "
                           "DOFs on this rank).",
                           rval.first, rval.second, report.true_residual, libMesh::Utility::enum_to_string(reason),
                           local_unconstrained);
    pytec3d::diagnostics::check_linear_solve(report.factorization, report.true_residual);
    if (reason < 0)
        pytec3d::logging::warn("The linear solver did not converge: {}.", libMesh::Utility::enum_to_string(reason));
}

void pyTEC3DLinearSystem::post_init()
{
    pyTEC3DDOFSystem<libMesh::LinearImplicitSystem>::post_init();

    auto *linear_solver = this->get_linear_solver();
    if (!linear_solver)
        pytec3d::logging::error_and_abort("pyTEC3DLinearSystem::post_init: PETSc linear solver not available.");
    if (this->prefix_with_name())
        linear_solver->init(this->prefix().c_str());
    else
        linear_solver->init();
    linear_solver->init_systems(*this);
    this->solver_config_ = std::make_unique<PetscSolverConfiguration>(
        this->get_pytec3d_equation_systems().get_solver_settings(), linear_solver);
    this->solver_config_->configure_solver();
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
    this->typed_equation_systems().step_report().system_matrix_updated |= update_matrix;
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

// ---- pyTEC3DNonLinearSystem -------------------------------------------------------------------------------

namespace
{
    /** @brief The Newton solver of libMesh with the linear solver of pyTEC3D (ReducedPetscLinearSolver). */
    class pyTEC3DNewtonSolver : public libMesh::NewtonSolver
    {
    public:
        explicit pyTEC3DNewtonSolver(sys_type &system)
            : libMesh::NewtonSolver(system)
        {
            _linear_solver = std::make_unique<ReducedPetscLinearSolver>(system.comm());
        }
    };
} // namespace

pyTEC3DNonLinearSystem::pyTEC3DNonLinearSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number)
    : pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>(es, name, number)
{
    auto newton_solver = std::make_unique<pyTEC3DNewtonSolver>(*this);
    newton_solver->verbose = true;
    newton_solver->require_residual_reduction = true;
    diff_solver = std::move(newton_solver);
    this->nonlinear_solver->matvec = &pyTEC3DNonLinearSystem::matvec_callback;
    this->nonlinear_solver->residual = &pyTEC3DNonLinearSystem::residual_callback;
}

pyTEC3DNonLinearSystem::~pyTEC3DNonLinearSystem() noexcept = default;

ReducedPetscLinearSolver &pyTEC3DNonLinearSystem::reduced_linear_solver()
{
    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    auto *solver = newton_solver ? dynamic_cast<ReducedPetscLinearSolver *>(&newton_solver->get_linear_solver()) : nullptr;
    if (!solver)
        pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem: the Newton linear solver is not a ReducedPetscLinearSolver.");
    return *solver;
}

void pyTEC3DNonLinearSystem::add_matrices()
{
    libMesh::ImplicitSystem::add_matrices();
    this->get_system_matrix().use_hash_table(true);
    this->add_matrix("linear_stiffness");
    this->get_matrix("linear_stiffness").use_hash_table(true);
}

void pyTEC3DNonLinearSystem::matvec_callback(const libMesh::NumericVector<libMesh::Number> & /*X*/,
                                             libMesh::NumericVector<libMesh::Number> *R,
                                             libMesh::SparseMatrix<libMesh::Number> *J,
                                             libMesh::NonlinearImplicitSystem &S)
{
    auto *system = dynamic_cast<pyTEC3DNonLinearSystem *>(&S);
    if (!system)
        pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::matvec_callback: system is not pyTEC3DNonLinearSystem.");
    if (R && J)
    {
        system->assemble_system();
        system->residual_ready_from_matvec_ = true;
        R->close();
        system->record_newton_iterate(R->l2_norm());
    }
    else if (J)
    {
        system->assemble_jacobian_only();
        system->residual_ready_from_matvec_ = false;
    }
    else if (R)
    {
        system->assemble_residual_only();
        system->residual_ready_from_matvec_ = false;
    }
}

void pyTEC3DNonLinearSystem::residual_callback(const libMesh::NumericVector<libMesh::Number> & /*X*/,
                                               libMesh::NumericVector<libMesh::Number> & /*R*/,
                                               libMesh::NonlinearImplicitSystem &S)
{
    auto *system = dynamic_cast<pyTEC3DNonLinearSystem *>(&S);
    if (!system)
        pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::residual_callback: system is not pyTEC3DNonLinearSystem.");
    if (system->residual_ready_from_matvec_)
    {
        system->residual_ready_from_matvec_ = false;
        return;
    }
    system->assemble_residual_only();
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
    this->typed_equation_systems().step_report().system_matrix_updated |= update_matrix;
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
        pytec3d::logging::info("Newton iteration: Jacobian and residual assembled.");
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

void pyTEC3DNonLinearSystem::record_newton_iterate(libMesh::Real residual)
{
    std::vector<NewtonIteration> &history = this->typed_equation_systems().step_report().newton_history;
    NewtonIteration iteration;
    iteration.residual = residual;
    if (!previous_iterate_)
        previous_iterate_ = this->solution->clone();
    else
    {
        // The step from the previous iterate, line search included.
        previous_iterate_->add(-1.0, *this->solution);
        iteration.step = previous_iterate_->l2_norm();
        const libMesh::Real norm = this->solution->l2_norm();
        iteration.relative_step = norm > 0.0 ? iteration.step / norm : iteration.step;
        if (full_newton_step_ > 0.0)
            iteration.step_length = iteration.step / full_newton_step_;
        *previous_iterate_ = *this->solution;
    }
    full_newton_step_ = -1.0;
    history.push_back(iteration);
}

void pyTEC3DNonLinearSystem::solve()
{
    StepReport &report = this->typed_equation_systems().step_report();
    pytec3d::logging::ScopedTimer timer("Nonlinear system solved", &report.seconds_solve);
    previous_iterate_.reset();
    full_newton_step_ = -1.0;
    const ReducedPetscLinearSolver::Counts before = this->reduced_linear_solver().counts();
    libMesh::NonlinearImplicitSystem::solve();
    report.full_factorizations += this->reduced_linear_solver().counts().full - before.full;
    report.numeric_factorizations += this->reduced_linear_solver().counts().numeric - before.numeric;
    // The last iterate: libMesh stops after the residual of the last step, without assembling the Jacobian again.
    if (previous_iterate_)
    {
        this->rhs->close();
        record_newton_iterate(this->rhs->l2_norm());
        if (report.newton_history.size() >= 2 &&
            report.newton_history.back().step == 0.0) // converged at the start of an iteration: no new iterate
            report.newton_history.pop_back();
    }
    previous_iterate_.reset();
    this->update_time_derivatives();
    if (auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get()))
        if (auto *petsc_solver = dynamic_cast<libMesh::PetscLinearSolver<libMesh::Number> *>(
                &newton_solver->get_linear_solver()))
            report.factorization = pytec3d::diagnostics::factorization_report(petsc_solver->ksp());
    pytec3d::diagnostics::check_linear_solve(report.factorization, -1.0);

    const unsigned int result = diff_solver->solve_result();
    std::string reasons;
    for (const auto &[flag, text] : std::initializer_list<std::pair<unsigned int, const char *>>{
             {libMesh::DiffSolver::CONVERGED_ABSOLUTE_RESIDUAL, "absolute residual"},
             {libMesh::DiffSolver::CONVERGED_RELATIVE_RESIDUAL, "relative residual"},
             {libMesh::DiffSolver::CONVERGED_ABSOLUTE_STEP, "absolute step"},
             {libMesh::DiffSolver::CONVERGED_RELATIVE_STEP, "relative step"},
             {libMesh::DiffSolver::DIVERGED_MAX_NONLINEAR_ITERATIONS, "maximum iterations reached"},
             {libMesh::DiffSolver::DIVERGED_BACKTRACKING_FAILURE, "line search failed"},
             {libMesh::DiffSolver::DIVERGED_LINEAR_SOLVER_FAILURE, "linear solver failed"},
             {libMesh::DiffSolver::DIVERGED_NO_REASON, "diverged"}})
        if (result & flag)
            reasons += (reasons.empty() ? "" : ", ") + std::string(text);
    report.newton_iterations = diff_solver->total_outer_iterations();
    report.linear_iterations = diff_solver->total_inner_iterations(); // libMesh counts them per solve
    report.newton_result = reasons.empty() ? "no convergence test met" : reasons;
    pytec3d::logging::info("Newton: {} iteration(s), {} linear iteration(s); {}.", report.newton_iterations,
                           report.linear_iterations, report.newton_result);
    if (result >= libMesh::DiffSolver::DIVERGED_NO_REASON)
        pytec3d::logging::warn("The Newton solver did not converge ({}); the solution is the last iterate.", reasons);
}

void pyTEC3DNonLinearSystem::post_init()
{
    pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>::post_init();

    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    if (!newton_solver)
        pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::post_init: Newton solver not available.");
    // The Newton settings: the tolerances and limits through the parameters libMesh copies into the solver at every
    // solve (NonlinearImplicitSystem::set_solver_parameters), the others on the solver itself.
    const NewtonSettings &newton = this->get_pytec3d_equation_systems().get_solver_settings().nonlinear;
    libMesh::Parameters &parameters = this->get_equation_systems().parameters;
    parameters.set<unsigned int>("nonlinear solver maximum iterations") = newton.max_iterations;
    parameters.set<libMesh::Real>("nonlinear solver absolute residual tolerance") = newton.absolute_residual_tolerance;
    parameters.set<libMesh::Real>("nonlinear solver relative residual tolerance") = newton.relative_residual_tolerance;
    parameters.set<libMesh::Real>("nonlinear solver absolute step tolerance") = newton.absolute_step_tolerance;
    parameters.set<libMesh::Real>("nonlinear solver relative step tolerance") = newton.relative_step_tolerance;
    newton_solver->require_residual_reduction = newton.require_residual_reduction;
    newton_solver->require_finite_residual = newton.require_finite_residual;
    newton_solver->brent_line_search = newton.brent_line_search;
    newton_solver->minsteplength = newton.minimum_step_length;
    newton_solver->linear_tolerance_multiplier = newton.linear_tolerance_multiplier;
    newton_solver->continue_after_max_iterations = newton.continue_after_max_iterations;
    newton_solver->continue_after_backtrack_failure = newton.continue_after_backtrack_failure;
    newton_solver->verbose = newton.verbose;
    newton_solver->linear_solution_monitor = std::make_unique<FullStepMonitor>(full_newton_step_);

    this->solver_config_ = std::make_unique<PetscSolverConfiguration>(
        this->get_pytec3d_equation_systems().get_solver_settings(), &newton_solver->get_linear_solver(), true);
    this->solver_config_->configure_solver();
}

// ---- explicit systems -------------------------------------------------------------------------------------

pyTEC3DExplicitSystem::pyTEC3DExplicitSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number)
    : pyTEC3DSystem<libMesh::ExplicitSystem>(es, name, number)
{
    this->remove_vector("RHS");
}

void pyTEC3DNodalSolutionSystem::assemble_system()
{
    pytec3d::logging::ScopedTimer timer("Nodal fields computed",
                                        &this->typed_equation_systems().step_report().seconds_nodal_fields);

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
}

void pyTEC3DElementSolutionSystem::assemble_system()
{
    pytec3d::logging::ScopedTimer timer("Element fields computed",
                                        &this->typed_equation_systems().step_report().seconds_element_fields);

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
}

template <typename SysType>
void pyTEC3DSystem<SysType>::init_variables()
{
    auto &equation_systems = this->typed_equation_systems();
    const ModelElementTypes &element_types = equation_systems.element_types();
    const auto type = pyTEC3DSystemBase::system_type_from_name(this->name());

    auto vars_to_subdomain_ids = VariablesMap::associate_variables_with_element_IDs(element_types, type);

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

        const auto &entry = element_types.get(*subdomain_ids.begin());
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

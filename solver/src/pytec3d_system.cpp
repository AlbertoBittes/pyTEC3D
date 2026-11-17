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
#include "element_builder.h"
#include "model_element_types.h"
#include "constraints.h"
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

#include <petscmat.h>

#include <exception>
#include <omp.h>


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

    using ElemRecordIndex = std::unordered_map<unsigned int, std::streampos>;

    NodalLoadInjectorElemByVarMap build_nodal_load_injector_map_for_system(
        const libMesh::MeshBase &mesh,
        const libMesh::ImplicitSystem &dof_system,
        const pytec3d::ModelArrays &arrays)
    {
        NodalLoadInjectorElemByVarMap nodal_load_map;

        const auto dof_nodal_load_map =
            VariablesMap::get_dof_variable_nodal_load_map_for_system(dof_system);

        for (const auto &entry : dof_nodal_load_map)
        {
            const std::string &dof_variable_name = entry.first;
            const std::string &nodal_load_array_name = entry.second;

            if (!dof_system.has_variable(dof_variable_name))
                continue;

            const pytec3d::RealArray *load_array = arrays.point_field(nodal_load_array_name);
            if (!load_array)
                continue;

            const unsigned int dof_var_num = dof_system.variable_number(dof_variable_name);
            const unsigned int sys_num = dof_system.number();
            const auto &dof_var = dof_system.variable(dof_var_num);
            const unsigned int n_components = load_array->n_components();

            // One element per loaded node (a node with a non-zero component), which adds all its components.
            std::map<libMesh::dof_id_type, libMesh::dof_id_type> local_injector_elem_by_node;
            for (const auto &elem : mesh.active_local_element_ptr_range())
            {
                if (!dof_var.active_on_subdomain(elem->subdomain_id()))
                    continue;

                for (unsigned int i = 0; i < elem->n_nodes(); ++i)
                {
                    const auto &node = elem->node_ref(i);
                    if (node.n_dofs(sys_num, dof_var_num) == 0)
                        continue;
                    const auto point = static_cast<std::size_t>(node.id());
                    bool loaded = false;
                    for (unsigned int c = 0; c < n_components && !loaded; ++c)
                        loaded = (*load_array)(point, c) != 0.0;
                    if (loaded)
                        local_injector_elem_by_node.emplace(node.id(), static_cast<libMesh::dof_id_type>(elem->id()));
                }
            }

            auto global_injector_elem_by_node = local_injector_elem_by_node;
            mesh.comm().set_union(global_injector_elem_by_node);
            auto &injector_elem_by_node = nodal_load_map[dof_variable_name];
            for (const auto &injector_entry : global_injector_elem_by_node)
                injector_elem_by_node[injector_entry.first] = injector_entry.second;
        }

        return nodal_load_map;
    }

}

namespace
{
    /** @brief Elements computed on the threads, then added in element order: the work of a chunk at a time in memory. */
    constexpr std::size_t kChunk = 4096;

    /**
     * @brief fn(builder, i) for i in [begin, end) on the threads of libMesh (n_threads()), each thread with its own
     * ElementBuilder (made by @p make). An exception of a thread is rethrown after the loop.
     */
    template <typename Make, typename Fn>
    void parallel_elements(std::size_t begin, std::size_t end, Make &&make, Fn &&fn)
    {
        std::exception_ptr failure;
        const int threads = std::max(1, static_cast<int>(libMesh::n_threads()));
#pragma omp parallel num_threads(threads)
        {
            ElementBuilder builder = make();
#pragma omp for schedule(dynamic, 16)
            for (std::size_t i = begin; i < end; ++i)
            {
                try
                {
                    fn(builder, i);
                }
                catch (...)
                {
#pragma omp critical(pytec3d_element_failure)
                    if (!failure)
                        failure = std::current_exception();
                }
            }
        }
        if (failure)
            std::rethrow_exception(failure);
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
    const bool settings_changed = solver_settings_changed_;
    if (settings_changed)
    {
        solver_settings_changed_ = false;
        // another KSP or PC type: a full factorization
        this->reduced_linear_solver().invalidate_pattern();
        pytec3d::logging::info("New linear solver settings: updating system matrix");
    }

    // the constraints are checked also with new settings: new constrained DOFs restrict the solve
    auto *bc_object = this->get_global_constraints_object();

    if (bc_object && bc_object->check_update_boundary_conditions())
    {
        pytec3d::logging::info("New constraints: updating system matrix");
        return true;
    }
    if (settings_changed)
        return true;

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
            // The same step size to the rounding of the step times: each time carries about eps |t|, so two dt computed
            // from them (e.g. 10.2802 - 10.2800 and 10.2804 - 10.2802) differ by a few eps |t| although the step is the
            // same; times accumulated as t += dt drift further (the relative term). Within this the factorization of the
            // previous step is kept: the matrix differs from it by less than the rounding of the data.
            const libMesh::Real eps = std::numeric_limits<libMesh::Real>::epsilon();
            const libMesh::Real dt_tolerance =
                std::max(8.0 * eps * std::max(std::abs(current_step->time), std::abs(old_step->time)),
                         1e-10 * std::max(std::abs(current_step->dt), std::abs(old_step->dt)));
            // Two static steps: the matrix (K, c_v = c_a = 0) does not depend on dt.
            const bool dt_matters = current_step->transient_effects || old_step->transient_effects;
            const bool dt_differs = dt_matters && std::abs(current_step->dt - old_step->dt) > dt_tolerance;
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
const std::vector<const libMesh::Elem *> &pyTEC3DDOFSystem<ImplicitSys>::local_elements()
{
    if (local_elements_.empty())
        for (const libMesh::Elem *elem : this->typed_equation_systems().get_mesh().active_local_element_ptr_range())
            local_elements_.push_back(elem);
    return local_elements_;
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::write_element_matrices()
{
    auto &equation_systems = this->typed_equation_systems();
    const auto &comm = equation_systems.get_mesh().comm();
    StepReport &report = equation_systems.step_report();
    if (element_matrices_written_)
    {
        pytec3d::logging::info("Reusing the operators of the model (K, C, M)");
        return;
    }
    report.element_matrices_computed = true;
    const auto &elements = local_elements();
    pytec3d::logging::ScopedTimer timer("Operators assembled (" + std::to_string(elements.size()) + " local elements)",
                                        &report.seconds_element_matrices);

    // K, C, M once, global and unconstrained (the constraints are applied to the global systems, ReducedPetscLinearSolver);
    // the operators of the static nodal and element fields kept in memory, applied at each output; the rhs bases of the
    // elements (the integrated test functions of their rhs terms), on which the loads of each step are evaluated.
    nodal_operators_.assign(elements.size(), {});
    element_operators_.assign(elements.size(), {});
    rhs_bases_.assign(elements.size(), {});
    rhs_cached_.assign(elements.size(), 0);
    libMesh::SparseMatrix<libMesh::Number> &K = this->effective_matrix();
    libMesh::SparseMatrix<libMesh::Number> *C = this->have_matrix("damping") ? &this->get_matrix("damping") : nullptr;
    libMesh::SparseMatrix<libMesh::Number> *M = this->have_matrix("mass") ? &this->get_matrix("mass") : nullptr;
    K.zero();
    if (C)
        C->zero();
    if (M)
        M->zero();
    bool local_damping = false, local_mass = false;

    struct Built
    {
        std::vector<libMesh::dof_id_type> dofs;
        std::vector<libMesh::DenseMatrix<libMesh::Number>> dof, nodal, element;
    };
    std::vector<Built> built;
    for (std::size_t begin = 0; begin < elements.size(); begin += kChunk)
    {
        const std::size_t end = std::min(elements.size(), begin + kChunk);
        built.clear();
        built.resize(end - begin);
        parallel_elements(begin, end, [&] { return ElementBuilder(&equation_systems); },
                          [&](ElementBuilder &builder, std::size_t i)
                          {
                              builder.set_element_data(elements[i]);
                              builder.build_element();
                              rhs_cached_[i] = builder.build_rhs_basis(rhs_bases_[i]) ? 1 : 0;
                              Built &out = built[i - begin];
                              this->get_dof_map().dof_indices(elements[i], out.dofs);
                              auto &m = builder.get_matrices();
                              out.dof = std::move(m[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES]);
                              out.nodal = std::move(m[ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES]);
                              out.element = std::move(m[ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES]);
                          });
        for (std::size_t i = begin; i < end; ++i)
        {
            Built &out = built[i - begin];
            K.add_matrix(out.dof[0], out.dofs);
            if (out.dof.size() > 1 && out.dof[1].m() > 0 && C)
            {
                C->add_matrix(out.dof[1], out.dofs);
                local_damping = true;
            }
            if (out.dof.size() > 2 && out.dof[2].m() > 0 && M)
            {
                M->add_matrix(out.dof[2], out.dofs);
                local_mass = true;
            }
            // A zero operator is stored empty (0 x 0, ElementBuilder::build_element()): no memory for it.
            nodal_operators_[i] = std::move(out.nodal);
            element_operators_[i] = std::move(out.element);
        }
    }
    K.close();
    if (C)
        C->close();
    if (M)
        M->close();
    comm.max(local_damping);
    comm.max(local_mass);
    has_damping_ = local_damping;
    has_mass_ = local_mass;
    keff_c_v_ = keff_c_a_ = 0.0; // effective_matrix() holds K
    element_matrices_written_ = true;
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::apply_constraints(bool rows_changed)
{
    this->reduced_linear_solver().set_constraints(this->get_dof_map(), rows_changed);
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::update_effective_matrix()
{
    // K_eff += (c_v - c_v') C + (c_a - c_a') M: the same pattern (all the matrices have the sparsity of the DofMap).
    const TimeScheme scheme = this->typed_equation_systems().time_scheme();
    const libMesh::Real c_v = has_damping_ && scheme.transient() ? scheme.c_v : 0.0;
    const libMesh::Real c_a = has_mass_ && scheme.transient() ? scheme.c_a : 0.0;
    libMesh::SparseMatrix<libMesh::Number> &keff = this->effective_matrix();
    if (c_v != keff_c_v_)
        keff.add(c_v - keff_c_v_, this->get_matrix("damping"));
    if (c_a != keff_c_a_)
        keff.add(c_a - keff_c_a_, this->get_matrix("mass"));
    keff_c_v_ = c_v;
    keff_c_a_ = c_a;
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::assemble_rhs(libMesh::NumericVector<libMesh::Number> &F)
{
    auto &equation_systems = this->typed_equation_systems();
    const auto *nodal_load_map = this->ensure_nodal_load_injector_map();
    const auto &elements = local_elements();
    F.zero();

    // The element right-hand sides on the threads, added in element order.
    struct Built
    {
        std::vector<libMesh::dof_id_type> dofs;
        libMesh::DenseVector<libMesh::Number> rhs;
    };
    std::vector<Built> built;
    for (std::size_t begin = 0; begin < elements.size(); begin += kChunk)
    {
        const std::size_t end = std::min(elements.size(), begin + kChunk);
        built.clear();
        built.resize(end - begin);
        parallel_elements(begin, end,
                          [&]
                          {
                              ElementBuilder builder(&equation_systems);
                              builder.set_nodal_load_map(nodal_load_map);
                              return builder;
                          },
                          [&](ElementBuilder &builder, std::size_t i)
                          {
                              builder.set_element_data(elements[i]);
                              // The loads now on the cached geometry of the rhs (its basis), else the full evaluation.
                              if (rhs_cached_.size() == elements.size() && rhs_cached_[i])
                                  builder.build_element_rhs_from_basis(rhs_bases_[i]);
                              else
                                  builder.build_element_rhs();
                              Built &out = built[i - begin];
                              this->get_dof_map().dof_indices(elements[i], out.dofs);
                              out.rhs = std::move(builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0]);
                          });
        for (std::size_t i = begin; i < end; ++i)
            F.add_vector(built[i - begin].rhs, built[i - begin].dofs);
    }
    F.close();

    // The history of the step: - C h_v - M h_a, h = hu u_n + hv v_n + ha a_n (time_step.h).
    const TimeScheme scheme = equation_systems.time_scheme();
    if (!scheme.transient())
        return;
    const auto history = [&](const char *matrix, libMesh::Real hu, libMesh::Real hv, libMesh::Real ha)
    {
        if (hu == 0.0 && hv == 0.0 && ha == 0.0)
            return;
        std::unique_ptr<libMesh::NumericVector<libMesh::Number>> h = this->solution->zero_clone();
        h->add(hu, *this->old_local_solution);
        h->add(hv, this->get_vector(pyTEC3DSystemBase::kOldVelocity));
        h->add(ha, this->get_vector(pyTEC3DSystemBase::kOldAcceleration));
        h->close();
        std::unique_ptr<libMesh::NumericVector<libMesh::Number>> t = this->solution->zero_clone();
        this->get_matrix(matrix).vector_mult(*t, *h);
        F.add(-1.0, *t);
    };
    if (has_damping_)
        history("damping", scheme.hu_v, scheme.hv_v, scheme.ha_v);
    if (has_mass_)
        history("mass", scheme.hu_a, scheme.hv_a, scheme.ha_a);
    F.close();
}

template <typename ImplicitSys>
const NodalLoadInjectorElemByVarMap *pyTEC3DDOFSystem<ImplicitSys>::ensure_nodal_load_injector_map()
{
    if (!nodal_load_injector_map_)
    {
        auto &equation_systems = this->typed_equation_systems();
        nodal_load_injector_map_ = std::make_unique<NodalLoadInjectorElemByVarMap>(
            build_nodal_load_injector_map_for_system(
                equation_systems.get_mesh(),
                *this,
                equation_systems.model().arrays()));
    }

    return nodal_load_injector_map_.get();
}

template <typename ImplicitSys>
void pyTEC3DDOFSystem<ImplicitSys>::apply_solver_settings(const SolverSettings &settings)
{
    if (!solver_config_)
        pytec3d::logging::error_and_abort("pyTEC3DDOFSystem::apply_solver_settings: the system is not initialized.");
    if (solver_config_->set_settings(settings))
        solver_settings_changed_ = true;
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
    // The system matrix holds K_eff; C and M (the transient element types) for its updates and the history of a step.
    libMesh::ImplicitSystem::add_matrices();
    const unsigned int time_order = this->typed_equation_systems().time_order();
    if (time_order >= 1)
        this->add_matrix("damping");
    if (time_order >= 2)
        this->add_matrix("mass");
}

ReducedPetscLinearSolver &pyTEC3DLinearSystem::reduced_linear_solver()
{
    auto *solver = dynamic_cast<ReducedPetscLinearSolver *>(this->linear_solver.get());
    if (!solver)
        pytec3d::logging::error_and_abort("pyTEC3DLinearSystem: the linear solver is not a ReducedPetscLinearSolver.");
    return *solver;
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
    this->solver_config_ = std::make_unique<PetscSolverConfiguration>(linear_solver);
}

void pyTEC3DLinearSystem::assemble_system()
{
    // K_eff from the operators of the model (assembled once); only the right-hand side is assembled at each step.
    const bool update_matrix = this->check_update_system_matrix();
    this->typed_equation_systems().step_report().system_matrix_updated |= update_matrix;
    if (update_matrix)
        this->update_effective_matrix();
    this->get_linear_solver()->reuse_preconditioner(!update_matrix);
    pytec3d::logging::info("Building the right-hand side");
    this->assemble_rhs(*this->rhs);
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
    // The Newton steps: homogeneous constraints, the residual projected by the assembly.
    this->reduced_linear_solver().set_homogeneous(true);
    // F - C h_v - M h_a of the step, assembled once per step (it follows the loads, not the iterate).
    this->add_vector("linear_rhs", false);
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
    // The system matrix holds the Jacobian; linear_stiffness K_eff; C and M as in the linear system.
    libMesh::ImplicitSystem::add_matrices();
    this->add_matrix("linear_stiffness");
    const unsigned int time_order = this->typed_equation_systems().time_order();
    if (time_order >= 1)
        this->add_matrix("damping");
    if (time_order >= 2)
        this->add_matrix("mass");
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
    auto &system = *this;
    auto &equation_systems = this->typed_equation_systems();
    const TimeScheme scheme = equation_systems.time_scheme();
    libMesh::SparseMatrix<libMesh::Number> &jacobian = system.get_system_matrix();
    libMesh::SparseMatrix<libMesh::Number> &keff = this->effective_matrix();

    if (assemble_jacobian)
    {
        const bool update_matrix = this->check_update_system_matrix();
        equation_systems.step_report().system_matrix_updated |= update_matrix;
        if (update_matrix)
            this->update_effective_matrix();
    }
    if (assemble_jacobian && assemble_residual)
        pytec3d::logging::info("Building global Jacobian and residual");
    else if (assemble_jacobian)
        pytec3d::logging::info("Building global Jacobian");
    else
        pytec3d::logging::info("Building global residual");

    // The law terms of the elements with laws, on the threads: J_laws and R_laws at the current solution.
    struct Built
    {
        std::vector<libMesh::dof_id_type> dofs;
        libMesh::DenseMatrix<libMesh::Number> tangent;
        libMesh::DenseVector<libMesh::Number> residual;
    };
    const std::vector<const libMesh::Elem *> &elements = this->law_elements();
    if (assemble_jacobian)
    {
        // J = K_eff + J_laws: the same pattern (the sparsity of the DofMap).
        jacobian.zero();
        jacobian.close();
        jacobian.add(1.0, keff);
    }
    if (assemble_residual)
    {
        // R = K_eff u - F (+ R_laws): the residual the Newton step u <- u - du drives to zero.
        keff.vector_mult(*system.rhs, *system.solution);
        system.rhs->add(-1.0, system.get_vector("linear_rhs"));
    }
    const bool picard = equation_systems.get_solver_settings().nonlinear.picard;
    std::vector<Built> built;
    for (std::size_t begin = 0; begin < elements.size(); begin += kChunk)
    {
        const std::size_t end = std::min(elements.size(), begin + kChunk);
        built.clear();
        built.resize(end - begin);
        parallel_elements(begin, end, [&] { return ElementBuilder(&equation_systems); },
                          [&](ElementBuilder &builder, std::size_t i)
                          {
                              const libMesh::Elem *elem = elements[i];
                              builder.set_element_data(elem);
                              builder.build_element_tangent_stiffness_residual(scheme, picard);
                              Built &out = built[i - begin];
                              system.get_dof_map().dof_indices(elem, out.dofs);
                              out.tangent = std::move(builder.get_matrices()[ElementBuilder::MatricesSet::TANGENT_STIFFNESS_MATRICES][0]);
                              out.residual = std::move(builder.get_vectors()[ElementBuilder::VectorsSet::RESIDUAL_VECTORS][0]);
                              // The stored rates of the law outputs under dt(): hv_d q'_n + ha_d q''_n (the builder adds
                              // the rest).
                              if (!scheme.transient())
                                  return;
                              if (const auto it = old_law_rates_.find(elem->id()); it != old_law_rates_.end())
                              {
                                  const std::array<std::array<libMesh::Real, 2>, 2> h = {{{scheme.hv_v, scheme.ha_v},
                                                                                          {scheme.hv_a, scheme.ha_a}}};
                                  for (unsigned int d = 0; d < 2; ++d)
                                  {
                                      if (it->second.v[d].size() == out.residual.size())
                                          out.residual.add(h[d][0], it->second.v[d]);
                                      if (it->second.a[d].size() == out.residual.size())
                                          out.residual.add(h[d][1], it->second.a[d]);
                                  }
                              }
                          });
        for (std::size_t i = begin; i < end; ++i)
        {
            const Built &out = built[i - begin];
            if (assemble_jacobian)
                jacobian.add_matrix(out.tangent, out.dofs);
            if (assemble_residual)
                system.rhs->add_vector(out.residual, out.dofs);
        }
    }
    if (assemble_jacobian)
        jacobian.close();
    if (assemble_residual)
    {
        system.rhs->close();
        // P^T R at the unconstrained DOFs, 0 at the constrained ones: the residual (and its norm) of the reduced system.
        this->reduced_linear_solver().project_residual(*system.rhs);
    }
    if (assemble_jacobian && assemble_residual)
        pytec3d::logging::info("Newton iteration: Jacobian and residual assembled.");
}

const std::vector<const libMesh::Elem *> &pyTEC3DNonLinearSystem::law_elements()
{
    if (!law_elements_ready_)
    {
        ElementBuilder builder(&this->typed_equation_systems());
        for (const libMesh::Elem *elem : this->local_elements())
        {
            builder.set_element_data(elem);
            if (builder.nonlinear_element())
                law_elements_.push_back(elem);
        }
        law_elements_ready_ = true;
    }
    return law_elements_;
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

void pyTEC3DNonLinearSystem::init_time_step(const TimeStepInput &time_step)
{
    pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>::init_time_step(time_step);
    if (const auto *step = this->current_time_step(); step && step->updated)
        old_law_rates_ = law_rates_;
}

void pyTEC3DNonLinearSystem::update_time_derivatives()
{
    pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>::update_time_derivatives();

    // q' = c_v q(u) + hu_v q(u_n) + hv_v q'_n + ha_v q''_n and q'' likewise (as the DOFs): the rates of the law outputs
    // under dt(), per element; zero on a static step.
    const TimeScheme scheme = this->typed_equation_systems().time_scheme();
    const libMesh::MeshBase &mesh = this->typed_equation_systems().get_mesh();
    ElementBuilder builder(&this->typed_equation_systems());
    law_rates_.clear();
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        if (!builder.rate_residual_element())
            continue;
        LawRates rates;
        const auto old = old_law_rates_.find(elem->id());
        for (unsigned int d = 0; d < 2; ++d)
        {
            const libMesh::DenseVector<libMesh::Number> q = builder.law_output(*this->current_local_solution, d + 1);
            const libMesh::DenseVector<libMesh::Number> q_n = builder.law_output(*this->old_local_solution, d + 1);
            const auto rate = [&](libMesh::Real c, libMesh::Real hu, libMesh::Real hv, libMesh::Real ha)
            {
                libMesh::DenseVector<libMesh::Number> x(q.size());
                if (!scheme.transient())
                    return x;
                x.add(c, q);
                x.add(hu, q_n);
                if (old != old_law_rates_.end() && old->second.v[d].size() == x.size())
                {
                    x.add(hv, old->second.v[d]);
                    x.add(ha, old->second.a[d]);
                }
                return x;
            };
            rates.v[d] = rate(scheme.c_v, scheme.hu_v, scheme.hv_v, scheme.ha_v);
            rates.a[d] = rate(scheme.c_a, scheme.hu_a, scheme.hv_a, scheme.ha_a);
        }
        law_rates_.emplace(elem->id(), std::move(rates));
    }
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
    this->assemble_rhs(this->get_vector("linear_rhs"));
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
    newton_solver->linear_solution_monitor = std::make_unique<FullStepMonitor>(full_newton_step_);

    this->solver_config_ = std::make_unique<PetscSolverConfiguration>(&newton_solver->get_linear_solver(), true);
}

void pyTEC3DNonLinearSystem::apply_solver_settings(const SolverSettings &settings)
{
    pyTEC3DDOFSystem<libMesh::NonlinearImplicitSystem>::apply_solver_settings(settings);

    auto *newton_solver = dynamic_cast<libMesh::NewtonSolver *>(this->diff_solver.get());
    if (!newton_solver)
        pytec3d::logging::error_and_abort("pyTEC3DNonLinearSystem::apply_solver_settings: Newton solver not available.");
    // The tolerances and limits through the parameters libMesh copies into the solver at every solve
    // (NonlinearImplicitSystem::set_solver_parameters), the others on the solver itself.
    const NewtonSettings &newton = settings.nonlinear;
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
}

// ---- explicit systems -------------------------------------------------------------------------------------

pyTEC3DExplicitSystem::pyTEC3DExplicitSystem(libMesh::EquationSystems &es, const std::string &name, unsigned int number)
    : pyTEC3DSystem<libMesh::ExplicitSystem>(es, name, number)
{
    this->remove_vector("RHS");
}

namespace
{
    /** @brief The nodal or element fields of the local elements, from their operators kept by the DOF system, on the
     * threads. */
    void build_field_system(pyTEC3DEquationSystems &eq_sys, const pyTEC3DDOFSystemBase::FieldOperators &operators, bool nodal)
    {
        const libMesh::MeshBase &mesh = eq_sys.get_mesh();
        std::vector<const libMesh::Elem *> elements;
        for (const libMesh::Elem *elem : mesh.active_local_element_ptr_range())
            elements.push_back(elem);
        if (operators.size() != elements.size())
            pytec3d::logging::error_and_abort("Field operators: not computed for the local elements of the mesh.");
        // The local array of the field system, written in place by the threads (each DOF by its only element).
        libMesh::System &system = eq_sys.get_system(pyTEC3DSystemBase::system_name_from_type(
            nodal ? pyTEC3DSystemBase::NODAL_SOLUTION : pyTEC3DSystemBase::ELEMENT_SOLUTION));
        auto &solution = libMesh::cast_ref<libMesh::PetscVector<libMesh::Number> &>(*system.solution);
        libMesh::Number *local_array = solution.get_array();
        parallel_elements(0, elements.size(), [&]
                          {
                              ElementBuilder builder(&eq_sys);
                              builder.set_field_output(local_array);
                              return builder;
                          },
                          [&](ElementBuilder &builder, std::size_t i)
                          {
                              builder.set_element_data(elements[i]);
                              builder.set_field_operators(&operators[i]);
                              if (nodal)
                                  builder.build_nodal_solution();
                              else
                                  builder.build_element_solution();
                          });
        solution.restore_array();
    }
} // namespace

void pyTEC3DNodalSolutionSystem::assemble_system()
{
    auto &eq_sys = this->typed_equation_systems();
    pytec3d::logging::ScopedTimer timer("Nodal fields computed", &eq_sys.step_report().seconds_nodal_fields);
    build_field_system(eq_sys, eq_sys.get_dof_solution_system().nodal_field_operators(), true);
    this->solution->close();
}

void pyTEC3DElementSolutionSystem::assemble_system()
{
    auto &eq_sys = this->typed_equation_systems();
    pytec3d::logging::ScopedTimer timer("Element fields computed", &eq_sys.step_report().seconds_element_fields);
    build_field_system(eq_sys, eq_sys.get_dof_solution_system().element_field_operators(), false);
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

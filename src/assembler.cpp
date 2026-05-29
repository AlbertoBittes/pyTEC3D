/**
 * @file assembler.cpp
 * @brief Implementation of the global FE matrix and RHS assembler.
 *
 * Reads rank-local .emat cache files, applies the generalised-theta time-
 * integration scheme, and scatter-adds constrained element contributions into
 * the global libMesh sparse system.
 */
#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/transient_system.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/dof_map.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/elem.h>
#include <libmesh/perf_log.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>

#include "assembler.h"
#include "apec_m_system.h"
#include "apec_m_equation_systems.h"
#include "apec_m_solver.h"
#include "element_builder.h"
#include "element_types.h"
#include "global_constraints_object.h"
#include "binary_reader.h"
#include "apec_m_logging.h"

using libMesh::dof_id_type;
using libMesh::DofMap;
using libMesh::MeshBase;
using libMesh::Number;
using libMesh::PerfLog;
using libMesh::PetscLinearSolver;
using libMesh::Real;
using libMesh::SparseMatrix;
using libMesh::System;

namespace
{
    /**
     * @brief Iterates over all active local elements, calling prepare_fn() then build_fn() for each.
     *
     * Centralises the element loop pattern shared by matrix, RHS, and solution builds.
     */
    template <typename PrepareFn, typename BuildFn>
    void run_for_elements(const MeshBase &mesh, ElementBuilder &builder, PrepareFn &&prepare_fn, BuildFn &&build_fn)
    {
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            builder.set_element_data(elem);
            prepare_fn();
            build_fn();
        }
    }

    /** @brief Returns the rank-local .emat cache directory (cwd/emat). */
    std::filesystem::path get_emat_folderpath_from_solver()
    {
        return std::filesystem::current_path() / "emat";
    }
} // namespace

AssemblerLinear::AssemblerLinear(ApecMEquationSystems *eq_systems)
    : FEParentObjectBase<ApecMEquationSystems>(eq_systems)
{
}

void AssemblerLinear::assemble()
{
    PerfLog perf_log("Matrix Assembly");

    auto eq_sys = this->get_parent_object();

    auto builder = ElementBuilder(this->get_parent_object());

    auto &system = eq_sys->get_system<ApecMLinearSystem>("DOF Solution");
    SparseMatrix<Number> &matrix = system.get_system_matrix();
    system.rhs->zero();

    const DofMap &dof_map = system.get_dof_map();
    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = system.comm();

    const bool update_matrix = this->check_update_system_matrix();
    std::vector<dof_id_type> dof_indices;
    std::vector<Real> old_solution_values;
    std::vector<Real> old_time_derivative_values;

    const auto time_step = system.time_step.get();
    const Real dt = time_step->dt;
    const Real theta = time_step->theta;
    const bool transient_effect = time_step->transient_effects;
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();
    const std::string rank_str = std::to_string(comm.rank());

    BinaryReader reader;
    auto streams = reader.open_input_streams(emat_folderpath,
                                             {"mtrx_" + rank_str + ".emat",
                                              "rhs_" + rank_str + ".emat"});
    auto &mtrx_input_file = *streams[0];
    auto &rhs_input_file = *streams[1];

    apec::logging::info("Building global matrix and RHS");

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices(elem, dof_indices);

        BinaryReader::read_element_matrices(builder, mtrx_input_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
        auto &Ke = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][0];
        auto &Ce = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][1];

        if (transient_effect)
        {
            Ce.scale(1 / (dt * theta));
            Ke += Ce;
        }
        BinaryReader::read_element_vectors(builder, rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        auto &rhs = builder.get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];

        if (transient_effect)
        {
            system.old_local_solution->get(dof_indices, old_solution_values);
            Ce.vector_mult_add(rhs, 1, old_solution_values);

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

bool AssemblerLinear::check_update_system_matrix()
{
    auto eq_sys = this->get_parent_object();

    auto &system = eq_sys->get_system<ApecMLinearSystem>("DOF Solution");
    auto bc_object = system.get_global_constraints_object();
    auto *solver = system.get_petsc_solver();

    if (bc_object && bc_object->check_update_boundary_conditions())
    {
        apec::logging::info("New constraints: updating system matrix");
        system.get_system_matrix().zero();
        if (solver)
            solver->reuse_preconditioner(false);
        return true;
    }

    bool update_matrix = false;
    const auto *current_step = system.time_step.get();
    const auto *old_step = system.old_time_step.get();

    if (!current_step)
    {
        update_matrix = true;
    }
    else if (current_step->updated)
    {
        if (!old_step)
            update_matrix = true;
        else
        {
            const Real current_dt = current_step->dt;
            const Real old_dt = old_step->dt;
            const Real current_theta = current_step->theta;
            const Real old_theta = old_step->theta;
            const bool current_transient_effect = current_step->transient_effects;
            const bool old_transient_effect = old_step->transient_effects;

            const bool dt_differs = std::abs(current_dt - old_dt) > std::numeric_limits<Real>::epsilon();
            const bool theta_differs = std::abs(current_theta - old_theta) > std::numeric_limits<Real>::epsilon();
            const bool transient_effect_differs = current_transient_effect != old_transient_effect;

            if (dt_differs || theta_differs || transient_effect_differs)
                update_matrix = true;
        }
    }

    if (update_matrix)
    {
        apec::logging::info("Updating system matrix for current time-step");
        system.get_system_matrix().zero();
        if (solver)
            solver->reuse_preconditioner(false);
    }
    else
    {
        apec::logging::info("Reusing system matrix from previous time step");
        if (solver)
            solver->reuse_preconditioner(true);
    }

    return update_matrix;
}

void AssemblerLinear::write_element_matrices_and_rhs()
{
    auto eq_sys = this->get_parent_object();

    auto builder = ElementBuilder(this->get_parent_object());

    auto &system = eq_sys->get_system<ApecMLinearSystem>("DOF Solution");
    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = system.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;

    const std::vector<std::string> filenames = {
        "mtrx_" + rank_str + ".emat",
        "nodal_sol_mtrx_" + rank_str + ".emat",
        "element_sol_mtrx_" + rank_str + ".emat",
        "base_rhs_" + rank_str + ".emat"};

    const bool emat_files_exist = std::all_of(filenames.begin(), filenames.end(), [&](const std::string &filename)
                                              { return std::filesystem::exists(emat_folderpath / filename); });

    if (!emat_files_exist)
    {
        apec::logging::info("Computing element matrices and RHS for current time-step");

        auto streams = reader.open_output_streams(emat_folderpath, filenames);
        auto &mtrx_output_file = *streams[0];
        auto &nodal_sol_mtrx_output_file = *streams[1];
        auto &element_sol_mtrx_output_file = *streams[2];
        auto &base_rhs_output_file = *streams[3];

        run_for_elements(mesh, builder, []() {}, [&]()
                         {
                             builder.build_element();

                             BinaryReader::write_element_matrices(builder, mtrx_output_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
                             BinaryReader::write_element_matrices(builder, nodal_sol_mtrx_output_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES);
                             BinaryReader::write_element_matrices(builder, element_sol_mtrx_output_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES);
                             BinaryReader::write_element_vectors(builder, base_rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS); });

        reader.close_file_streams(streams);
    }
    else
    {
        apec::logging::info("Reusing element matrices from previous time step");
    }

    const std::vector<std::string> inp_filenames = {"base_rhs_" + rank_str + ".emat"};
    const std::vector<std::string> out_filenames = {"rhs_" + rank_str + ".emat"};

    auto streams_in = reader.open_input_streams(emat_folderpath, inp_filenames);
    auto streams_out = reader.open_output_streams(emat_folderpath, out_filenames);

    auto &base_rhs_input_file = *streams_in[0];
    auto &rhs_output_file = *streams_out[0];

    apec::logging::info("Applying loads to elements RHS");

    builder.clear_nodal_loads_dof_cache();

    run_for_elements(mesh, builder, [&]()
                     { BinaryReader::read_element_vectors(builder, base_rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS); }, [&]()
                     {
                         builder.build_element_rhs();
                         BinaryReader::write_element_vectors(builder, rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS); });

    reader.close_file_streams(streams_in);
    reader.close_file_streams(streams_out);
}

void AssemblerLinear::write_nodal_solution()
{
    apec::logging::info("Computing nodal solution");

    auto eq_sys = this->get_parent_object();

    auto builder = ElementBuilder(this->get_parent_object());

    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = mesh.comm();
    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"nodal_sol_mtrx_" + rank_str + ".emat"};
    auto inp_streams = reader.open_input_streams(emat_folderpath, inp_filenames);
    auto &nodal_solution_mtrx_input_file = *inp_streams[0];

    run_for_elements(mesh, builder, [&]()
                     { BinaryReader::read_element_matrices(builder, nodal_solution_mtrx_input_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES); }, [&]()
                     { builder.build_nodal_solution(); });

    reader.close_file_streams(inp_streams);

    const auto &nodal_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::NODAL_SOLUTION);
    auto &nodal_system = eq_sys->get_system<System>(nodal_name);

    nodal_system.solution->close();
    nodal_system.update();
}

void AssemblerLinear::write_element_solution()
{
    apec::logging::info("Computing elemental solution");

    auto eq_sys = this->get_parent_object();

    auto builder = ElementBuilder(this->get_parent_object());

    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = mesh.comm();
    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver();

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"element_sol_mtrx_" + rank_str + ".emat"};
    auto inp_streams = reader.open_input_streams(emat_folderpath, inp_filenames);
    auto &element_solution_mtrx_input_file = *inp_streams[0];

    run_for_elements(mesh, builder, [&]()
                     { BinaryReader::read_element_matrices(builder, element_solution_mtrx_input_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES); }, [&]()
                     { builder.build_element_solution(); });

    reader.close_file_streams(inp_streams);

    const auto &element_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::ELEMENT_SOLUTION);
    auto &element_system = eq_sys->get_system<System>(element_name);

    element_system.solution->close();
    element_system.update();
}

bool AssemblerLinear::check_update_element_matrices()
{
    // Element matrices are recomputed only when the .emat cache files are absent.
    // For linear problems with time-invariant material properties this is correct;
    // nonlinear (ferromagnetic) support will require a proper update check here.
    return false;
}

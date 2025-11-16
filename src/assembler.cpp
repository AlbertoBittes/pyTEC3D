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
#include <libmesh/perf_log.h>

#include <algorithm>
#include <cmath>
#include <filesystem>

#include "assembler.h"
#include "apec_m_system_objects.h"
#include "apec_m_eq_system_object.h"
#include "element_builder.h"
#include "element_types.h"
#include "global_constraints_object.h"
#include "binary_reader.h"

using namespace libMesh;

AssemblerLinear::AssemblerLinear(ApecMEquationSystems *eq_systems,
                                 std::shared_ptr<ContextView> ctx_view)
    : FEMServiceWithOldContext<TimeStep, Session>(ctx_view),
      FEParentObjectBase<ApecMEquationSystems>(eq_systems)
{
}

std::unique_ptr<ElementBuilder> AssemblerLinear::make_element_builder() const
{
    auto eq_sys = this->get_parent_object();
    if (!eq_sys)
        return nullptr;

    auto analysis = eq_sys->get_analysis();
    if (!analysis)
        return nullptr;

    auto analysis_ctx = analysis->get_context_view();
    if (!analysis_ctx)
        return nullptr;

    auto builder_view = analysis_ctx->template project_view<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>();

    return std::make_unique<ElementBuilder>(eq_sys, builder_view);
}

void AssemblerLinear::assemble()
{
    PerfLog perf_log("Matrix Assembly");

    auto eq_sys = this->get_parent_object();
    if (!eq_sys)
        return;

    auto builder = make_element_builder();
    if (!builder)
    {
        libMesh::out << "AssemblerLinear: unable to create element builder for assembly." << std::endl;
        return;
    }

    auto session = this->context<Session>();
    if (!session)
        return;

    auto time_step = this->context<TimeStep>();
    if (!time_step)
        return;

    auto time_substep = time_step->get_active_substep();
    if (!time_substep)
        return;

    auto &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    SparseMatrix<Number> &matrix = system.get_system_matrix();
    system.rhs->zero();

    const DofMap &dof_map = system.get_dof_map();
    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = system.comm();

    const bool update_matrix = this->check_update_system_matrix();
    std::vector<dof_id_type> dof_indices;
    std::vector<Real> old_solution_values;
    std::vector<Real> old_time_derivative_values;

    const Real dt = time_substep->get_dt();
    const Real theta = time_substep->get_theta();
    const std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");
    const std::string rank_str = std::to_string(comm.rank());

    BinaryReader reader;
    auto streams = reader.open_file_streams(emat_folderpath,
                                            {"mtrx_" + rank_str + ".emat",
                                             "rhs_" + rank_str + ".emat"},
                                            FileMode::Read);

    auto *mtrx_input_file = dynamic_cast<std::ifstream *>(streams[0].get());
    auto *rhs_input_file = dynamic_cast<std::ifstream *>(streams[1].get());
    if (!mtrx_input_file || !rhs_input_file)
    {
        reader.close_file_streams(streams);
        libMesh::out << "AssemblerLinear: failed to open element matrix or rhs file." << std::endl;
        return;
    }

    libMesh::out << "Building system global matrix and RHS" << std::endl;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices(elem, dof_indices);

        BinaryReader::read_element_matrices(*builder, *mtrx_input_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
        auto &Ke = builder->get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][0];
        auto &Ce = builder->get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][1];

        Ce.scale(1 / (dt * theta));

        if (update_matrix)
        {
            if (time_substep->get_transient_effect())
                Ke += Ce;

            auto constrained_indices = dof_indices;
            dof_map.constrain_element_matrix(Ke, constrained_indices);
            matrix.add_matrix(Ke, constrained_indices);
        }

        BinaryReader::read_element_vectors(*builder, *rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        auto &rhs = builder->get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];

        if (time_substep->get_transient_effect())
        {
            system.old_local_solution->get(dof_indices, old_solution_values);
            Ce.vector_mult_add(rhs, 1, old_solution_values);

            if (theta != 1.0)
            {
                system.older_local_solution->get(dof_indices, old_time_derivative_values);
                Ce.vector_mult_add(rhs, dt * (1.0 - theta), old_time_derivative_values);
            }
        }

        dof_map.constrain_element_vector(rhs, dof_indices);
        system.rhs->add_vector(rhs, dof_indices);
    }

    reader.close_file_streams(streams);
    libMesh::out << "Finished building system global matrix and RHS" << std::endl;
}
// That's it.  We don't need to do anything else to the
// PerfLog.  When it goes out of scope (at this function return)
// it will print its log to the screen. Pretty easy, huh?

// in a new time step we could have a different dt, different bcs, different time integration effects
// Here add the flag for ferromagnetic
bool AssemblerLinear::check_update_system_matrix()
{
    auto eq_sys = this->get_parent_object();
    if (!eq_sys)
        return false;

    auto &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    GlobalConstraintsObject *bc_object = dynamic_cast<GlobalConstraintsObject *>(&system.get_constraint_object());

    auto current_time_step = this->context<TimeStep>();

    auto current_time_substep = current_time_step ? current_time_step->get_active_substep() : nullptr;

    auto old_time_step = this->old_context<TimeStep>();

    auto old_time_substep = old_time_step ? old_time_step->get_active_substep() : nullptr;

    auto *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());
    const Real rel_eps = 1e-9;

    bool update_matrix = false;

    if (current_time_substep->get_index() == 0)
    {
        // If first substep and BCs have changed, we need to update the matrix
        // This will also include automatically the case of very first substep
        if (bc_object->check_update_boundary_conditions())
            update_matrix = true;

        // If first substep but time integration parameters have changed, we need to update the matrix
        else
        {
            Real diff_dt = std::abs(current_time_substep->get_dt() - old_time_substep->get_dt());
            Real max_old_new_dt = std::max(std::abs(current_time_substep->get_dt()), std::abs(old_time_substep->get_dt()));
            if (diff_dt > rel_eps * max_old_new_dt ||
                current_time_substep->get_transient_effect() != old_time_substep->get_transient_effect() ||
                current_time_substep->get_theta() != old_time_substep->get_theta())
            {
                update_matrix = true;
            }
        }
    }

    if (update_matrix)
    {
        libMesh::out << "Updating system matrix for current time step" << std::endl;
        system.get_system_matrix().zero();
        if (solver)
            solver->reuse_preconditioner(false);
    }
    else
    {
        libMesh::out << "Reusing system matrix from previous time step" << std::endl;
        if (solver)
            solver->reuse_preconditioner(true);
    }

    return update_matrix;
}

void AssemblerLinear::write_element_matrices_and_rhs()
{
    auto eq_sys = this->get_parent_object();
    if (!eq_sys)
        return;

    auto builder = make_element_builder();
    if (!builder)
    {
        libMesh::out << "AssemblerLinear: unable to create element builder for writing element matrices." << std::endl;
        return;
    }

    auto session = this->context<Session>();
    if (!session)
        return;

    auto &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto &dof_map = system.get_dof_map();
    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = system.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    BinaryReader reader;

    std::vector<dof_id_type> dof_indices;

    if (this->check_update_element_matrices())
    {
        libMesh::out << "Recomputing element matrices and vectors for current time step" << std::endl;
        const std::vector<std::string> filenames = {
            "mtrx_" + rank_str + ".emat",
            "nodal_sol_mtrx_" + rank_str + ".emat",
            "element_sol_mtrx_" + rank_str + ".emat",
            "base_rhs_" + rank_str + ".emat"};

        auto streams = reader.open_file_streams(emat_folderpath, filenames, FileMode::Write);
        auto *mtrx_output_file = dynamic_cast<std::ofstream *>(streams[0].get());
        auto *nodal_sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[1].get());
        auto *element_sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[2].get());
        auto *base_rhs_output_file = dynamic_cast<std::ofstream *>(streams[3].get());

        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            dof_map.dof_indices(elem, dof_indices);

            builder->set_element_data(elem);
            builder->build_element();

            BinaryReader::write_element_matrices(*builder, *mtrx_output_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
            BinaryReader::write_element_matrices(*builder, *nodal_sol_mtrx_output_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES);
            BinaryReader::write_element_matrices(*builder, *element_sol_mtrx_output_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES);
            BinaryReader::write_element_vectors(*builder, *base_rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        }

        reader.close_file_streams(streams);
    }
    else
    {
        libMesh::out << "Element matrices will be reused from previous time step" << std::endl;
    }

    const std::vector<std::string> inp_filenames = {"base_rhs_" + rank_str + ".emat"};
    const std::vector<std::string> out_filenames = {"rhs_" + rank_str + ".emat"};

    auto streams_in = reader.open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);
    auto streams_out = reader.open_file_streams(emat_folderpath, out_filenames, FileMode::Write);

    auto *base_rhs_input_file = dynamic_cast<std::ifstream *>(streams_in[0].get());
    auto *rhs_output_file = dynamic_cast<std::ofstream *>(streams_out[0].get());

    libMesh::out << "Applying loads to RHS" << std::endl;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices(elem, dof_indices);

        builder->set_element_data(elem);

        BinaryReader::read_element_vectors(*builder, *base_rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        builder->build_element_rhs();
        BinaryReader::write_element_vectors(*builder, *rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS);
    }

    reader.close_file_streams(streams_in);
    reader.close_file_streams(streams_out);
}

void AssemblerLinear::write_solution_vectors()
{
    libMesh::out << "Computing derived quantities..." << std::endl;

    auto eq_sys = this->get_parent_object();
    if (!eq_sys)
        return;

    auto builder = make_element_builder();
    if (!builder)
    {
        libMesh::out << "AssemblerLinear: unable to create element builder for solution vectors." << std::endl;
        return;
    }

    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = mesh.comm();
    auto session = this->context<Session>();
    if (!session)
        return;

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"nodal_sol_mtrx_" + rank_str + ".emat",
                                                    "element_sol_mtrx_" + rank_str + ".emat"};
    auto inp_streams = reader.open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    auto *nodal_solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());
    auto *element_solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[1].get());

    if (!nodal_solution_mtrx_input_file || !element_solution_mtrx_input_file)
    {
        reader.close_file_streams(inp_streams);
        libMesh::out << "AssemblerLinear: failed to open solution matrix file." << std::endl;
        return;
    }

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder->set_element_data(elem);
        BinaryReader::read_element_matrices(*builder, *nodal_solution_mtrx_input_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES);
        BinaryReader::read_element_matrices(*builder, *element_solution_mtrx_input_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES);
        builder->build_solution_vectors();
    }

    reader.close_file_streams(inp_streams);

    const auto &nodal_name = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, Globals::NODAL_SOLUTION).name;
    const auto &element_name = Globals::get_info_by_key(Globals::SYSTEM_INFO_MAP, Globals::ELEMENT_SOLUTION).name;

    auto &nodal_system = eq_sys->get_system<System>(nodal_name);
    auto &element_system = eq_sys->get_system<System>(element_name);

    nodal_system.solution->close();
    element_system.solution->close();
    nodal_system.update();
    element_system.update();
}

bool AssemblerLinear::check_update_element_matrices()
{
    if (this->context<TimeStep>()->get_active_substep()->get_index() == 0)
    {
        auto old_ts = this->old_context<TimeStep>();
        if (old_ts == nullptr)
            return true;
        else
        {
            // DUMMY for NOW: always return true if old context exists
            // Materials for now cannot change in time, so no need to check that
            return true;
        }
    }
    return false;
}

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

using namespace libMesh;
namespace py = pybind11;

namespace
{
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

    py::dict get_time_step_dict_from_solver(ApecMEquationSystems *eq_sys)
    {
        auto *solver_parent = eq_sys->get_parent_object();
        py::object apec_object = solver_parent->get_apec_object();
        return apec_object.attr("time_step").cast<py::dict>();
    }

    std::filesystem::path get_emat_folderpath_from_solver(ApecMEquationSystems *eq_sys)
    {
        auto *solver_parent = eq_sys->get_parent_object();
        py::object apec_object = solver_parent->get_apec_object();
        const std::string active_mesh_name = py::cast<std::string>(apec_object.attr("active_mesh_name"));
        return std::filesystem::current_path() / active_mesh_name / "emat";
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

    const auto time_step = system.time_step.get();
    const Real dt = time_step->dt;
    const Real theta = py::cast<Real>(time_step->time_step_dict["theta"]);
    const bool transient_effect = py::cast<bool>(time_step->time_step_dict["transient_effects"]);
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver(eq_sys);
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

        // for (int n = 0; n < elem->n_nodes(); n++)
        // {
        //     const auto node_id = elem->node_id(n);
        //     auto node_ptr = elem->node_ptr(n);
        //     dof_id_type dof_id = node_ptr->dof_number(system.number(), 0, 0);
        //     int a = 1;
        // }
        BinaryReader::read_element_matrices(builder, *mtrx_input_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
        auto &Ke = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][0];
        auto &Ce = builder.get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][1];

        // auto check_id = elem->id();
        // auto first_step = (time_step->get_step() == 1 && time_substep->get_index() == 0);
        // if (elem->id() == check_id && mesh.comm().size() > 1 && first_step)
        // {
        //     auto filename = "element_Ke_elem_par_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     Ke.print(file_stream);
        // }
        // else if (elem->id() == check_id && first_step)
        // {
        //     auto filename = "element_Ke_elem_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     Ke.print(file_stream);
        // }
        // if (elem->id() == check_id && mesh.comm().size() > 1 && first_step)
        // {
        //     auto filename = "element_Ce_elem_par_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     Ce.print(file_stream);
        // }
        // else if (elem->id() == check_id && first_step)
        // {
        //     auto filename = "element_Ce_elem_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     Ce.print(file_stream);
        // }

        if (transient_effect)
        {
            Ce.scale(1 / (dt * theta));
            Ke += Ce;
        }
        BinaryReader::read_element_vectors(builder, *rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS);
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

        // if (elem->id() == check_id && mesh.comm().size() > 1 && first_step)
        // {
        //     auto filename = "element_matrix_elem_par_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     Ke.print(file_stream);
        // }
        // else if (elem->id() == check_id && first_step)
        // {
        //     auto filename = "element_matrix_elem_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     Ke.print(file_stream);
        // }

        // if (elem->id() == check_id && mesh.comm().size() > 1 && first_step)
        // {
        //     auto filename = "element_rhs_elem_par_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     rhs.print(file_stream);
        // }
        // else if (elem->id() == check_id && first_step)
        // {
        //     auto filename = "element_rhs_elem_" + std::to_string(elem->id()) + ".txt";
        //     std::ofstream file_stream(filename);
        //     rhs.print(file_stream);
        // }
    }

    reader.close_file_streams(streams);
    libMesh::out << "Finished building system global matrix and RHS" << std::endl;
}
// That's it.  We don't need to do anything else to the
// PerfLog.  When it goes out of scope (at this function return)
// it will print its log to the screen. Pretty easy, huh?

// in a new time step we could have a different dt, different bcs, different time integration effects
// Here add the flag for ferromagnetic
// 1. If nonlinear ferromagnetic, always update
// 2. If new bcs, update
// 3. If new time integration scheme wrt last substep, update
bool AssemblerLinear::check_update_system_matrix()
{
    auto eq_sys = this->get_parent_object();

    auto &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    GlobalConstraintsObject *bc_object = dynamic_cast<GlobalConstraintsObject *>(&system.get_constraint_object());
    auto *solver = dynamic_cast<PetscLinearSolver<Number> *>(system.get_linear_solver());

    if (bc_object->check_update_boundary_conditions())
    {
        libMesh::out << "Updating system matrix for current time step" << std::endl;
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
            const Real current_theta = py::cast<Real>(current_step->time_step_dict["theta"]);
            const Real old_theta = py::cast<Real>(old_step->time_step_dict["theta"]);
            const bool current_transient_effect = py::cast<bool>(current_step->time_step_dict["transient_effects"]);
            const bool old_transient_effect = py::cast<bool>(old_step->time_step_dict["transient_effects"]);

            const bool dt_differs = std::abs(current_dt - old_dt) > std::numeric_limits<Real>::epsilon();
            const bool theta_differs = std::abs(current_theta - old_theta) > std::numeric_limits<Real>::epsilon();

            if ((dt_differs || theta_differs) &&
                current_transient_effect == old_transient_effect)
                update_matrix = true;
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

    auto builder = ElementBuilder(this->get_parent_object());

    auto &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto &dof_map = system.get_dof_map();
    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = system.comm();

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver(eq_sys);

    BinaryReader reader;

    std::vector<dof_id_type> dof_indices;

    const std::vector<std::string> filenames = {
        "mtrx_" + rank_str + ".emat",
        "nodal_sol_mtrx_" + rank_str + ".emat",
        "element_sol_mtrx_" + rank_str + ".emat",
        "base_rhs_" + rank_str + ".emat"};

    const bool emat_files_exist = std::all_of(filenames.begin(), filenames.end(), [&](const std::string &filename)
                                              { return std::filesystem::exists(emat_folderpath / filename); });

    if (!emat_files_exist)
    {
        libMesh::out << "Recomputing element matrices and vectors because one or more .emat files are missing" << std::endl;

        auto streams = reader.open_file_streams(emat_folderpath, filenames, FileMode::Write);
        auto *mtrx_output_file = dynamic_cast<std::ofstream *>(streams[0].get());
        auto *nodal_sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[1].get());
        auto *element_sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[2].get());
        auto *base_rhs_output_file = dynamic_cast<std::ofstream *>(streams[3].get());

        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            builder.set_element_data(elem);

            builder.build_element();

            BinaryReader::write_element_matrices(builder, *mtrx_output_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
            BinaryReader::write_element_matrices(builder, *nodal_sol_mtrx_output_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES);
            BinaryReader::write_element_matrices(builder, *element_sol_mtrx_output_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES);
            BinaryReader::write_element_vectors(builder, *base_rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS);
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

    builder.clear_nodal_loads_dof_cache();

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);

        BinaryReader::read_element_vectors(builder, *base_rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        builder.build_element_rhs();
        BinaryReader::write_element_vectors(builder, *rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS);
    }

    reader.close_file_streams(streams_in);
    reader.close_file_streams(streams_out);
}

void AssemblerLinear::write_nodal_solution()
{
    libMesh::out << "Computing derived quantities..." << std::endl;

    auto eq_sys = this->get_parent_object();

    auto builder = ElementBuilder(this->get_parent_object());

    builder.set_compute_nodal_solution_vectors(true);

    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = mesh.comm();
    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver(eq_sys);

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"nodal_sol_mtrx_" + rank_str + ".emat"};
    auto inp_streams = reader.open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    auto *nodal_solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());

    if (!nodal_solution_mtrx_input_file)
    {
        reader.close_file_streams(inp_streams);
        libMesh::out << "AssemblerLinear: failed to open solution matrix file." << std::endl;
        return;
    }

    run_for_elements(mesh, builder, [&]()
                     { BinaryReader::read_element_matrices(builder, *nodal_solution_mtrx_input_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES); }, [&]()
                     { builder.build_nodal_solution(); });

    reader.close_file_streams(inp_streams);

    const auto &nodal_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::NODAL_SOLUTION);
    auto &nodal_system = eq_sys->get_system<System>(nodal_name);

    nodal_system.solution->close();
    nodal_system.update();
}

void AssemblerLinear::write_element_solution()
{
    libMesh::out << "Computing derived quantities..." << std::endl;

    auto eq_sys = this->get_parent_object();

    auto builder = ElementBuilder(this->get_parent_object());

    builder.set_compute_nodal_solution_vectors(false);

    const MeshBase &mesh = eq_sys->get_mesh();
    const auto &comm = mesh.comm();
    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = get_emat_folderpath_from_solver(eq_sys);

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"element_sol_mtrx_" + rank_str + ".emat"};
    auto inp_streams = reader.open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    auto *element_solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());

    if (!element_solution_mtrx_input_file)
    {
        reader.close_file_streams(inp_streams);
        libMesh::out << "AssemblerLinear: failed to open solution matrix file." << std::endl;
        return;
    }

    run_for_elements(mesh, builder, [&]()
                     { BinaryReader::read_element_matrices(builder, *element_solution_mtrx_input_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES); }, [&]()
                     { builder.build_element_solution(); });

    reader.close_file_streams(inp_streams);

    const auto &element_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(ApecMSystemBase::ELEMENT_SOLUTION);
    auto &element_system = eq_sys->get_system<System>(element_name);

    element_system.solution->close();
    element_system.update();
}

bool AssemblerLinear::check_update_element_matrices()
{
    // auto ts = this->context<TimeStep>();
    // // Assuming only linear properties for now
    // // You must recompute element matrices only if it is the first substep of the time step, first PIecard operation
    // // No need to recompute if not first Picard iteration or first substep, since for linear problems the element matrices do not change in time
    // // If first substep of the time step
    // if (ts->get_active_substep()->get_index() == 0)
    // {
    //     auto old_ts = ts->get_previous_time_step();
    //     // If not first picard iteration, no need to recompute element matrices
    //     if (ts->get_active_substep()->get_number_picard_iterations() > 0)
    //         return false;
    //     // If very first substep, need to recompute element matrices
    //     if (old_ts == nullptr)
    //         return true;
    //     // else
    //     // {
    //     //     // DUMMY for NOW: always return true if old context exists
    //     //     // Materials for now cannot change in time, so no need to check that
    //     //     return true;
    //     // }
    // }
    // return false;

    // Change: compute matrices if not present as file.
    // In the future: check for nonlinear elements, dummy for now
    return false;
}

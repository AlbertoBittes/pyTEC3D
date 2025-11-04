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

void AssemblerLinear::assemble()
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    PerfLog perf_log("Matrix Assembly");
    // Get a constant reference to the mesh object.
    // Get a reference to the mesh object.
    TransientSystem<LinearImplicitSystem> *system = dynamic_cast<TransientSystem<LinearImplicitSystem> *>(this->get_parent_object());

    ApecMEquationSystems *eq_sys = dynamic_cast<ApecMEquationSystems *>(&system->get_equation_systems());

    const MeshBase &mesh = system->get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving

    const auto &comm = system->comm();
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap &dof_map = parent_object->get_dof_map();

    auto session = this->context<Session>();

    std::shared_ptr<ElementBuilder> elem_builder = eq_sys->get_element_builder();
    if (!elem_builder)
        return;

    auto time_substep = this->context<TimeStep>()->get_active_substep();

    Real dt = time_substep->get_dt();
    Real theta = time_substep->get_theta();
    // if (this->check_compute_matrix())
    // {
    //     this->build_write_element_matrices();
    // }
    // else if (time_substep->get_index() == 0)
    // {
    //     this->update_rhs();
    // }
    std::vector<dof_id_type> dof_indices;
    std::vector<Real> old_solution_values;
    std::vector<Real> old_time_derivative_values;

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.

    // The global system matrix
    SparseMatrix<Number> &matrix = system->get_system_matrix();
    bool update_matrix = check_update_matrix(time_substep);
    unsigned int n_dofs;

    std::filesystem::path filepath;
    std::string mtrx_vec_filename;

    auto emat_folderpath = session->get_path("emat_folderpath");
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Here add the flag for ferromagnetic

    const auto rank_str = std::to_string(comm.rank());

    std::vector<std::string> filenames = {
        "mtrx_" + rank_str + ".temp",
        "rhs_" + rank_str + ".temp"};

    // Open the main files
    BinaryReader reader;
    auto streams = reader.open_file_streams(emat_folderpath, filenames, FileMode::Read);

    // Downcast to use as ifstream
    std::ifstream *mtrx_input_file = dynamic_cast<std::ifstream *>(streams[0].get());
    std::ifstream *rhs_input_file = dynamic_cast<std::ifstream *>(streams[1].get());

    // Matrix changes when we have nonlinearities or different dt

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations
    // Also, note that here we call heterogenously_constrain_element_matrix_and_vector
    // to impose a inhomogeneous Dirichlet boundary conditions.

    // The element matrix and right-hand-side are now built
    // for this element.  Add them to the global matrix and
    // right-hand-side vector.  The SparseMatrix::add_matrix()
    // and NumericVector::add_vector() members do this for us.
    // Start logging the insertion of the local (element)
    // matrix and vector into the global matrix and vector
    // LOG_SCOPE_WITH("matrix insertion", "", perf_log);
    libMesh::out << "Building system global matrix and RHS" << std::endl;
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // std::shared_ptr<ElementType> elem_type = element_types_map[elem->get_extra_integer(etype_num_integer)];
        // // Start logging the shape function initialization.
        // // This is done through a simple function call with
        // // the name of the event to log.
        // conductive = elem_type->is_conductive();
        // auto elem_id = elem->id();
        // this can be done at the beginning only
        dof_map.dof_indices(elem, dof_indices);
        // Cache the number of degrees of freedom on this element, for
        // use as a loop bound later.  We use cast_int to explicitly
        // convert from size() (which may be 64-bit) to unsigned int
        // (which may be 32-bit but which is definitely enough to count
        // *local* degrees of freedom.
        n_dofs = cast_int<unsigned int>(dof_indices.size());
        // Read the matrices and vector from file

        // if (elem_id == 0)
        // {
        //     // libMesh::out << "Time step changed from " << dt2 << " to " << dt1 << std::endl;
        //     libMesh::out << "Computing new system matrix" << std::endl;
        // }

        BinaryReader::read_element_matrices(*elem_builder, *mtrx_input_file, 0);

        auto &Ke = elem_builder->get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][0];
        auto &Ce = elem_builder->get_matrices()[ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES][1];

        Ce.scale(1 / (time_substep->get_dt() * theta));

        if (update_matrix)
        {
            if (time_substep->get_transient_effect())
                Ke += Ce;

            auto dof_indices_buf = dof_indices;
            dof_map.constrain_element_matrix(Ke, dof_indices_buf);
            matrix.add_matrix(Ke, dof_indices_buf);
        }

        BinaryReader::read_element_vectors(*elem_builder, *rhs_input_file, 0);
        auto &rhs = elem_builder->get_vectors()[ElementBuilder::VectorsSet::RHS_VECTORS][0];

        // elem_builder->apply_load();

        // HERE YOU HAVE TO RAMP THE LOAD !!!

        // dof_map.constrain_element_vector(FFe, dof_indices);
        // system.rhs->add_vector(FFe, dof_indices);
        // here I constrain but I modify later, wrong
        // }
        if (time_substep->get_transient_effect())
        {
            system->old_local_solution->get(dof_indices, old_solution_values);

            Ce.vector_mult_add(rhs, 1, old_solution_values);

            if (theta != 1.0)
            {
                system->older_local_solution->get(dof_indices, old_time_derivative_values);
                Ce.vector_mult_add(rhs, dt * (1 - theta), old_time_derivative_values);
            }
            // auto elem_id = elem->id();
            // if (elem_id == 169)
            // {
            //     std::ofstream ffe_file("rhs_parallel_" + std::to_string(elem_id) + "_" + std::to_string(time_substep->get_index()) + ".txt");
            //     if (ffe_file.is_open())
            //     {
            //         rhs.print_scientific(ffe_file);
            //     }
            //     ffe_file.close();
            //     std::ofstream Ke_file("matrix_parallel_" + std::to_string(elem_id) + "_" + std::to_string(time_substep->get_index()) + ".txt");
            //     if (Ke_file.is_open())
            //     {
            //         Ke.print_scientific(Ke_file);
            //     }
            //     Ke_file.close();
            // }
        }

        dof_map.constrain_element_vector(rhs, dof_indices);
        system->rhs->add_vector(rhs, dof_indices);
    }
    reader.close_file_streams(streams);
    this->prev_substep = time_substep;
    libMesh::out << "Finished building system global matrix and RHS" << std::endl;
    libMesh::out << "Solving linear system..." << std::endl;
    // Matrix changes when we have nonlinearities or different dt

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).

    // If this assembly program were to be used on an adaptive mesh,
    // we would have to apply any hanging node constraint equations
    // Also, note that here we call heterogenously_constrain_element_matrix_and_vector
    // to impose a inhomogeneous Dirichlet boundary conditions.

    // The element matrix and right-hand-side are now built
    // for this element.  Add them to the global matrix and
    // right-hand-side vector.  The SparseMatrix::add_matrix()
    // and NumericVector::add_vector() members do this for us.
    // Start logging the insertion of the local (element)
    // matrix and vector into the global matrix and vector
    // LOG_SCOPE_WITH("matrix insertion", "", perf_log);
}
// That's it.  We don't need to do anything else to the
// PerfLog.  When it goes out of scope (at this function return)
// it will print its log to the screen. Pretty easy, huh?

// in a new time step we could have a different dt, different bcs, different time integration effects
// Here add the flag for ferromagnetic
bool AssemblerLinear::check_update_matrix(Substep *time_substep)
{
    bool update_matrix = false;

    auto *system = dynamic_cast<TransientSystem<LinearImplicitSystem> *>(this->get_parent_object());
    if (!system)
        return false;

    auto *solver = dynamic_cast<PetscLinearSolver<Number> *>(system->get_linear_solver());
    const Real rel_eps = 1e-9;

    if (prev_substep == nullptr)
    {
        update_matrix = true;
        if (solver)
            solver->reuse_preconditioner(false);
    }
    else if (time_substep->get_index() > 0)
    {
        update_matrix = false;
        if (solver)
            solver->reuse_preconditioner(true);
    }
    else
    {
        Real diff_dt = std::abs(time_substep->get_dt() - prev_substep->get_dt());
        Real max_old_new_dt = std::max(std::abs(time_substep->get_dt()), std::abs(prev_substep->get_dt()));
        if (diff_dt > rel_eps * max_old_new_dt ||
            time_substep->get_transient_effect() != prev_substep->get_transient_effect() ||
            time_substep->get_theta() != prev_substep->get_theta())
        {
            update_matrix = true;
            if (solver)
                solver->reuse_preconditioner(false);
        }
        else
        {
            update_matrix = false;
            if (solver)
                solver->reuse_preconditioner(true);
        }
    }

    auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&system->get_equation_systems());
    if (auto builder = equation_systems ? equation_systems->get_element_builder() : nullptr)
        builder->update_matrix = update_matrix;

    if (update_matrix)
        libMesh::out << "Updating system matrix for current time step" << std::endl;
    else
        libMesh::out << "Reusing system matrix from previous time step" << std::endl;

    prev_substep = time_substep;
    return update_matrix;
}

void AssemblerLinear::write_element_matrices_and_rhs(ElementBuilder &builder)
{
    auto equation_systems = builder.get_parent_object();
    if (!equation_systems)
        return;

    const MeshBase &mesh = equation_systems->get_mesh();
    const auto &comm = equation_systems->comm();
    auto session = this->context<Session>();
    if (!session)
        return;

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    BinaryReader reader;

    if (builder.needs_matrix_update())
    {
        libMesh::out << "Recomputing element matrices and vectors for current time step" << std::endl;
        const std::vector<std::string> filenames = {
            "mtrx_" + rank_str + ".temp",
            "sol_mtrx_" + rank_str + ".temp",
            "base_rhs_" + rank_str + ".temp"};

        auto streams = reader.open_file_streams(emat_folderpath, filenames, FileMode::Write);
        auto *mtrx_output_file = dynamic_cast<std::ofstream *>(streams[0].get());
        auto *sol_mtrx_output_file = dynamic_cast<std::ofstream *>(streams[1].get());
        auto *base_rhs_output_file = dynamic_cast<std::ofstream *>(streams[2].get());

        auto *main_dof_map = builder.dof_maps[Globals::DOF_SOLUTION_SYSTEM];

        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            if (main_dof_map)
                main_dof_map->dof_indices(elem, builder.dof_indices);

            builder.set_element_data(elem);
            builder.build_element();

            BinaryReader::write_element_matrices(builder, *mtrx_output_file, ElementBuilder::MatricesSet::DOF_SOLUTION_MATRICES);
            BinaryReader::write_element_matrices(builder, *sol_mtrx_output_file, ElementBuilder::MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES);
            BinaryReader::write_element_vectors(builder, *base_rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        }

        reader.close_file_streams(streams);
    }
    else
    {
        libMesh::out << "Element matrices will be reused from previous time step" << std::endl;
    }

    const std::vector<std::string> inp_filenames = {"base_rhs_" + rank_str + ".temp"};
    const std::vector<std::string> out_filenames = {"rhs_" + rank_str + ".temp"};

    auto streams_in = reader.open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);
    auto streams_out = reader.open_file_streams(emat_folderpath, out_filenames, FileMode::Write);

    auto *base_rhs_input_file = dynamic_cast<std::ifstream *>(streams_in[0].get());
    auto *rhs_output_file = dynamic_cast<std::ofstream *>(streams_out[0].get());

    libMesh::out << "Applying loads to RHS" << std::endl;

    auto *main_dof_map = builder.dof_maps[Globals::DOF_SOLUTION_SYSTEM];

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        if (main_dof_map)
            main_dof_map->dof_indices(elem, builder.dof_indices);

        builder.set_element_data(elem);

        BinaryReader::read_element_vectors(builder, *base_rhs_input_file, ElementBuilder::VectorsSet::RHS_VECTORS);
        builder.build_element_rhs();
        BinaryReader::write_element_vectors(builder, *rhs_output_file, ElementBuilder::VectorsSet::RHS_VECTORS);
    }

    reader.close_file_streams(streams_in);
    reader.close_file_streams(streams_out);
}

void AssemblerLinear::write_solution_vectors(ElementBuilder &builder)
{
    auto equation_systems = builder.get_parent_object();
    if (!equation_systems)
        return;

    const MeshBase &mesh = equation_systems->get_mesh();
    const auto &comm = mesh.comm();
    auto session = this->context<Session>();
    if (!session)
        return;

    const std::string rank_str = std::to_string(comm.rank());
    const std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    BinaryReader reader;
    const std::vector<std::string> inp_filenames = {"sol_mtrx_" + rank_str + ".temp"};
    auto inp_streams = reader.open_file_streams(emat_folderpath, inp_filenames, FileMode::Read);

    auto *solution_mtrx_input_file = dynamic_cast<std::ifstream *>(inp_streams[0].get());

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        builder.set_element_data(elem);
        BinaryReader::read_element_matrices(builder, *solution_mtrx_input_file, ElementBuilder::MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES);
        builder.build_solution_vectors();
    }

    solution_mtrx_input_file->close();
    reader.close_file_streams(inp_streams);

    builder.systems[Globals::NODAL_SOLUTION]->solution->close();
    builder.systems[Globals::ELEMENT_SOLUTION]->solution->close();
    builder.systems[Globals::NODAL_SOLUTION]->update();
    builder.systems[Globals::ELEMENT_SOLUTION]->update();
}

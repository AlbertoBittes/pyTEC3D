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

#include "assembler.h"
#include "apec_m_system_objects.h"
#include "element_builder.h"
#include "element_types.h"
#include "system_builder.h"

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
        "rhs_" + rank_str + ".temp",
        "rhs_" + std::to_string(comm.rank()) + "_prev.temp"};

    // Open the main files
    auto streams = session->open_file_streams(emat_folderpath, filenames, FileMode::Read);

    // Downcast to use as ifstream
    std::ifstream *mtrx_input_file = dynamic_cast<std::ifstream *>(streams[0].get());
    std::ifstream *rhs_input_file = dynamic_cast<std::ifstream *>(streams[1].get());
    std::ifstream *rhs_prev_input_file = dynamic_cast<std::ifstream *>(streams[2].get());

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
    auto &vecs = elem_builder->get_vectors()[1];
    vecs.resize(1);
    auto &dummy_prev_rhs = vecs[0];
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

        elem_builder->read_element_matrices(*mtrx_input_file, 0);

        auto &Ke = elem_builder->get_matrices()[0][0];
        auto &Ce = elem_builder->get_matrices()[0][1];

        Ce.scale(1 / (time_substep->get_dt() * theta));

        if (update_matrix)
        {
            if (time_substep->get_transient_effect())
                Ke += Ce;

            dof_map.constrain_element_matrix(Ke, dof_indices);
            matrix.add_matrix(Ke, dof_indices);
        }

        elem_builder->read_element_vectors(*rhs_input_file, 0);
        auto &rhs = elem_builder->get_vectors()[0][0];

        if (time_substep->get_ramped_load())
        {
            if (rhs_prev_input_file == nullptr)
            {
                dummy_prev_rhs.resize(n_dofs);
            }
            else
            {
                elem_builder->read_element_vectors(*rhs_prev_input_file, 1);
            }
            elem_builder->ramp_vectors(0, 1, 0, 0);
        }

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
            // if (elem_id == 969)
            // {
            //     std::ofstream ffe_file("../simulations/two_wires_3D/matrix/FFe_matrix" + std::to_string(time_step.get_step()) + ".txt");
            //     if (ffe_file.is_open())
            //     {
            //         FFe.print_scientific(ffe_file);
            //     }
            // }
        }

        dof_map.constrain_element_vector(rhs, dof_indices);
        system->rhs->add_vector(rhs, dof_indices);
    }
    session->close_file_streams(streams);
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

    ApecMLinearSystem *system = dynamic_cast<ApecMLinearSystem *>(this->get_parent_object());

    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(system->get_linear_solver());

    GlobalConstraintsObject *bc_object = dynamic_cast<GlobalConstraintsObject *>(&system->get_constraint_object());

    const Real rel_eps = 1e-9;

    if (prev_substep == nullptr)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Or if element matrices have changed
    else
    {

        if (time_substep->get_index() > 0)
        {
            update_matrix = false;
            solver->reuse_preconditioner(true);
        }
        else
        {
            Real diff_dt = std::abs(time_substep->get_dt() - prev_substep->get_dt());
            Real max_old_new_dt = std::max(std::abs(time_substep->get_dt()), std::abs(prev_substep->get_dt()));
            if (diff_dt > rel_eps * max_old_new_dt ||
                bc_object->get_update_boundary_conditions() ||
                time_substep->get_transient_effect() != prev_substep->get_transient_effect())
            {
                update_matrix = true;
                solver->reuse_preconditioner(false);
            }
            else
            {
                update_matrix = false;
                solver->reuse_preconditioner(true);
            }
        }
    }
    if (update_matrix)
        libMesh::out << "Updating system matrix for current time step" << std::endl;
    else
        libMesh::out << "Reusing system matrix from previous time step" << std::endl;
    return update_matrix;
}

// This you will need to improve it by making it a method of equations systems to save time on the resizing and retrieving of solution
void AssemblerSolution::compute_derived_quantities()
{
    ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

    ApecMEquationSystems *eq_sys = dynamic_cast<ApecMEquationSystems *>(&solution_system->get_equation_systems());

    std::shared_ptr<ElementBuilder> elem_builder = eq_sys->get_element_builder();
    // ApecMLinearSystem &system = eq_sys->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    const MeshBase &mesh = solution_system->get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = solution_system->comm();

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    // const DofMap &dof_map = system.get_dof_map();
    const DofMap &solution_dof_map = solution_system->get_dof_map();

    auto session = this->context<Session>();

    // std::shared_ptr<ElementBuilder> elem_builder = eq_sys->get_element_builder();

    auto sys_name = solution_system->name();

    const auto rank_str = std::to_string(comm.rank());

    std::filesystem::path emat_folderpath = session->get_path("emat_folderpath");

    std::vector<std::string> filenames = {};

    unsigned int vec_set = 0;
    if (sys_name == "NodalSolutionSystem")
    {
        filenames.push_back("nodal_sol_vec_" + rank_str + ".temp");
        vec_set = 2;
        libMesh::out << "Setting nodal solution results " << std::endl;
    }
    else if (sys_name == "ElementSolutionSystem")
    {
        filenames.push_back("elem_sol_vec_" + rank_str + ".temp");
        vec_set = 3;
        libMesh::out << "Setting element solution results " << std::endl;
    }
    auto streams = session->open_file_streams(emat_folderpath, filenames, FileMode::Read);
    std::ifstream *solution_input_file = dynamic_cast<std::ifstream *>(streams[0].get());

    std::vector<dof_id_type> dof_indices;

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        if (elem->infinite())
            continue;

        solution_dof_map.dof_indices(elem, dof_indices);

        elem_builder->read_element_vectors(*solution_input_file, vec_set);

        auto &solution_vector = elem_builder->get_vectors()[vec_set][0];

        for (unsigned int i = 0; i < dof_indices.size(); i++)
        {
            dof_id_type dof_index = dof_indices[i];
            if ((solution_system->solution->first_local_index() <= dof_index) &&
                (dof_index < solution_system->solution->last_local_index()))
                solution_system->solution->set(dof_index, solution_vector(i));
        }
    }
    solution_system->solution->close();
    solution_system->update();
}
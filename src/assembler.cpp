#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "session.h" // Include the header file for Session
#include "apec_m_linear_system.h"
#include "libmesh/string_to_enum.h"
#include "fem_data_context.h"
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

#include "libmesh/petsc_linear_solver.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include "libmesh/enum_solver_package.h"
#include "assembler.h"
#include <map> // Include for std::map
#include <cmath>
#include <stdexcept>

using namespace libMesh;

// // We now define the matrix assembly function for the
// // Poisson system.  We need to first compute element
// // matrices and right-hand sides, and then take into
// // account the boundary conditions.
void AssemblerLinear::assemble()
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    PerfLog perf_log("Matrix Assembly");
    // Get a constant reference to the mesh object.
    // Get a reference to the mesh object.
    ApecMLinearSystem *system = dynamic_cast<ApecMLinearSystem *>(this->get_parent_object());

    const MeshBase &mesh = system->get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving

    const auto &comm = system->comm();
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap &dof_map = parent_object->get_dof_map();

    auto context = std::dynamic_pointer_cast<AssemblerContext>(this->get_context());

    auto session = context->get_session();

    ElementBuilder *elem_builder = system->get_element_builder();

    auto time_substep = context->get_time_step()->get_active_substep();

    if (this->check_compute_matrix())
    {
        this->build_write_element_matrices();
    }
    else if (time_substep->get_index() == 0)
    {
        this->update_rhs();
    }
    std::vector<dof_id_type> dof_indices;
    std::vector<Real> old_solution_values;
    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.

    // The global system matrix
    SparseMatrix<Number> &matrix = system->get_system_matrix();
    bool update_matrix = check_update_matrix();
    unsigned int n_dofs;

    std::filesystem::path filepath;
    std::string mtrx_vec_filename;

    auto emat_folderpath = session->get_emat_folderpath();
    bool old_rhs = true;
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Here add the flag for ferromagnetic

    mtrx_vec_filename = "mtrx_" + std::to_string(comm.rank()) + ".temp";
    filepath = emat_folderpath / mtrx_vec_filename;
    std::ifstream mtrx_input_file(filepath);
    if (!mtrx_input_file.is_open())
    {
        libmesh_error_msg("Failed to open element data file for reading.");
    }

    mtrx_vec_filename = "rhs_" + std::to_string(comm.rank()) + ".temp";
    filepath = emat_folderpath / mtrx_vec_filename;
    std::ifstream rhs_input_file(filepath);
    if (!rhs_input_file.is_open())
    {
        libmesh_error_msg("Failed to open element data file for reading.");
    }

    std::ifstream rhs_prev_input_file;
    if (time_substep->get_parent_object()->get_step() == 1)
    {
        old_rhs = false;
    }
    else
    {
        mtrx_vec_filename = "rhs_" + std::to_string(comm.rank()) + "_prev.temp";
        filepath = emat_folderpath / mtrx_vec_filename;
        rhs_prev_input_file.open(filepath);
        if (!rhs_prev_input_file.is_open())
        {
            libmesh_error_msg("Failed to open element data file for reading.");
        }
    }

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

        elem_builder->read_element_matrices(mtrx_input_file, 0);

        auto Ke = elem_builder->get_matrices()[0][0];
        auto Ce = elem_builder->get_matrices()[0][1];

        Ce.scale(1 / time_substep->get_dt());

        if (update_matrix)
        {
            if (time_substep->get_transient_effect())
                Ke += Ce;

            dof_map.constrain_element_matrix(Ke, dof_indices);
            matrix.add_matrix(Ke, dof_indices);
        }

        elem_builder->read_element_vectors(rhs_input_file, 0);
        auto rhs = elem_builder->get_vectors()[0][0];

        if (time_substep->get_ramped_load())
        {
            if (old_rhs)
                elem_builder->read_element_vectors(rhs_prev_input_file, 1);
            else
            {
                elem_builder->get_vectors()[0][1].resize(rhs.size());
                elem_builder->get_vectors()[0][1].zero();
            }
            elem_builder->ramp_vectors(time_substep);
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
    mtrx_input_file.close();
    rhs_input_file.close();
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

void AssemblerLinear::build_write_element_matrices()
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    PerfLog perf_log("Matrix Assembly");

    // Get a constant reference to the mesh object.
    const MeshBase &mesh = eq_sys->get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = eq_sys->comm();
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.

    auto session = eq_sys->parameters.get<std::shared_ptr<Session>>("session");

    auto element_types_collector = eq_sys->parameters.get<std::shared_ptr<ElementTypesCollector>>("element_types_collector");
    ElementTypesMap &element_types_map = element_types_collector->elementTypes;

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.

    std::filesystem::path filepath;
    std::string mtrx_vec_filename;

    unsigned int etype_num_index = mesh.get_elem_integer_index("ansys_etype");
    unsigned int mat_num_index = mesh.get_elem_integer_index("ansys_material_type");
    unsigned int real_constants_num_index = mesh.get_elem_integer_index("ansys_real_constant");

    dof_id_type etype_num;
    dof_id_type material_num;
    dof_id_type real_constant_num;

    auto emat_folderpath = session->get_emat_folderpath();

    // if (old_time_substep == nullptr) // Here add the check for nonlinear elements
    // {
    // This function must be called prior to a check for the need of recomputing the matrix:
    // e.g. new materials, new real constants, nonlinear elements.

    mtrx_vec_filename = "mtrx_" + std::to_string(comm.rank()) + ".temp";
    filepath = emat_folderpath / mtrx_vec_filename;
    std::ofstream mtrx_output_file(filepath);
    if (!mtrx_output_file.is_open())
    {
        libmesh_error_msg("Failed to open element data file for writing.");
    }

    mtrx_vec_filename = "der_mtrx_" + std::to_string(comm.rank()) + ".temp";
    filepath = emat_folderpath / mtrx_vec_filename;
    std::ofstream der_mtrx_output_file(filepath);

    if (!der_mtrx_output_file.is_open())
    {
        libmesh_error_msg("Failed to open element data file for writing.");
    }

    mtrx_vec_filename = "rhs_" + std::to_string(comm.rank()) + ".temp";
    filepath = emat_folderpath / mtrx_vec_filename;
    if (std::filesystem::exists(filepath))
    {
        std::filesystem::path new_filepath = emat_folderpath / ("rhs_" + std::to_string(comm.rank()) + "_prev.temp");
        std::filesystem::rename(filepath, new_filepath);
    }
    std::ofstream rhs_output_file(filepath);
    if (!rhs_output_file.is_open())
    {
        libmesh_error_msg("Failed to open element data file for writing.");
    }

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.

        etype_num = elem->get_extra_integer(etype_num_index);
        material_num = elem->get_extra_integer(mat_num_index);
        real_constant_num = elem->get_extra_integer(real_constants_num_index);

        std::shared_ptr<ElementType> elem_type = element_types_map[etype_num];

        // Do this elsewhere, logically must no stay here
        // load_step->get_materials()->compute_material_properties(elem->id(), conductive);
        elem_builder->set_element_data(elem, elem_type, material_num, real_constant_num);

        elem_builder->prepare_element_matrix_building();

        elem_builder->build_element_matrix();

        elem_builder->build_element_rhs();
        // Here the function to write the matrices and vectors down to file
        elem_builder->writeElementMatrices(mtrx_output_file, elem_builder->get_matrices());

        elem_builder->writeElementMatrices(der_mtrx_output_file, elem_builder->get_matrices_derived());

        elem_builder->writeElementVectors(rhs_output_file, elem_builder->get_vectors());
    }
    mtrx_output_file.close();
    der_mtrx_output_file.close();
    rhs_output_file.close();
}

void AssemblerLinear::update_rhs()
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    PerfLog perf_log("Matrix Assembly");
    // Get a constant reference to the mesh object.
    auto eq_sys = this->get_parent_object();

    const MeshBase &mesh = eq_sys->get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    const auto &comm = eq_sys->comm();
    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.

    auto session = eq_sys->parameters.get<std::shared_ptr<Session>>("session");

    auto element_types_collector = eq_sys->parameters.get<std::shared_ptr<ElementTypesCollector>>("element_types_collector");
    ElementTypesMap &element_types_map = element_types_collector->elementTypes;

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.

    std::filesystem::path filepath;
    std::string mtrx_vec_filename;

    unsigned int etype_num_index = mesh.get_elem_integer_index("ansys_etype");

    unsigned int etype_num;

    auto emat_folderpath = session->get_emat_folderpath();

    // if (old_time_substep == nullptr) // Here add the check for nonlinear elements
    // {
    // This function must be called prior to a check for the need of recomputing the matrix:
    // e.g. new materials, new real constants, nonlinear elements.

    mtrx_vec_filename = "rhs_" + std::to_string(comm.rank()) + ".temp";
    filepath = emat_folderpath / mtrx_vec_filename;
    if (std::filesystem::exists(filepath))
    {
        std::filesystem::path new_filepath = emat_folderpath / ("rhs_" + std::to_string(comm.rank()) + "_prev.temp");
        std::filesystem::rename(filepath, new_filepath);
    }
    std::ofstream rhs_output_file(filepath);
    if (!rhs_output_file.is_open())
    {
        libmesh_error_msg("Failed to open element data file for writing.");
    }

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.

        etype_num = elem->get_extra_integer(etype_num_index);

        std::shared_ptr<ElementType> elem_type = element_types_map[etype_num];

        // Do this elsewhere, logically must no stay here
        // load_step->get_materials()->compute_material_properties(elem->id(), conductive);
        elem_builder->set_element_type(elem_type);

        elem_builder->prepare_element_rhs_building();

        elem_builder->build_element_rhs();
        // Here the function to write the matrices and vectors down to file

        elem_builder->writeElementVectors(rhs_output_file, elem_builder->get_vectors());
    }
    rhs_output_file.close();
}
// in a new time step we could have a different dt, different bcs, different time integration effects
// Here add the flag for ferromagnetic
bool AssemblerLinear::check_update_matrix(TimeSubstep *old_time_substep)
{
    bool update_matrix = false;

    const Real relative_epsilon = 1e-9;
    auto eq_sys = this->get_parent_object();
    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(eq_sys->get_system<TransientLinearImplicitSystem>("ApecMLinearSystem").get_linear_solver());
    auto load_step = eq_sys->parameters.get<std::shared_ptr<LoadStepper>>("load_stepper")->get_active_load_step();
    auto prev_load_step = load_step->get_parent_object()->get_old_load_step();
    auto time_step = load_step->get_time_step();
    auto active_time_substep = time_step->get_active_substep();
    auto old_time_substep = time_step->get_previous_substep();

    if (old_time_substep == nullptr)
    {
        if (prev_load_step != nullptr)
        {
            old_time_substep = &(prev_load_step->get_time_step()->get_substeps().back());
        }
        else
        {
            old_time_substep = nullptr;
        }
    }
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Here add the flag for ferromagnetic
    if (old_time_substep == nullptr)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    else
    {

        if (active_time_substep->get_index() > 0)
        {
            update_matrix = false;
            solver->reuse_preconditioner(true);
        }
        else
        {
            Real diff_dt = std::abs(active_time_substep->get_dt() - old_time_substep->get_dt());
            Real max_old_new_dt = std::max(std::abs(active_time_substep->get_dt()), std::abs(old_time_substep->get_dt()));
            if (diff_dt > rel_eps * max_old_new_dt ||
                node_dof_constraints->get_update_boundary_conditions() ||
                active_time_substep->get_transient_effect() != old_time_substep->get_transient_effect())
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
    return update_matrix;
}

bool AssemblerLinear::check_compute_matrix()
{
    // bool compute_matrix = false;
    auto eq_sys = this->get_parent_object();
    auto load_step = eq_sys->parameters.get<std::shared_ptr<LoadStepper>>("load_stepper")->get_active_load_step();
    auto prev_load_step = load_step->get_parent_object()->get_old_load_step();
    auto active_time_substep = load_step->get_time_step()->get_active_substep();
    auto old_time_substep = load_step->get_time_step()->get_previous_substep();
    bool update_matrix;
    auto solver = eq_sys->get_system<TransientLinearImplicitSystem>("ApecMLinearSystem").get_linear_solver();
    Real relative_epsilon = 1e-8;
    if (prev_load_step == nullptr && active_time_substep->get_index() == 0)
    {
        // Here add the check for nonlinear elements
        // This function must be called prior to a check for the need of recomputing the matrix:
        // e.g. new materials, new real constants, nonlinear elements.
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    else if (old_time_substep == nullptr && prev_load_step != nullptr)
    {
        old_time_substep = &(prev_load_step->get_time_step()->get_substeps().back());
    }
    // else if (active_time_substep == nullptr && prev_load_step != nullptr)
    // {
    //     // Here add the check for nonlinear elements
    //     // This function must be called prior to a check for the need of recomputing the matrix:
    //     // e.g. new materials, new real constants, nonlinear elements.
    //     update_matrix = true;
    //     solver->reuse_preconditioner(false);
    // }
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Here add the flag for ferromagnetic
    else
    {

        if (active_time_substep->get_index() > 0)
        {
            update_matrix = false;
            solver->reuse_preconditioner(true);
        }
        else
        {
            NodeBCObject *node_dof_constraints = dynamic_cast<NodeBCObject *>(&eq_sys->get_system(0).get_constraint_object());
            auto current_material_object = load_step->get_materials();
            auto old_material_object = prev_load_step->get_materials();
            auto current_real_constant_object = load_step->get_real_constants();
            auto old_real_constant_object = prev_load_step->get_real_constants();
            if (std::abs(active_time_substep->get_dt() - old_time_substep->get_dt()) >
                    relative_epsilon * std::max(std::abs(active_time_substep->get_dt()),
                                                std::abs(old_time_substep->get_dt())) ||
                node_dof_constraints->get_update_boundary_conditions() ||
                active_time_substep->get_transient_effect() != old_time_substep->get_transient_effect() /*||
                current_material_object->get_data_map() != old_material_object->get_data_map() ||
                current_real_constant_object->get_data_map() != old_real_constant_object->get_data_map())*/
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
    return update_matrix;
}
// void AssemblerLinear::compute_derived_quantities()
// {
//     const MeshBase &mesh = eq_sys.get_mesh();
//     const unsigned int dim = mesh.mesh_dimension();

//     ExplicitSystem &derived_quantities_system = eq_sys.get_system<ExplicitSystem>("DerivedQuantitiesSystem");
//     ExplicitSystem &derived_quantities_system_elem = eq_sys.get_system<ExplicitSystem>("DerivedQuantitiesSystemElem");
//     TransientLinearImplicitSystem &system = eq_sys.get_system<TransientLinearImplicitSystem>("ApecMLinearSystem");

//     const DofMap &dof_map = system.get_dof_map();
//     const DofMap &derived_quantities_system_dof_map = derived_quantities_system.get_dof_map();
//     const DofMap &derived_quantities_system_elem_dof_map = derived_quantities_system_elem.get_dof_map();

//     std::vector<unsigned int> A_num_var;
//     std::vector<unsigned int> B_num_var;
//     std::vector<unsigned int> J_num_var;
//     unsigned int v_var;

//     A_num_var.push_back(system.variable_number("A_z"));
//     B_num_var.push_back(derived_quantities_system.variable_number("B_x"));
//     B_num_var.push_back(derived_quantities_system.variable_number("B_y"));
//     J_num_var.push_back(derived_quantities_system.variable_number("J_elem_z"));

//     if (dim == 3)
//     {
//         A_num_var.push_back(system.variable_number("A_x"));
//         A_num_var.push_back(system.variable_number("A_y"));
//         v_var = system.variable_number("V");

//         B_num_var.push_back(derived_quantities_system.variable_number("B_z"));
//         J_num_var.push_back(derived_quantities_system.variable_number("J_elem_x"));
//         J_num_var.push_back(derived_quantities_system.variable_number("J_elem_y"));
//     }

//     std::shared_ptr<TimeStepper> timeStepper = eq_sys.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper");
//     Real dt = timeStepper->get_time_substep_dt();
//     bool trans_effects = timeStepper->get_time_steps()[timeStepper->get_current_time_step()].get_transient_effect();

//     // Storage for the stress dof indices on each element
//     std::vector<std::vector<dof_id_type>> A_dof_indices_var;
//     std::vector<std::vector<dof_id_type>> B_dof_indices_var;
//     std::vector<std::vector<dof_id_type>> J_dof_indices_var;
//     std::vector<dof_id_type> A_scalar_dof_indices_var;
//     std::vector<dof_id_type> B_scalar_dof_indices_var;
//     std::vector<dof_id_type> J_scalar_dof_indices_var;
//     std::vector<dof_id_type> V_indices_var;

//     std::shared_ptr<LoadsObject> loads_object = eq_sys.parameters.get<std::shared_ptr<LoadsObject>>("loads_object");

//     DenseVector<Number> Ae;
//     DenseVector<Number> Ae_dt;
//     DenseVector<Number> Ve;
//     DenseVector<Number> Ve_dt;
//     DenseVector<Number> Be;
//     // DenseVector<Number> Fe;
//     DenseVector<Number> Be_q;
//     DenseVector<Number> J_elem;
//     DenseVector<Number> J_elem_A;
//     DenseVector<Number> J_elem_V;
//     // DenseVector<Number> B_elem;
//     // DenseVector<Number> F_elem;
//     char n_J_dofs;
//     if (dim == 2)
//     {
//         n_J_dofs = 1;
//         J_elem.resize(n_J_dofs);
//         // B_elem.resize(2);
//         // F_elem.resize(2);
//     }
//     else
//     {
//         n_J_dofs = 3;
//         J_elem.resize(n_J_dofs);
//         B_elem.resize(n_J_dofs);
//         F_elem.resize(n_J_dofs);
//     }
//     char n_A_dofs;
//     char n_V_dofs;
//     char n_B_dofs;
//     bool conductive;
//     // std::set<subdomain_id_type> *conductive_subdomain = eq_sys.parameters.get<std::set<subdomain_id_type> *>("conductive_subdomain");
//     unsigned int local_elem_idx = 0;
//     for (const auto &elem : mesh.active_local_element_ptr_range())
//     {
//         loads_object->local_elem_id = local_elem_idx;

//         if (elem->subdomain_id() == 1)
//             conductive = true;
//         else
//             conductive = false;

//         if (!trans_effects)
//             conductive = false;

//         loads_object->compute_material_properties(conductive);

//         for (unsigned int var = 0; var < A_num_var.size(); var++)
//         {
//             dof_map.dof_indices(elem, A_scalar_dof_indices_var, A_num_var[var]);
//             A_dof_indices_var.push_back(A_scalar_dof_indices_var);
//         }
//         for (unsigned int var = 0; var < B_num_var.size(); var++)
//         {
//             dof_map.dof_indices(elem, A_scalar_dof_indices_var, A_num_var[var]);
//             A_dof_indices_var.push_back(A_scalar_dof_indices_var);
//         }
//         for (unsigned int var = 0; var < A_num_var.size(); var++)
//         {
//             dof_map.dof_indices(elem, A_scalar_dof_indices_var, A_num_var[var]);
//             A_dof_indices_var.push_back(A_scalar_dof_indices_var);
//         }
//         dof_map.dof_indices(elem, A_indices_var, a_var);
//         n_A_dofs = std::size(A_indices_var);
//         Ae.resize(n_A_dofs);
//         Ae_dt.resize(n_A_dofs);

//         derived_quantities_system_dof_map.dof_indices(elem, B_indices_var, b_n_var);
//         n_B_dofs = std::size(B_indices_var);
//         Be_q.resize(n_B_dofs);
//         Be.resize(n_B_dofs);

//         derived_quantities_system_elem_dof_map.dof_indices(elem, J_indices_var, j_elem_var);
//         B_elem.zero();
//         F_elem.zero();
//         if (conductive)
//         {
//             J_elem.zero();
//             if (dim == 3)
//             {
//                 dof_map.dof_indices(elem, V_indices_var, v_var);
//                 n_V_dofs = std::size(V_indices_var);
//                 Ve.resize(n_V_dofs);
//                 Ve_dt.resize(n_V_dofs);
//                 J_elem_A.zero();
//                 J_elem_V.zero();
//             }
//         }

//         for (std::size_t j = 0; j < n_A_dofs; j++)
//         {
//             Ae(j) = system.current_solution(A_indices_var[j]);
//             if (conductive)
//             {
//                 Ae_dt(j) = (Ae(j) - system.old_solution(A_indices_var[j])) / dt;

//                 if (dim == 3)
//                 {
//                     if (j < n_V_dofs)
//                     {
//                         Ve_dt(j) = (system.current_solution(V_indices_var[j]) - system.old_solution(V_indices_var[j])) / dt;
//                     }
//                 }
//             }
//         }

//         loads_object->aae_matrix_derived_map[local_elem_idx].vector_mult(Be_q, Ae);

//         // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
//         // Note that this gives a discontinuous stress plot on element boundaries, which is
//         // appropriate. We then also get the von Mises stress from the projected stress tensor.
//         loads_object->me_matrix_derived_map[local_elem_idx].cholesky_solve(Be_q, Be);

//         unsigned int j_qp_size;
//         if (conductive)
//         {
//             if (dim == 2)
//             {
//                 j_qp_size = loads_object->jae_matrix_derived_map[local_elem_idx].m();
//                 J_elem_A.resize(j_qp_size);
//                 loads_object->jae_matrix_derived_map[local_elem_idx].vector_mult(J_elem_A, Ae);
//                 for (char k = 0; k < n_J_dofs; k++)
//                     for (char qp = 0; qp < j_qp_size; qp++)
//                         J_elem(k) -= J_elem_A(qp + k) / loads_object->rs_matrix[k] / (j_qp_size / n_J_dofs);
//             }
//             else
//             {
//                 j_qp_size = loads_object->jae_matrix_derived_map[local_elem_idx].m();
//                 J_elem_A.resize(j_qp_size);
//                 J_elem_V.resize(j_qp_size);
//                 loads_object->jae_matrix_derived_map[local_elem_idx].vector_mult(J_elem_A, Ae_dt);
//                 loads_object->jve_matrix_derived_map[local_elem_idx].vector_mult(J_elem_V, Ve_dt);
//                 for (unsigned int k = 0; k < n_J_dofs; k++)
//                     for (unsigned int qp = 0; qp < j_qp_size / n_J_dofs; qp++)
//                         J_elem(k) += (-J_elem_A(qp * 3 + k) - J_elem_V(qp * 3 + k)) / loads_object->rs_matrix[k] / (j_qp_size / n_J_dofs);
//             }
//         }
//         else
//         {
//             for (unsigned int k = 0; k < n_J_dofs; k++)
//             {
//                 J_elem(k) = loads_object->compute_js(k);
//             }
//         }

//         dof_id_type dof_index;
//         // Loop over B_indices_var_values
//         for (unsigned int i = 0; i < n_B_dofs; i++)
//         {
//             dof_index = B_indices_var[i];
//             if ((derived_quantities_system.solution->first_local_index() <= dof_index) &&
//                 (dof_index < derived_quantities_system.solution->last_local_index()))
//                 derived_quantities_system.solution->set(dof_index, Be(i));
//         }

//         for (unsigned int i = 0; i < n_J_dofs; i++)
//         {
//             dof_index = J_indices_var[i];
//             if ((derived_quantities_system_elem.solution->first_local_index() <= dof_index) &&
//                 (dof_index < derived_quantities_system_elem.solution->last_local_index()))
//                 derived_quantities_system_elem.solution->set(dof_index, J_elem(i));
//         }
//         // Should call close and update when we set vector entries directly
//         local_elem_idx++;
//     }
//     derived_quantities_system.solution->close();
//     derived_quantities_system.update();
//     derived_quantities_system_elem.solution->close();
//     derived_quantities_system_elem.update();
// }
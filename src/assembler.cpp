#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
#include "../include/session.h" // Include the header file for Session
#include "libmesh/string_to_enum.h"
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
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
#include "vtkDataArray.h"
#include "../include/json.hpp"
#include "../include/assembler.h"
#include "../include/element_builder.h" // Include the header file for ElementBuilder
#include <map>                          // Include for std::map
#include <vtkCellData.h>
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
    const MeshBase &mesh = eq_sys.get_mesh();

    // Get a reference to the LinearImplicitSystem we are solving
    TransientLinearImplicitSystem &system = eq_sys.get_system<TransientLinearImplicitSystem>("PlasmaArcLinearSystem");

    // A reference to the DofMap object for this system.  The DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the DofMap
    // in future examples.
    const DofMap &dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first variable in the system.

    // this doesn't work in 2D
    const unsigned int a_var = system.variable_number("A");
    unsigned int n_A_dofs;
    std::vector<dof_id_type> A_indices_var;

    std::string fe_family = eq_sys.parameters.get<Session *>("session")->config["Finite Elements Settings"]["fe_family"].get<std::string>();
    std::string formulation = eq_sys.parameters.get<Session *>("session")->config["Finite Elements Settings"]["formulation"].get<std::string>();
    std::string order = eq_sys.parameters.get<Session *>("session")->config["Finite Elements Settings"]["approx_order"].get<std::string>();
    unsigned int dim = eq_sys.get_mesh().mesh_dimension();

    FEType fe_type_vec(Utility::string_to_enum<Order>(order), LAGRANGE_VEC);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr<FEVectorBase> fe_vec(FEVectorBase::build(eq_sys.get_mesh().mesh_dimension(), fe_type_vec));

    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(eq_sys.get_mesh().mesh_dimension(), fe_type_vec.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe_vec->attach_quadrature_rule(&qrule);

    // Get a constant reference to the Finite Element type
    // for the first variable in the system.
    FEType fe_type_scalar(Utility::string_to_enum<Order>(order), LAGRANGE);
    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr<FEBase> fe_scalar(FEBase::build(dim, fe_type_scalar));

    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule_scalar(dim, fe_type_scalar.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe_scalar->attach_quadrature_rule(&qrule_scalar);
    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    // We begin with the element Jacobian * quadrature weight at each
    // integration point.
    // const std::vector<Real> &JxW = fe->get_JxW();

    // // The element shape functions evaluated at the quadrature points.
    // const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // // The element shape function gradients for the velocity
    // // variables evaluated at the quadrature points.
    // const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // const Real time_factor = eq_sys.parameters.get<Real>("time_factor");
    // const unsigned int time_step = eq_sys.parameters.get<unsigned int>("time_step");

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe". More detail is in example 3.
    DenseMatrix<Number> KKe;
    DenseMatrix<Number> CCe;
    DenseMatrix<Number> Ke;
    DenseMatrix<Number> Ce;

    DenseVector<Number> FFe;
    DenseVector<Number> Fe;
    std::vector<Real> old_solution_values;
    auto timeStepper = eq_sys.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper");
    // check correct time stepper get
    // Get the cell data from the vtkGrid

    // Get the materials properties
    auto element_data_manager = eq_sys.parameters.get<std::shared_ptr<ElementDataManager>>("element_data_manager");
    element_data_manager->update_time_step();
    const TimeStep &time_step = timeStepper->get_time_steps()[timeStepper->get_current_time_step()];
    ElementBuilder elem_builder(Ke, Ce, Fe, element_data_manager, fe_vec, fe_scalar);

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();
    check_update_matrix();
    unsigned int n_dofs;
    dof_id_type local_elem_id = 0;
    bool conductive;
    std::set<subdomain_id_type> *conductive_subdomain = eq_sys.parameters.get<std::set<subdomain_id_type> *>("conductive_subdomain");

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.
        if (conductive_subdomain->find(elem->subdomain_id()) != conductive_subdomain->end())
            conductive = true;
        else
            conductive = false;
        element_data_manager->local_elem_id = local_elem_id;
        if (timeStepper->get_cumulative_substep() == 1) // Here add the check for nonlinear elements
        {
            // this can be done at the beginning only
            dof_map.dof_indices(elem, dof_indices);
            element_data_manager->dof_indices[local_elem_id] = dof_indices;
            // Cache the number of degrees of freedom on this element, for
            // use as a loop bound later.  We use cast_int to explicitly
            // convert from size() (which may be 64-bit) to unsigned int
            // (which may be 32-bit but which is definitely enough to count
            // *local* degrees of freedom.
            n_dofs = cast_int<unsigned int>(dof_indices.size());
        }
        else
        {
            dof_indices = element_data_manager->dof_indices[local_elem_id];
            n_dofs = cast_int<unsigned int>(dof_indices.size());
        }
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        if (timeStepper->get_cumulative_substep() == 1) // Here add the check for nonlinear elements
        {
            // this can be done at the beginning only
            // dof_map.dof_indices(elem, dof_indices);
            // this->insert_dof_indices(elem_id, dof_indices);

            // // Cache the number of degrees of freedom on this element, for
            // // use as a loop bound later.  We use cast_int to explicitly
            // // convert from size() (which may be 64-bit) to unsigned int
            // // (which may be 32-bit but which is definitely enough to count
            // // *local* degrees of freedom.
            // const unsigned int n_dofs =
            //     cast_int<unsigned int>(dof_indices.size());
            // Compute the element-specific data for the current
            // element.  This involves computing the location of the
            // quadrature points (q_point) and the shape functions
            // (phi, dphi) for the current element.
            if (dim == 3)
                fe_vec->reinit(elem);
            if (conductive || dim == 2)
                fe_scalar->reinit(elem);

            element_data_manager->compute_material_properties(conductive);
            // With one variable, we should have the same number of degrees
            // of freedom as shape functions.
            // libmesh_assert_equal_to(n_dofs, phi.size());

            Ke.resize(n_dofs, n_dofs);
            if (conductive)
                Ce.resize(n_dofs, n_dofs);

            elem_builder.build_element_matrix(elem->n_nodes(), conductive, formulation, dim, fe_family);

            DenseMatrix<Number> AAe;
            DenseMatrix<Number> Me;

            DenseMatrix<Number> JAe;
            DenseMatrix<Number> JVe;
            dof_map.dof_indices(elem, A_indices_var, a_var);
            // here maybe you will want to add the A dof indices to a map
            n_A_dofs = std::size(A_indices_var);
            if (dim == 2)
            {
                AAe.resize(2 * n_A_dofs, n_A_dofs);
                Me.resize(2 * n_A_dofs, 2 * n_A_dofs);
                if (conductive)
                    JAe.resize(fe_scalar->n_quadrature_points(), n_A_dofs);
                // JVe.resize(1, n_A_dofs);
            }
            else
            {
                AAe.resize(n_A_dofs, n_A_dofs);
                Me.resize(n_A_dofs, n_A_dofs);
                if (conductive)
                {
                    JAe.resize(3 * fe_scalar->n_quadrature_points(), n_A_dofs);
                    // fix this / 3
                    JVe.resize(3 * fe_scalar->n_quadrature_points(), n_A_dofs / 3);
                }
            }

            elem_builder.build_element_matrix_derived(AAe, Me, JAe, JVe, conductive, formulation, dim, fe_family, n_A_dofs);
            element_data_manager->aae_matrix_derived_map[local_elem_id] = AAe;
            element_data_manager->me_matrix_derived_map[local_elem_id] = Me;
            if (conductive)
            {
                element_data_manager->jae_matrix_derived_map[local_elem_id] = JAe;
                if (dim == 3)
                    element_data_manager->jve_matrix_derived_map[local_elem_id] = JVe;
            }

            element_data_manager->ke_matrix_map[local_elem_id] = Ke;
            if (conductive)
            {
                element_data_manager->ce_matrix_map[local_elem_id] = Ce;
                // if (elem_id == 969)
                // {
                //     std::ofstream ke_file("../simulations/two_wires_3D/matrix/Ke_matrix"
                //                           ".txt");
                //     std::ofstream ce_file("../simulations/two_wires_3D/matrix/Ce_matrix.txt");

                //     if (ke_file.is_open() && ce_file.is_open())
                //     {
                //         Ke.print_scientific(ke_file);
                //         Ce.print_scientific(ce_file);
                //     }
                // }
            }
            else
            {
                Fe.resize(n_dofs);
                // Here remove conductive from list of parameters
                elem_builder.build_element_rhs(elem->n_nodes(), conductive, formulation, dim, fe_family);
                element_data_manager->rhs_map[local_elem_id] = Fe;
            }
        }
        // in a new time step we could have a different dt, different bcs, different time integration effects
        // Here add the flag for ferromagnetic
        if (update_matrix)
        {
            // if (elem_id == 0)
            // {
            //     // libMesh::out << "Time step changed from " << dt2 << " to " << dt1 << std::endl;
            //     libMesh::out << "Computing new system matrix" << std::endl;
            // }
            KKe.resize(n_dofs, n_dofs);
            KKe += element_data_manager->ke_matrix_map[local_elem_id];
            if (conductive)
            {
                if (time_step.get_transient_effect())
                {
                    CCe.resize(n_dofs, n_dofs);
                    CCe += element_data_manager->ce_matrix_map[local_elem_id];
                    CCe.scale(1 / time_step.get_substeps()[timeStepper->get_current_substep()].get_dt());
                    KKe += CCe;
                    // if (elem_id == 969)
                    // {
                    //     std::ofstream kke_file("../simulations/two_wires_3D/matrix/Kke_matrix.txt");
                    //     if (kke_file.is_open())
                    //     {
                    //         KKe.print_scientific(kke_file);
                    //     }
                    // }
                }
            }
            dof_map.constrain_element_matrix(KKe, dof_indices);
            matrix.add_matrix(KKe, dof_indices);
        }
        FFe.resize(n_dofs);
        if (!conductive)
        {
            FFe += element_data_manager->rhs_map[local_elem_id];
            if (dim == 2)
            {
                FFe.scale(element_data_manager->compute_js(0));
            }
            // here understand if it is better to copy the object or to sum from 0
            else
            {
                Real js[3];
                for (unsigned int var = 0; var < 3; var++)
                {
                    js[var] = element_data_manager->compute_js(var);
                    for (unsigned int y = 0; y < elem->n_nodes(); y++)
                    {
                        FFe(y * 3 + var) *= js[var];
                    }
                }
            }
            // dof_map.constrain_element_vector(FFe, dof_indices);
            // system.rhs->add_vector(FFe, dof_indices);
        }
        // here I constrain but I modify later, wrong
        // }
        if (conductive && time_step.get_transient_effect())
        {
            system.old_local_solution->get(dof_indices, old_solution_values);

            element_data_manager->ce_matrix_map[local_elem_id].vector_mult(FFe, old_solution_values);
            FFe.scale(1 / time_step.get_substeps()[timeStepper->get_current_substep()].get_dt());
            // if (elem_id == 969)
            // {
            //     std::ofstream ffe_file("../simulations/two_wires_3D/matrix/FFe_matrix" + std::to_string(time_step.get_step()) + ".txt");
            //     if (ffe_file.is_open())
            //     {
            //         FFe.print_scientific(ffe_file);
            //     }
            // }
        }
        dof_map.constrain_element_vector(FFe, dof_indices);
        system.rhs->add_vector(FFe, dof_indices);
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
        local_elem_id++;
    }
    // That's it.  We don't need to do anything else to the
    // PerfLog.  When it goes out of scope (at this function return)
    // it will print its log to the screen. Pretty easy, huh?
}

void AssemblerLinear::check_update_matrix()
{
    update_matrix = false;

    const Real relative_epsilon = 1e-9;
    Real dt1;
    Real dt2;
    PetscLinearSolver<Number> *solver = dynamic_cast<PetscLinearSolver<Number> *>(eq_sys.get_system<TransientLinearImplicitSystem>("PlasmaArcLinearSystem").get_linear_solver());
    auto timeStepper = eq_sys.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper");

    dt1 = timeStepper->get_time_substep_dt();
    dt2 = timeStepper->get_prev_time_substep_dt();
    // in a new time step we could have a different dt, different bcs, different time integration effects
    // Here add the flag for ferromagnetic
    if (timeStepper->get_cumulative_substep() == 1)
    {
        update_matrix = true;
        solver->reuse_preconditioner(false);
    }
    else
    {
        if (timeStepper->get_current_substep() > 0)
        {
            update_matrix = false;
            solver->reuse_preconditioner(true);
        }
        else
        {
            if (std::abs(dt1 - dt2) > relative_epsilon * std::max(std::abs(dt1), std::abs(dt2)) ||
                updated_bcs ||
                timeStepper->get_time_steps()[timeStepper->get_current_time_step()].get_transient_effect() != timeStepper->get_time_steps()[timeStepper->get_current_time_step() - 1].get_transient_effect())
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
    return;
}

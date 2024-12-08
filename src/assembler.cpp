#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
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
#include <map> // Include for std::map
#include <vtkCellData.h>

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
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(0);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

    // A 2nd order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    // We begin with the element Jacobian * quadrature weight at each
    // integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    const Real time_factor = eq_sys.parameters.get<Real>("time_factor");
    const unsigned int time_step = eq_sys.parameters.get<unsigned int>("time_step");

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe". More detail is in example 3.
    DenseMatrix<Number> Ke;
    DenseMatrix<Number> Ce;
    DenseVector<Number> Fe;

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides
    DenseSubVector<Number> F[4] = {DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe), DenseSubVector<Number>(Fe)};

    // References to velocity-velocity coupling blocks
    DenseSubMatrix<Number> KAA[3][3] = {{DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
                                        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
                                        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}};

    // References to velocity-pressure coupling blocks
    DenseSubMatrix<Number> KAV[3] = {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)};
    DenseSubMatrix<Number> KVA[3] = {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)};
    DenseSubMatrix<Number> KVV(Ke);

    // References to velocity-velocity coupling blocks
    DenseSubMatrix<Number> CAA[3][3] = {{DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)},
                                        {DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)},
                                        {DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)}};

    // References to velocity-pressure coupling blocks
    DenseSubMatrix<Number> CAV[3] = {DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)};
    DenseSubMatrix<Number> CVA[3] = {DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce), DenseSubMatrix<Number>(Ce)};
    DenseSubMatrix<Number> CVV(Ce);
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();

    std::vector<Real> old_solution_values;

    // You will have to fix this for transient problems
    std::string array_name = "Js_" + std::to_string(time_step);
    vtkDataArray *jsArray = cell_data->GetArray(array_name.c_str());

    array_name = "ansys_material_type";
    vtkDataArray *matArray = cell_data->GetArray(array_name.c_str());
    if (!jsArray)
    {
        // Get the value from the "Js" array using the element ID
        libmesh_error_msg("Error: Js_" + std::to_string(time_step) + " array not found in vtkCellData.");
    }
    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.  See example 3 for a discussion of the
    // element iterators.
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.
        perf_log.push("elem init");

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        // Cache the number of degrees of freedom on this element, for
        // use as a loop bound later.  We use cast_int to explicitly
        // convert from size() (which may be 64-bit) to unsigned int
        // (which may be 32-bit but which is definitely enough to count
        // *local* degrees of freedom.
        const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        const unsigned int subdomain_id = elem->subdomain_id();

        dof_id_type el_id = elem->id();

        // Retrieve the material number for the current element
        const unsigned int material_number = static_cast<unsigned int>(matArray->GetTuple1(el_id)); // Retrieve the material number using GetTuple1

        // Access the corresponding material properties from materials
        const MaterialsProperties &mat_properties = materials.at(material_number);

        // Create a diagonal matrix with murx, mury, and murz as the diagonal terms
        std::vector<std::vector<double>> mu_matrix(3, std::vector<double>(3, 0.0));
        mu_matrix[0][0] = mat_properties.murx;
        mu_matrix[1][1] = mat_properties.mury;
        mu_matrix[2][2] = mat_properties.murz;

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to(n_dofs, phi.size());

        Ke.resize(n_dofs, n_dofs);

        Ce.resize(n_dofs, n_dofs);

        Fe.resize(n_dofs);

        old_solution_values.resize(n_dofs);

        perf_log.pop("elem init");
        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        if (subdomain_id == 1 && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
        {
            // Create a diagonal matrix with rsvx, rsvy, and rsvz as the diagonal terms
            std::vector<std::vector<double>> rs_matrix(3, std::vector<double>(3, 0.0));
            rs_matrix[0][0] = mat_properties.rsvx;
            rs_matrix[1][1] = mat_properties.rsvy;
            rs_matrix[2][2] = mat_properties.rsvz;

            libmesh_error_msg("Error: Conducting elements with 3D nodal elements are not implemented.");
            build_solid97_cond();
            // Ke.resize(n_dofs, n_dofs);

            // Fe.resize(n_dofs);

            // // Stop logging the shape function initialization.
            // // If you forget to stop logging an event the PerfLog
            // // object will probably catch the error and abort.
            // perf_log.pop("elem init");

            // // Now we will build the element matrix.  This involves
            // // a double loop to integrate the test functions (i) against
            // // the trial functions (j).
            // //
            // // We have split the numeric integration into two loops
            // // so that we can log the matrix and right-hand-side
            // // computation separately.
            // //
            // // Now start logging the element matrix computation
            // perf_log.push("Ke");

            // // Here I will have to differentiate based on the formulation, transient, dimension, fe_family, etc.
            // for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
            // {
            //     // Values to hold the old solution & its gradient.
            //     Number u_old = 0.;
            //     Gradient grad_u_old;

            //     // Compute the old solution & its gradient.
            //     for (std::size_t l = 0; l < phi.size(); l++)
            //     {
            //         u_old += phi[l][qp] * system.old_solution(dof_indices[l]);
            //     }
            //     for (unsigned int i = 0; i != n_dofs; i++)
            //         for (unsigned int j = 0; j != n_dofs; j++)
            //             // The matrix contribution
            //             Ke(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp] +

            //                                    .5 * dt * (
            //                                                  // Convection term
            //                                                  (velocity * dphi[j][qp]) * phi[i][qp] +

            //                                                  // Diffusion term
            //                                                  0.01 * (dphi[i][qp] * dphi[j][qp])));
            // }
            // // Stop logging the matrix computation
            // perf_log.pop("Ke");

            // // Now we build the element right-hand-side contribution.
            // // This involves a single loop in which we integrate the
            // // "forcing function" in the PDE against the test functions.
            // //
            // // Start logging the right-hand-side computation
            // perf_log.push("Fe");
            // // Here I will have to differentiate based on the formulation, transient, dimension, fe_family, etc.
            // dof_id_type el_id = elem->id();
            // if (dim == 3)
            // {
            //     Real tuple[3];
            //     jsArray->GetTuple(el_id, tuple);
            //     RealGradient js(4 * libMesh::pi * 1e-7 * tuple[0],
            //                     4 * libMesh::pi * 1e-7 * tuple[1],
            //                     4 * libMesh::pi * 1e-7 * tuple[2]);
            //     for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
            //     {
            //         // Add the RHS contribution
            //         for (unsigned int i = 0; i != n_dofs; i++)
            //             Fe(i) += JxW[qp] * fxy * phi[i][qp];
            //     }
            // }

            // // Stop logging the right-hand-side computation
            // perf_log.pop("Fe");
        }
        if (subdomain_id == 0 && formulation == "AV" && dim == 3 && fe_family == "LAGRANGE")
        {
            build_solid97();
        }
        if (subdomain_id == 1 && formulation == "AV" && fe_family == "NEDELEC_ONE")
        {
            // Create a diagonal matrix with rsvx, rsvy, and rsvz as the diagonal terms
            std::vector<std::vector<double>> rs_matrix(3, std::vector<double>(3, 0.0));
            rs_matrix[0][0] = mat_properties.rsvx;
            rs_matrix[1][1] = mat_properties.rsvy;
            rs_matrix[2][2] = mat_properties.rsvz;
            build_solid236_cond();
        }
        if (subdomain_id == 0 && formulation == "AV" && fe_family == "NEDELEC_ONE")
        {
            build_solid236();
        }
        if (subdomain_id == 1 && formulation == "AV" && dim == 2)
        {
            // Create a diagonal matrix with rsvx, rsvy, and rsvz as the diagonal terms
            std::vector<std::vector<double>> rs_matrix(3, std::vector<double>(3, 0.0));
            rs_matrix[0][0] = mat_properties.rsvx;
            rs_matrix[1][1] = mat_properties.rsvy;
            rs_matrix[2][2] = mat_properties.rsvz;
            build_plane13_cond();
        }
        if (subdomain_id == 0 && formulation == "AV" && dim == 2)
        {

            build_plane13(Ke, Fe, elem, dof_indices, dof_map, el_id, qrule, phi, dphi, JxW, perf_log, jsArray, mu_matrix[0][0], n_dofs);
        }

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        // Also, note that here we call heterogenously_constrain_element_matrix_and_vector
        // to impose a inhomogeneous Dirichlet boundary conditions.
        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The SparseMatrix::add_matrix()
        // and NumericVector::add_vector() members do this for us.
        // Start logging the insertion of the local (element)
        // matrix and vector into the global matrix and vector
        // LOG_SCOPE_WITH("matrix insertion", "", perf_log);

        matrix.add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }

    // That's it.  We don't need to do anything else to the
    // PerfLog.  When it goes out of scope (at this function return)
    // it will print its log to the screen. Pretty easy, huh?
}

// Temporary placeholder for build_solid97_cond
void AssemblerLinear::build_solid97_cond(DenseMatrix<Number> &Ke,
                                         DenseVector<Number> &Fe,
                                         DenseMatrix<Number> &Ce,
                                         DenseSubMatrix<Number> **KAA,
                                         DenseSubMatrix<Number> *KAV,
                                         DenseSubMatrix<Number> *KVA,
                                         DenseSubMatrix<Number> KVV,
                                         DenseSubMatrix<Number> **CAA,
                                         DenseSubMatrix<Number> *CAV,
                                         DenseSubMatrix<Number> *CVA,
                                         DenseSubMatrix<Number> CVV,
                                         DenseSubVector<Number> *F,
                                         const Elem *elem,
                                         const std::vector<dof_id_type> &dof_indices,
                                         const DofMap &dof_map,
                                         dof_id_type el_id,
                                         unsigned int n_dofs,
                                         Real dt,
                                         bool time_factor,
                                         TransientLinearImplicitSystem &system,
                                         QGauss &qrule,
                                         const std::vector<std::vector<Real>> &phi,
                                         const std::vector<std::vector<RealGradient>> &dphi,
                                         const std::vector<Real> &JxW,
                                         PerfLog &perf_log,
                                         std::vector<std::vector<double>> &sigma_matrix,
                                         std::vector<std::vector<double>> &mu_matrix,
                                         std::vector<Real> old_solution_values)
{
    // Numeric ids corresponding to each variable in the system
    // const unsigned int ax_var = system.variable_number("AX");
    // const unsigned int ay_var = system.variable_number("AY");
    // const unsigned int az_var = system.variable_number("AZ");
    const unsigned int num_vars = 4;
    std::vector<std::vector<dof_id_type>> dof_indices_var(num_vars);

    // std::vector<dof_id_type> dof_indices_ax;
    // std::vector<dof_id_type> dof_indices_ay;
    // std::vector<dof_id_type> dof_indices_az;
    // dof_map.dof_indices(elem, dof_indices_ax, ax_var);
    // dof_map.dof_indices(elem, dof_indices_ay, ay_var);
    // dof_map.dof_indices(elem, dof_indices_az, az_var);
    // const unsigned int n_ax_dofs = dof_indices_ax.size();
    // const unsigned int n_ay_dofs = dof_indices_ay.size();
    // const unsigned int n_az_dofs = dof_indices_az.size();

    for (unsigned int var = 0; var < num_vars; var++)
        dof_map.dof_indices(elem, dof_indices_var[var], var);

    const unsigned int n_var_dofs = dof_indices_var[0].size();

    for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
        for (unsigned int var_j = 0; var_j < num_vars - 1; var_j++)
        {
            KAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);
            CAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);
        }

    for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
    {
        KAV[var_i].reposition(var_i * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
        KVA[var_i].reposition((num_vars - 1) * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);
        CAV[var_i].reposition(var_i * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
        CVA[var_i].reposition((num_vars - 1) * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);
    }
    for (unsigned int var = 0; var < num_vars; var++)
        F[var].reposition(var * n_var_dofs, n_var_dofs);

    KVV.reposition((num_vars - 1) * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
    CVV.reposition((num_vars - 1) * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    perf_log.push("Ke");
    Real CAV_value;
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        // Real js[3];
        // jsArray->GetTuple(el_id, js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                {
                    // Here you will put the time factor as well
                    CAV_value = JxW[qp] * sigma_matrix[k][k] * phi[dof_i][qp] * dphi[dof_j][qp](k); // Not sure
                    CAV[k](dof_i, dof_j) += CAV_value;
                    CVA[k](dof_j, dof_i) += CAV_value;
                    CVV(dof_i, dof_j) += JxW[qp] * dphi[dof_i][qp](k) * sigma_matrix[k][k] * dphi[dof_j][qp](k);
                    KAA[k][k](dof_i, dof_j) += JxW[qp] / (4 * libMesh::pi * 1e-7 * mu_matrix[k][k]) *
                                               (dphi[dof_i][qp] * dphi[dof_j][qp]);
                    CAA[k][k](dof_i, dof_j) += JxW[qp] * sigma_matrix[k][k] *
                                               (phi[dof_i][qp] * phi[dof_j][qp]);
                }

        Ce.scale(1 / dt);
        Ke += Ce;
        // NumberVectorValue U_old;
        // for (unsigned int k = 0; k < num_vars; k++)
        // {
        //     for (unsigned int dof_i = 0; dof_i < n_var_dofs; k++)
        //     {
        //         U_old(k) += phi[dof_i][qp] * system.old_solution(dof_indices_var[k][dof_i]);
        //     }
        // }
        // Now we will build the element matrix.  This involves
        // a double loop to integrate the test functions (i) against
        // the trial functions (j).
        //
        // We have split the numeric integration into two loops
        // so that we can log the matrix and right-hand-side
        // computation separately.
        //

        // Stop logging the matrix computation
        perf_log.pop("Ke");

        // Now we build the element right-hand-side contribution.
        // This involves a single loop in which we integrate the
        // "forcing function" in the PDE against the test functions.
        //
        // Start logging the right-hand-side computation
        perf_log.push("Fe");
        // // Retrieve the tuple values
        // for (unsigned int k = 0; k < num_vars; k++)
        //     for (unsigned int dof_i = 0; dof_i < n_var_dofs; k++)
        //     {
        //         F[k](dof_i) += JxW[qp] * (js[k] * +U_old(k) * * / dt);
        //     }
        system.old_local_solution->get(dof_indices, old_solution_values);

        Ce.vector_mult(Fe, old_solution_values);

        // // Stop logging the right-hand-side computation
        perf_log.pop("Fe");
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}

// Temporary placeholder for build_solid97
void AssemblerLinear::build_solid97(DenseMatrix<Number> &Ke,
                                    DenseVector<Number> &Fe,
                                    DenseMatrix<Number> &Ce,
                                    DenseSubMatrix<Number> **KAA,
                                    DenseSubVector<Number> *F,
                                    const Elem *elem,
                                    const std::vector<dof_id_type> &dof_indices,
                                    const DofMap &dof_map,
                                    dof_id_type el_id,
                                    unsigned int n_dofs,
                                    Real dt,
                                    bool time_factor,
                                    TransientLinearImplicitSystem &system,
                                    QGauss &qrule,
                                    const std::vector<std::vector<Real>> &phi,
                                    const std::vector<std::vector<RealGradient>> &dphi,
                                    const std::vector<Real> &JxW,
                                    PerfLog &perf_log,
                                    std::vector<std::vector<double>> &sigma_matrix,
                                    std::vector<std::vector<double>> &mu_matrix,
                                    vtkDataArray *jsArray,
                                    std::vector<Real> old_solution_values)
{
    // Numeric ids corresponding to each variable in the system
    // const unsigned int ax_var = system.variable_number("AX");
    // const unsigned int ay_var = system.variable_number("AY");
    // const unsigned int az_var = system.variable_number("AZ");
    const unsigned int num_vars = 4;
    std::vector<std::vector<dof_id_type>> dof_indices_var(num_vars);

    // std::vector<dof_id_type> dof_indices_ax;
    // std::vector<dof_id_type> dof_indices_ay;
    // std::vector<dof_id_type> dof_indices_az;
    // dof_map.dof_indices(elem, dof_indices_ax, ax_var);
    // dof_map.dof_indices(elem, dof_indices_ay, ay_var);
    // dof_map.dof_indices(elem, dof_indices_az, az_var);
    // const unsigned int n_ax_dofs = dof_indices_ax.size();
    // const unsigned int n_ay_dofs = dof_indices_ay.size();
    // const unsigned int n_az_dofs = dof_indices_az.size();

    for (unsigned int var = 0; var < num_vars; var++)
        dof_map.dof_indices(elem, dof_indices_var[var], var);

    const unsigned int n_var_dofs = dof_indices_var[0].size();

    for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
        for (unsigned int var_j = 0; var_j < num_vars - 1; var_j++)
            KAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);

    for (unsigned int var = 0; var < num_vars; var++)
        F[var].reposition(var * n_var_dofs, n_var_dofs);

    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    perf_log.push("Ke");
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        Real js[3];
        jsArray->GetTuple(el_id, js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                {
                    KAA[k][k](dof_i, dof_j) += JxW[qp] / (4 * libMesh::pi * 1e-7 * mu_matrix[k][k]) *
                                               (dphi[dof_i][qp] * dphi[dof_j][qp]);
                    F[k](dof_i) += JxW[qp] * js[k] * phi[dof_i][qp];
                }

        // Now we will build the element matrix.  This involves
        // a double loop to integrate the test functions (i) against
        // the trial functions (j).
        //
        // We have split the numeric integration into two loops
        // so that we can log the matrix and right-hand-side
        // computation separately.
        //

        // Stop logging the matrix computation
        perf_log.pop("Ke");

        // Now we build the element right-hand-side contribution.
        // This involves a single loop in which we integrate the
        // "forcing function" in the PDE against the test functions.
        //
        // Start logging the right-hand-side computation
        perf_log.push("Fe");
        // // Retrieve the tuple values
        // for (unsigned int k = 0; k < num_vars; k++)
        //     for (unsigned int dof_i = 0; dof_i < n_var_dofs; k++)
        //     {
        //         F[k](dof_i) += JxW[qp] * (js[k] * +U_old(k) * * / dt);
        //     }

        // // Stop logging the right-hand-side computation
        perf_log.pop("Fe");
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}

// Temporary placeholder for build_solid236
void AssemblerLinear::build_solid236()
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_solid236_cond
void AssemblerLinear::build_solid236_cond()
{
    libmesh_error_msg("build_solid236_cond is not implemented yet.");
}

// Temporary placeholder for build_plane13
void AssemblerLinear::build_plane13(DenseMatrix<Number> &Ke,
                                    DenseVector<Number> &Fe,
                                    const Elem *elem,
                                    const std::vector<dof_id_type> &dof_indices,
                                    const DofMap &dof_map,
                                    dof_id_type el_id,
                                    QGauss &qrule,
                                    const std::vector<std::vector<Real>> &phi,
                                    const std::vector<std::vector<RealGradient>> &dphi,
                                    const std::vector<Real> &JxW,
                                    PerfLog &perf_log,
                                    vtkDataArray *jsArray,
                                    Real mu,
                                    unsigned int n_dofs)
{
    // Stop logging the shape function initialization.
    // If you forget to stop logging an event the PerfLog
    // object will probably catch the error and abort.
    perf_log.pop("elem init");

    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    perf_log.push("Ke");

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        for (unsigned int i = 0; i != n_dofs; i++)
            for (unsigned int j = 0; j != n_dofs; j++)
                Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]) / (4 * libMesh::pi * 1e-7 * mu);

    // Stop logging the matrix computation
    perf_log.pop("Ke");

    // Now we build the element right-hand-side contribution.
    // This involves a single loop in which we integrate the
    // "forcing function" in the PDE against the test functions.
    //
    // Start logging the right-hand-side computation
    perf_log.push("Fe");

    Real js = jsArray->GetTuple1(el_id);

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        // Add the RHS contribution
        // Check this, mu_matrix should be put in the Ke, and requires some attention, this is temporary
        for (unsigned int i = 0; i != n_dofs; i++)
            Fe(i) += JxW[qp] * js * phi[i][qp];
    }

    // Stop logging the right-hand-side computation
    perf_log.pop("Fe");
}

// Temporary placeholder for build_plane13_cond
void AssemblerLinear::build_plane13_cond(DenseMatrix<Number> &Ke,
                                         DenseVector<Number> &Fe,
                                         DenseMatrix<Number> &Ce,
                                         DenseSubMatrix<Number> **KAA,
                                         DenseSubMatrix<Number> *KAV,
                                         DenseSubMatrix<Number> *KVA,
                                         DenseSubMatrix<Number> KVV,
                                         DenseSubMatrix<Number> **CAA,
                                         DenseSubMatrix<Number> *CAV,
                                         DenseSubMatrix<Number> *CVA,
                                         DenseSubMatrix<Number> CVV,
                                         DenseSubVector<Number> *F,
                                         const Elem *elem,
                                         const std::vector<dof_id_type> &dof_indices,
                                         const DofMap &dof_map,
                                         dof_id_type el_id,
                                         unsigned int n_dofs,
                                         Real dt,
                                         bool time_factor,
                                         TransientLinearImplicitSystem &system,
                                         QGauss &qrule,
                                         const std::vector<std::vector<Real>> &phi,
                                         const std::vector<std::vector<RealGradient>> &dphi,
                                         const std::vector<Real> &JxW,
                                         PerfLog &perf_log,
                                         Real sigma,
                                         Real mu,
                                         std::vector<Real> old_solution_values)
{
    // Numeric ids corresponding to each variable in the system
    const unsigned int num_vars = 2;
    std::vector<std::vector<dof_id_type>> dof_indices_var(num_vars);

    for (unsigned int var = 0; var < num_vars; var++)
        dof_map.dof_indices(elem, dof_indices_var[var], var);

    const unsigned int n_var_dofs = dof_indices_var[0].size();

    for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
        for (unsigned int var_j = 0; var_j < num_vars - 1; var_j++)
        {
            KAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);
            CAA[var_i][var_j].reposition(var_i * n_var_dofs, var_j * n_var_dofs, n_var_dofs, n_var_dofs);
        }

    for (unsigned int var_i = 0; var_i < num_vars - 1; var_i++)
    {
        KAV[var_i].reposition(var_i * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
        KVA[var_i].reposition((num_vars - 1) * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);
        CAV[var_i].reposition(var_i * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
        CVA[var_i].reposition((num_vars - 1) * n_var_dofs, var_i * n_var_dofs, n_var_dofs, n_var_dofs);
    }
    for (unsigned int var = 0; var < num_vars; var++)
        F[var].reposition(var * n_var_dofs, n_var_dofs);

    KVV.reposition((num_vars - 1) * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
    CVV.reposition((num_vars - 1) * n_var_dofs, (num_vars - 1) * n_var_dofs, n_var_dofs, n_var_dofs);
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    perf_log.push("Ke");
    Real CAV_value;
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        // Real js[3];
        // jsArray->GetTuple(el_id, js);
        // // Here check that js is 0

        for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
            {
                // Here you will put the time factor as well
                CAV_value = JxW[qp] * sigma * phi[dof_i][qp] * dphi[dof_j][qp]; // Not sure
                CAV[0](dof_i, dof_j) += CAV_value;
                CVA[0](dof_j, dof_i) += CAV_value;
                CVV(dof_i, dof_j) += JxW[qp] * sigma * dphi[dof_i][qp] * dphi[dof_j][qp];
                KAA[0][0](dof_i, dof_j) += JxW[qp] / (4 * libMesh::pi * 1e-7 * mu) *
                                           (dphi[dof_i][qp] * dphi[dof_j][qp]);
                CAA[0][0](dof_i, dof_j) += JxW[qp] * sigma *
                                           (phi[dof_i][qp] * phi[dof_j][qp]);
            }

        Ce.scale(1 / dt);
        Ke += Ce;
        // Now we will build the element matrix.  This involves
        // a double loop to integrate the test functions (i) against
        // the trial functions (j).
        //
        // We have split the numeric integration into two loops
        // so that we can log the matrix and right-hand-side
        // computation separately.
        //

        // Stop logging the matrix computation
        perf_log.pop("Ke");

        // Now we build the element right-hand-side contribution.
        // This involves a single loop in which we integrate the
        // "forcing function" in the PDE against the test functions.
        //
        // Start logging the right-hand-side computation
        perf_log.push("Fe");
        // // Retrieve the tuple values
        // for (unsigned int k = 0; k < num_vars; k++)
        //     for (unsigned int dof_i = 0; dof_i < n_var_dofs; k++)
        //     {
        //         F[k](dof_i) += JxW[qp] * (js[k] * +U_old(k) * * / dt);
        //     }
        system.old_local_solution->get(dof_indices, old_solution_values);

        Ce.vector_mult(Fe, old_solution_values);

        // // Stop logging the right-hand-side computation
        perf_log.pop("Fe");
    }
    // libmesh_error_msg("build_solid97 is not implemented yet.");
}

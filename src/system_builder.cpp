#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
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
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/gmsh_io.h"
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

// To impose Dirichlet boundary conditions
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/zero_function.h"
#include "../include/global_resources.h"
#include "vtkDataArray.h"

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
    LinearImplicitSystem &system = eq_sys.get_system<LinearImplicitSystem>("LinearTransient");

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

    Real time_factor;
    bool timint = true;
    if (timint)
    {
        time_factor = 1.0;
    }
    else
    {
        time_factor = 0.0;
    }
    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe". More detail is in example 3.
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // Do the preparation of the matrixes and vectors based on number of variables.
    DenseSubMatrix<Number>
        K00(Ke), K01(Ke), K02(Ke), K03(Ke),
        K10(Ke), K11(Ke), K12(Ke), K13(Ke),
        K20(Ke), K21(Ke), K22(Ke), K23(Ke),
        K30(Ke), K31(Ke), K32(Ke), K33(Ke);

    DenseSubVector<Number>
        F0(Fe),
        F1(Fe),
        F2(Fe),
        F3(Fe);

    // References to momentum equation right hand sides
    std::reference_wrapper<DenseSubVector<Number>> F[4] = {F0, F1, F2, F3};

    // References to velocity-velocity coupling blocks
    std::reference_wrapper<DenseSubMatrix<Number>> K[3][3] = {{K00, K01, K02},
                                                              {K10, K11, K12},
                                                              {K20, K21, K22}};

    // References to velocity-pressure coupling blocks
    std::reference_wrapper<DenseSubMatrix<Number>> B[3] = {K03, K13, K23};
    std::reference_wrapper<DenseSubMatrix<Number>> BT[3] = {K30, K31, K32};

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe". More detail is in example 3.
    DenseMatrix<Number> Ke_no_cond;
    DenseVector<Number> Fe_no_cond;

    // Do the preparation of the matrixes and vectors based on number of variables.
    DenseSubMatrix<Number>
        K00_no_cond(Ke_no_cond), K01_no_cond(Ke_no_cond), K02_no_cond(Ke_no_cond),
        K10_no_cond(Ke_no_cond), K11_no_cond(Ke_no_cond), K12_no_cond(Ke_no_cond),
        K20_no_cond(Ke_no_cond), K21_no_cond(Ke_no_cond), K22_no_cond(Ke_no_cond);

    DenseSubVector<Number>
        F0_no_cond(Fe_no_cond),
        F1_no_cond(Fe_no_cond),
        F2_no_cond(Fe_no_cond);
    // References to momentum equation right hand sides
    std::reference_wrapper<DenseSubVector<Number>> F_no_cond[3] = {F0_no_cond, F1_no_cond, F2_no_cond};

    // References to velocity-velocity coupling blocks
    std::reference_wrapper<DenseSubMatrix<Number>> K_no_cond[3][3] = {{K00_no_cond, K01_no_cond, K02_no_cond},
                                                                      {K10_no_cond, K11_no_cond, K12_no_cond},
                                                                      {K20_no_cond, K21_no_cond, K22_no_cond}};

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_0;
    std::vector<dof_id_type> dof_indices_1;
    std::vector<dof_id_type> dof_indices_2;
    std::vector<dof_id_type> dof_indices_3;
    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();

    // You will have to fix this for transient problems
    vtkDataArray *jsArray = cell_data->GetArray("Js");
    if (!jsArray)
    {
        // Get the value from the "Js" array using the element ID
        libmesh_error_msg("Error: 'Js' array not found in vtkCellData.");
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

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to(n_dofs, phi.size());

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        if (elem->subdomain_id() == 1)
        {
            Ke.resize(n_dofs, n_dofs);

            Fe.resize(n_dofs);

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

            // Here I will have to differentiate based on the formulation, transient, dimension, fe_family, etc.
            for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
            {
                // Values to hold the old solution & its gradient.
                Number u_old = 0.;
                Gradient grad_u_old;

                // Compute the old solution & its gradient.
                for (std::size_t l = 0; l < phi.size(); l++)
                {
                    u_old += phi[l][qp] * system.old_solution(dof_indices[l]);

                    // This will work,
                    // grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
                    // but we can do it without creating a temporary like this:
                    grad_u_old.add_scaled(dphi[l][qp], system.old_solution(dof_indices[l]));
                }
                for (unsigned int i = 0; i != n_dofs; i++)
                    for (unsigned int j = 0; j != n_dofs; j++)
                        // The matrix contribution
                        Ke(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp] +

                                               .5 * dt * (
                                                             // Convection term
                                                             (velocity * dphi[j][qp]) * phi[i][qp] +

                                                             // Diffusion term
                                                             0.01 * (dphi[i][qp] * dphi[j][qp])));
            }
            // Stop logging the matrix computation
            perf_log.pop("Ke");

            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            //
            // Start logging the right-hand-side computation
            perf_log.push("Fe");
            // Here I will have to differentiate based on the formulation, transient, dimension, fe_family, etc.
            dof_id_type el_id = elem->id();
            if (dim == 3)
            {
                Real tuple[3];
                jsArray->GetTuple(el_id, tuple);
                RealGradient js(4 * libMesh::pi * 1e-7 * tuple[0],
                                4 * libMesh::pi * 1e-7 * tuple[1],
                                4 * libMesh::pi * 1e-7 * tuple[2]);
                for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                {
                    // Add the RHS contribution
                    for (unsigned int i = 0; i != n_dofs; i++)
                        Fe(i) += JxW[qp] * fxy * phi[i][qp];
                }
            }

            // Stop logging the right-hand-side computation
            perf_log.pop("Fe");

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
            LOG_SCOPE_WITH("matrix insertion", "", perf_log);

            matrix.add_matrix(Ke, dof_indices);
            system.rhs->add_vector(Fe, dof_indices);
        }
    }

    // That's it.  We don't need to do anything else to the
    // PerfLog.  When it goes out of scope (at this function return)
    // it will print its log to the screen. Pretty easy, huh?
}

void SystemBuilder::createSystem()
{
    transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
    ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
    fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
    formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
    order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
    dim = equation_systems.get_mesh().mesh_dimension();

    if (ferromagnetic)
    {
        system_name = "NonlinearTransient";
        // Add TransientSystem<NonlinearImplicitSystem> as System and return its pointer
        equation_systems.add_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
    }
    else
    {
        system_name = "LinearTransient";
        // Add TransientSystem<LinearImplicitSystem> as System and return its pointer
        equation_systems.add_system<TransientSystem<LinearImplicitSystem>>(system_name);
    }
};

void SystemBuilder::setVariables()
{
    System &system = equation_systems.get_system(system_name);

    if (formulation == "AV")
    {
        system.add_variable("AZ", Utility::string_to_enum<Order>(order),
                            Utility::string_to_enum<FEFamily>(fe_family));
        if (dim == 3 && fe_family == "LAGRANGE")
        {
            system.add_variable("AX", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family));
            system.add_variable("AY", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family));
        }
        if (transient)
        {
            std::set<subdomain_id_type> conductive_subdomain;

            conductive_subdomain.clear();
            conductive_subdomain.insert(1);
            system.add_variable("V",
                                Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family),
                                &conductive_subdomain);
        }
        AssemblerLinear *assembler = new AssemblerLinear(equation_systems, cell_data, config);
        system.attach_assemble_object(*assembler);
    }
    else
    {
        libmesh_error_msg("not implemented");
    }
    // Add more conditions as needed for other combinations
}
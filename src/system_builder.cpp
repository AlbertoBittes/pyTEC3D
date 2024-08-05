#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
#include "libmesh/string_to_enum.h"
#include "../include/loads_bc_mat_object.h"
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
void AssemblerStatic2D::assemble()
{
    // Declare a performance log.  Give it a descriptive
    // string to identify what part of the code we are
    // logging, since there may be many PerfLogs in an
    // application.
    PerfLog perf_log("Matrix Assembly");
    // Get a constant reference to the mesh object.
    const MeshBase &mesh = eq_sys.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem &system = eq_sys.get_system<LinearImplicitSystem>("Linear");

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
    QGauss qrule(dim, SECOND);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // Boundary integration requires one quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim - 1, SECOND);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    // We begin with the element Jacobian * quadrature weight at each
    // integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe". More detail is in example 3.
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();
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

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
            for (unsigned int i = 0; i != n_dofs; i++)
                for (unsigned int j = 0; j != n_dofs; j++)
                    Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);

        // Stop logging the matrix computation
        perf_log.pop("Ke");

        // Now we build the element right-hand-side contribution.
        // This involves a single loop in which we integrate the
        // "forcing function" in the PDE against the test functions.
        //
        // Start logging the right-hand-side computation
        perf_log.push("Fe");
        dof_id_type el_id = elem->id();
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            Real fxy;
            Point centroid;
            centroid = elem->vertex_average();

            fxy = 4 * libMesh::pi * 1e-7 * jsArray->GetTuple1(el_id);
            // Add the RHS contribution
            for (unsigned int i = 0; i != n_dofs; i++)
                Fe(i) += JxW[qp] * fxy * phi[i][qp];
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

    // That's it.  We don't need to do anything else to the
    // PerfLog.  When it goes out of scope (at this function return)
    // it will print its log to the screen. Pretty easy, huh?
}

// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void AssemblerStatic3D::assemble()
{
    // Get a constant reference to the mesh object.
    const MeshBase &mesh = eq_sys.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem &system = eq_sys.get_system<LinearImplicitSystem>("Linear");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const DofMap &dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(0);

    // Build a Finite Element object of the specified type.
    // Note that FEVectorBase is a typedef for the templated FE
    // class.
    std::unique_ptr<FEVectorBase> fe(FEVectorBase::build(dim, fe_type));

    // A 2*p+1 order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, fe_type.default_quadrature_order());

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    std::unique_ptr<FEVectorBase> fe_face(FEVectorBase::build(dim, fe_type));

    // Boundary integration requires one quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim - 1, fe_type.default_quadrature_order());

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    // Notice the shape functions are a vector rather than a scalar.
    const std::vector<std::vector<RealGradient>> &phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points. Notice that the shape function gradients are a tensor.
    const std::vector<std::vector<RealTensor>> &dphi = fe->get_dphi();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();

    vtkDataArray *jsArray = cell_data->GetArray("Js");
    if (!jsArray)
    {
        // Get the value from the "Js" array using the element ID
        libmesh_error_msg("Error: 'Js' array not found in vtkCellData.");
    }
    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
    // Element iterators are a nice way to iterate through all the
    // elements, or all the elements that have some property.  The
    // iterator el will iterate from the first to the last element on
    // the local processor.  The iterator end_el tells us when to stop.
    // It is smart to make this one const so that we don't accidentally
    // mess it up!  In case users later modify this program to include
    // refinement, we will be safe and will only consider the active
    // elements; hence we use a variant of the active_elem_iterator.
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        dof_id_type elem_id = elem->id();

        Real tuple[3];
        jsArray->GetTuple(elem_id, tuple);
        RealGradient js(4 * libMesh::pi * 1e-7 * tuple[0],
                        4 * libMesh::pi * 1e-7 * tuple[1],
                        4 * libMesh::pi * 1e-7 * tuple[2]);

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);
        const unsigned int n_dofs = cast_int<unsigned int>(dof_indices.size());
        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        Ke.resize(n_dofs, n_dofs);

        Fe.resize(n_dofs);

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test functions (i) against
            // the trial functions (j).
            for (std::size_t i = 0; i < phi.size(); i++)
                for (std::size_t j = 0; j < phi.size(); j++)
                    Ke(i, j) += JxW[qp] * dphi[i][qp].contract(dphi[j][qp]);

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.

            for (std::size_t i = 0; i < phi.size(); i++)
                Fe(i) += JxW[qp] * js * phi[i][qp];
        }

        // We have now finished the quadrature point loop,
        // and have therefore applied all the boundary conditions.

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.heterogenously_constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        matrix.add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }
}

void SystemBuilder::createSystem()
{
    transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
    ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
    fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
    formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
    order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
    dim = equation_systems.get_mesh().mesh_dimension();

    if (transient)
    {
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
    }
    else
    {
        if (ferromagnetic)
        {
            system_name = "Nonlinear";
            // Add NonlinearImplicitSystem as System and return its pointer
            equation_systems.add_system<NonlinearImplicitSystem>(system_name);
        }
        else
        {
            system_name = "Linear";
            // Add LinearImplicitSystem as System and return its pointer
            equation_systems.add_system<LinearImplicitSystem>(system_name);
        }
        // Assume the rest of the method continues to handle other cases based on the attributes
    }
};

void SystemBuilder::setVariables()
{
    System &system = equation_systems.get_system(system_name);
    std::string fe_fam = fe_family;
    if (fe_family == "LAGRANGE" && dim == 3)
    {
        fe_fam = "LAGRANGE_VEC";
    }

    // Example combinations - replace with actual logic and variable names
    if (fe_family == "LAGRANGE" && dim == 2 && formulation == "AV")
    {
        system.add_variable("AZ", Utility::string_to_enum<Order>(order),
                            Utility::string_to_enum<FEFamily>(fe_family));
        if (system_name == "NonlinearTransient" || system_name == "LinearTransient")
        {
            system.add_variable("V", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_family));
        }
        AssemblerStatic2D *assembler = new AssemblerStatic2D(equation_systems, cell_data);
        system.attach_assemble_object(*assembler);
    }
    else if (fe_family == "LAGRANGE" && dim == 3 && formulation == "AV")
    {
        system.add_variable("A", Utility::string_to_enum<Order>(order),
                            Utility::string_to_enum<FEFamily>(fe_fam));
        // system.add_variable("AY", Utility::string_to_enum<Order>(order),
        //                     Utility::string_to_enum<FEFamily>(fe_family));
        // system.add_variable("AZ", Utility::string_to_enum<Order>(order),
        //                     Utility::string_to_enum<FEFamily>(fe_family));
        if (system_name == "NonlinearTransient" || system_name == "LinearTransient")
        {
            system.add_variable("V", Utility::string_to_enum<Order>(order),
                                Utility::string_to_enum<FEFamily>(fe_fam));
        }
        AssemblerStatic3D *assembler = new AssemblerStatic3D(equation_systems, cell_data);
        system.attach_assemble_object(*assembler);
    }
    // else if (fe_family == "family1" && dim == 3 && formulation == "form1")
    // {
    //     system.add_variable("variable3", SECOND, LAGRANGE);
    // }
    // else if (fe_family == "family1" && dim == 3 && formulation == "form2")
    // {
    //     system.add_variable("variable4", SECOND, LAGRANGE);
    // }
    // else if (fe_family == "family2" && dim == 2 && formulation == "form1")
    // {
    //     system.add_variable("variable5", FIRST, MONOMIAL);
    // }
    // else if (fe_family == "family2" && dim == 3 && formulation == "form2")
    // {
    //     system.add_variable("variable6", SECOND, MONOMIAL);
    // }
    else
    {
        libmesh_error_msg("not implemented");
    }
    // Add more conditions as needed for other combinations
}

// System &SystemBuilder::getSystemByType()
// {

//     if (system_name == "NonlinearTransient")
//     {
//         return equation_systems.get_system<TransientSystem<NonlinearImplicitSystem>>(system_name);
//     }
//     else if (system_name == "LinearTransient")
//     {
//         return equation_systems.get_system<TransientSystem<LinearImplicitSystem>>(system_name);
//     }
//     else if (system_name == "Nonlinear")
//     {
//         return equation_systems.get_system<NonlinearImplicitSystem>(system_name);
//     }
//     else if (system_name == "Linear")
//     {
//         return equation_systems.get_system<LinearImplicitSystem>(system_name);
//     }
//     return nullptr;
// }
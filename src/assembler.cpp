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
#include "../include/element_builder.h" // Include the header file for ElementBuilder
#include <map>                          // Include for std::map
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
    // for the first variable in the system.
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
    DenseMatrix<Number> Ke;
    DenseMatrix<Number> Ce;
    DenseVector<Number> Fe;

    const TimeStep &ts = time_stepper.get_time_steps()[0];

    ElementBuilder elem_builder(Ke, Ce, Fe, materials, fe, cell_data, ts, dof_map, qrule);

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();

    // You will have to fix this for transient problems
    // std::string array_name = "Js_" + std::to_string(time_step);
    // vtkDataArray *jsArray = cell_data->GetArray(array_name.c_str());
    // array_name = "ansys_material_type";
    // vtkDataArray *matArray = cell_data->GetArray(array_name.c_str());
    // if (!jsArray)
    // {
    //     // Get the value from the "Js" array using the element ID
    //     libmesh_error_msg("Error: Js_" + std::to_string(time_step) + " array not found in vtkCellData.");
    // }
    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.  See example 3 for a discussion of the
    // element iterators.
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.

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

        dof_id_type elem_id = elem->id();

        // Retrieve the material number for the current element
        //    const unsigned int material_number = static_cast<unsigned int>(matArray->GetTuple1(el_id)); // Retrieve the material number using GetTuple1

        // Access the corresponding material properties from materials
        // const MaterialsProperties &mat_properties = materials.at(material_number);

        // Create a diagonal matrix with murx, mury, and murz as the diagonal terms
        // std::vector<std::vector<double>> mu_matrix(3, std::vector<double>(3, 0.0));
        // mu_matrix[0][0] = mat_properties.murx;
        // mu_matrix[1][1] = mat_properties.mury;
        // mu_matrix[2][2] = mat_properties.murz;
        elem_builder.compute_material_properties(elem_id, subdomain_id);
        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        // libmesh_assert_equal_to(n_dofs, phi.size());

        Ke.resize(n_dofs, n_dofs);

        Ce.resize(n_dofs, n_dofs);

        Fe.resize(n_dofs);

        elem_builder.old_solution_values.resize(n_dofs);

        system.old_local_solution->get(dof_indices, elem_builder.old_solution_values);

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

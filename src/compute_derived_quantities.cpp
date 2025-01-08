// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/explicit_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/utility.h"
#include "libmesh/dense_matrix.h"
#include "../include/compute_derived_quantities.h"
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
#include <cmath>
#include <stdexcept>
using namespace libMesh;

// Post-process the solution to compute stresses
void AssemblerDerivedQuantitiesLinear::compute_derived_quantities(EquationSystems &es)
{
    {
        // Declare a performance log.  Give it a descriptive
        // string to identify what part of the code we are
        // logging, since there may be many PerfLogs in an
        // application.
        PerfLog perf_log("Matrix Assembly");
        // Get a constant reference to the mesh object.
        const MeshBase &mesh = es.get_mesh();

        // Get a reference to the LinearImplicitSystem we are solving
        TransientLinearImplicitSystem &system = es.get_system<TransientLinearImplicitSystem>("PlasmaArcLinearSystem");
        ExplicitSystem &derived_quantities = es.get_system<ExplicitSystem>("DerivedQuantitiesSystem");
        // A reference to the DofMap object for this system.  The DofMap
        // object handles the index translation from node and element numbers
        // to degree of freedom numbers.  We will talk more about the DofMap
        // in future examples.
        const DofMap &dof_map = system.get_dof_map();

        // Get a constant reference to the Finite Element type
        // for the first variable in the system.
        FEType fe_type_A = dof_map.variable_type(0);

        // Build a Finite Element object of the specified type.  Since the
        // FEBase::build() member dynamically creates memory we will
        // store the object as a std::unique_ptr<FEBase>.  This can be thought
        // of as a pointer that will clean up after itself.
        std::unique_ptr<FEVectorBase> fe_A(FEVectorBase::build(dim, fe_type_A));

        // A 2nd order Gauss quadrature rule for numerical integration.
        QGauss qrule(dim, fe_type_A.default_quadrature_order());

        // Tell the finite element object to use our quadrature rule.
        fe_A->attach_quadrature_rule(&qrule);

        // Get a constant reference to the Finite Element type
        // for the first variable in the system.
        FEType fe_type_V = dof_map.variable_type(1);

        // Build a Finite Element object of the specified type.  Since the
        // FEBase::build() member dynamically creates memory we will
        // store the object as a std::unique_ptr<FEBase>.  This can be thought
        // of as a pointer that will clean up after itself.
        std::unique_ptr<FEBase> fe_V(FEBase::build(dim, fe_type_V));

        // A 2nd order Gauss quadrature rule for numerical integration.
        QGauss qrule_V(dim, fe_type_V.default_quadrature_order());

        // Tell the finite element object to use our quadrature rule.
        fe_V->attach_quadrature_rule(&qrule_V);
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

        auto ts = eq_sys.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper");
        // check correct time stepper get
        // Get the cell data from the vtkGrid
        vtkCellData *cell_data = eq_sys.parameters.get<vtkCellData *>("cell_data");
        // check that cell data has been correctly passed
        if (cell_data == nullptr)
        {
            libmesh_error_msg("Error: cell data not found in vtkCellData.");
        }
        // Get the materials properties
        MaterialsMap *materials = eq_sys.parameters.get<MaterialsMap *>("materials_properties");
        const TimeStep &time_step = ts->get_time_steps()[ts->get_current_time_step()];
        ElementBuilder elem_builder(Ke, Ce, Fe, *materials, fe_A, fe_V, cell_data, time_step, qrule);

        // Allocate memory for 3D nodal elements matrixes with transient effects
        // I can have up to 4 variables per node/edge
        // References to momentum equation right hand sides

        // This vector will hold the degree of freedom indices for
        // the element.  These define where in the global system
        // the element degrees of freedom get mapped.
        std::vector<dof_id_type> dof_indices;

        // The global system matrix
        SparseMatrix<Number> &matrix = system.get_system_matrix();
        dof_id_type elem_id;
        subdomain_id_type subdomain_id;
        const Real relative_epsilon = 1e-9;
        Real dt1;
        Real dt2;

        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            // Start logging the shape function initialization.
            // This is done through a simple function call with
            // the name of the event to log.
            subdomain_id = elem->subdomain_id();
            elem_id = elem->id();
            unsigned int n_dofs;

            if (ts->get_cumulative_substep() == 1) // Here add the check for nonlinear elements
            {
                // this can be done at the beginning only
                dof_map.dof_indices(elem, dof_indices);
                this->insert_dof_indices(elem_id, dof_indices);

                // Cache the number of degrees of freedom on this element, for
                // use as a loop bound later.  We use cast_int to explicitly
                // convert from size() (which may be 64-bit) to unsigned int
                // (which may be 32-bit but which is definitely enough to count
                // *local* degrees of freedom.
                n_dofs = cast_int<unsigned int>(dof_indices.size());
            }
            else
            {
                dof_indices = this->get_dof_indices(elem_id);
                n_dofs = cast_int<unsigned int>(dof_indices.size());
            }
            // Get the degree of freedom indices for the
            // current element.  These define where in the global
            // matrix and right-hand-side this element will
            // contribute to.
            if (ts->get_cumulative_substep() == 1) // Here add the check for nonlinear elements
            {
                // this can be done at the beginning only
                dof_map.dof_indices(elem, dof_indices);
                this->insert_dof_indices(elem_id, dof_indices);

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
                fe_A->reinit(elem);
                fe_V->reinit(elem);

                elem_builder.compute_material_properties(elem_id, subdomain_id);
                // With one variable, we should have the same number of degrees
                // of freedom as shape functions.
                // libmesh_assert_equal_to(n_dofs, phi.size());

                Ke.resize(n_dofs, n_dofs);
                if (subdomain_id == 1)
                    Ce.resize(n_dofs, n_dofs);

                elem_builder.build_element_matrix(elem, dof_map, subdomain_id, formulation, dim, fe_family);

                this->insert_ke_matrix(elem_id, Ke);
                if (subdomain_id == 1)
                {
                    this->insert_ce_matrix(elem_id, Ce);
                    // Open a file stream for printing
                    std::ofstream ke_file("../simulations/single_element_cond/Ke_matrix.txt");
                    std::ofstream ce_file("../simulations/single_element_cond/Ce_matrix.txt");

                    if (ke_file.is_open() && ce_file.is_open())
                    {
                        Ke.print_scientific(ke_file);
                        Ce.print_scientific(ce_file);
                    }
                }
                if (subdomain_id == 0)
                {
                    Fe.resize(n_dofs);
                    elem_builder.build_element_rhs(elem, dof_map, subdomain_id, formulation, dim, fe_family);
                    this->insert_rhs_vector(elem_id, Fe);
                }
            }

            dt1 = ts->get_time_substep_dt();
            dt2 = ts->get_prev_time_substep_dt();
            // in a new time step we could have a different dt, different bcs, different time integration effects
            // Here add the flag for ferromagnetic
            if (std::abs(dt1 - dt2) > relative_epsilon * std::max(std::abs(dt1), std::abs(dt2)) || ts->get_current_substep() == 0)
            {
                if (elem_id == 0)
                {
                    // libMesh::out << "Time step changed from " << dt2 << " to " << dt1 << std::endl;
                    libMesh::out << "Computing new system matrix" << std::endl;
                }
                KKe.resize(n_dofs, n_dofs);
                KKe += this->get_ke_matrix(elem_id);

                if (subdomain_id == 1)
                {
                    if (time_step.get_transient_effect())
                    {
                        CCe.resize(n_dofs, n_dofs);
                        CCe += this->get_ce_matrix(elem_id);
                        CCe.scale(1 / time_step.get_substeps()[ts->get_current_substep()].get_dt());
                        KKe += CCe;
                        std::ofstream kke_file("../simulations/single_element_cond/Kke_matrix.txt");

                        if (kke_file.is_open())
                        {
                            Ke.print_scientific(kke_file);
                        }
                    }
                }
                dof_map.constrain_element_matrix(KKe, dof_indices);
                matrix.add_matrix(KKe, dof_indices);
            }

            FFe.resize(n_dofs);

            // dof_map.constrain_element_vector(FFe, dof_indices);
            // system.rhs->add_vector(FFe, dof_indices);
        }
        // here I constrain but I modify later, wrong
        // }
        if (subdomain_id == 1 && time_step.get_transient_effect())
        {
            system.old_local_solution->get(dof_indices, old_solution_values);
            this->get_ce_matrix(elem_id).vector_mult(FFe, old_solution_values);
            FFe.scale(1 / time_step.get_substeps()[ts->get_current_substep()].get_dt());
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
    }
    // That's it.  We don't need to do anything else to the
    // PerfLog.  When it goes out of scope (at this function return)
    // it will print its log to the screen. Pretty easy, huh?
}
const MeshBase &mesh = es.get_mesh();
const unsigned int dim = mesh.mesh_dimension();
LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>(0);

const unsigned int a_var = system.variable_number("AZ");

const DofMap &dof_map = system.get_dof_map();

FEType fe_type = dof_map.variable_type(a_var);

std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

QGauss qrule(dim, fe_type.default_quadrature_order());

fe->attach_quadrature_rule(&qrule);

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
const std::vector<Real> &dxidy = fe->get_dxidy();
const std::vector<Real> &detady = fe->get_detady();
const std::vector<Real> &dxidx = fe->get_dxidx();
const std::vector<Real> &detadx = fe->get_detadx();

const std::vector<std::vector<Real>> &dphidxi = fe->get_dphidxi();
const std::vector<std::vector<Real>> &dphideta = fe->get_dphideta();

// Also, get a reference to the ExplicitSystem
magnetic_field_system.add_variable("BX", FIRST, L2_LAGRANGE);
magnetic_field_system.add_variable("BY", FIRST, L2_LAGRANGE);
magnetic_field_system.init();
unsigned int bx_n_var = magnetic_field_system.variable_number("BX");
unsigned int by_n_var = magnetic_field_system.variable_number("BY");
// magnetic_field_system.add_variable("BY", FIRST, L2_LAGRANGE);
const DofMap &B_dof_map = magnetic_field_system.get_dof_map();

// Storage for the stress dof indices on each element
std::vector<dof_id_type> Bx_indices_var;
std::vector<dof_id_type> By_indices_var;
std::vector<dof_id_type> A_indices_var;
DenseMatrix<Number> Kxe;
DenseVector<Number> Ae;
DenseVector<Number> Bxe_q;
DenseMatrix<Number> Kye;
DenseVector<Number> Bye_q;
DenseVector<Number> Bxe;
DenseVector<Number> Bye;
for (const auto &elem : mesh.active_local_element_ptr_range())
{
    dof_map.dof_indices(elem, A_indices_var);

    const unsigned int n_var_dofs = A_indices_var.size();
    const unsigned int n_dofs = cast_int<unsigned int>(n_var_dofs);

    fe->reinit(elem);

    Kxe.resize(n_dofs, n_dofs);
    Ae.resize(n_dofs);
    Bxe_q.resize(n_dofs);
    Kye.resize(n_dofs, n_dofs);
    Bye_q.resize(n_dofs);
    Bxe.resize(n_dofs);
    Bye.resize(n_dofs);
    for (std::size_t j = 0; j < n_dofs; j++)
        Ae(j) = system.current_solution(A_indices_var[j]);
    // Build a vector containing dxixy, detady
    std::vector<std::vector<Real>> jacob_derivatives_bx = {dxidy, detady};
    std::vector<std::vector<Real>> jacob_derivatives_by = {dxidx, detadx};
    std::vector<std::vector<std::vector<Real>>> phi_derivatives = {dphidxi, dphideta};
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        // Row is variable u1, u2, or u3, column is x, y, or z
        for (std::size_t i = 0; i < n_dofs; i++)
            for (std::size_t j = 0; j < n_dofs; j++)
                for (std::size_t l = 0; l < dim; l++)
                {
                    Kxe(i, j) += JxW[qp] * phi[i][qp] * (jacob_derivatives_bx[l][qp] * phi_derivatives[l][j][qp]);
                    Kye(i, j) -= JxW[qp] * phi[i][qp] * (jacob_derivatives_by[l][qp] * phi_derivatives[l][j][qp]);
                }
    }

    Kxe.vector_mult(Bxe_q, Ae);
    Kye.vector_mult(Bye_q, Ae);
    B_dof_map.dof_indices(elem, Bx_indices_var, bx_n_var);
    B_dof_map.dof_indices(elem, By_indices_var, by_n_var);

    // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
    // Note that this gives a discontinuous stress plot on element boundaries, which is
    // appropriate. We then also get the von Mises stress from the projected stress tensor.
    DenseMatrix<Real> Me_bx(n_dofs, n_dofs);
    DenseMatrix<Real> Me_by(n_dofs, n_dofs);
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
        for (unsigned int i = 0; i < n_dofs; i++)
            for (unsigned int j = 0; j < n_dofs; j++)
            {
                Me_bx(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                Me_by(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
            }
    }

    Me_bx.cholesky_solve(Bxe_q, Bxe);
    Me_by.cholesky_solve(Bye_q, Bye);

    dof_id_type dof_index;
    // Loop over B_indices_var_values
    for (unsigned int i = 0; i < n_dofs; i++)
    {
        dof_index = Bx_indices_var[i];
        Real val;
        if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
            (dof_index < magnetic_field_system.solution->last_local_index()))
            magnetic_field_system.solution->set(dof_index, Bxe(i));
        dof_index = By_indices_var[i];
        if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
            (dof_index < magnetic_field_system.solution->last_local_index()))
            magnetic_field_system.solution->set(dof_index, Bye(i));
    }
}
// Should call close and update when we set vector entries directly
magnetic_field_system.solution->close();
magnetic_field_system.update();
}

// Post-process the solution to compute stresses
void compute_magnetic_field_static_3D(EquationSystems &es)
{
    const MeshBase &mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    ExplicitSystem &magnetic_field_system = es.add_system<ExplicitSystem>("MagneticFieldSystem");
    LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>(0);

    const unsigned int a_var = system.variable_number("A");

    const DofMap &dof_map = system.get_dof_map();

    FEType fe_type = dof_map.variable_type(a_var);

    std::unique_ptr<FEVectorBase> fe(FEVectorBase::build(dim, fe_type));

    QGauss qrule(dim, fe_type.default_quadrature_order());

    fe->attach_quadrature_rule(&qrule);

    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &curl_phi = fe->get_curl_phi();

    // Also, get a reference to the ExplicitSystem
    magnetic_field_system.add_variable("B", FIRST, L2_LAGRANGE_VEC);
    magnetic_field_system.init();
    unsigned int b_n_var = magnetic_field_system.variable_number("B");
    const DofMap &B_dof_map = magnetic_field_system.get_dof_map();

    // Storage for the stress dof indices on each element
    std::vector<dof_id_type> B_indices_var;
    std::vector<dof_id_type> A_indices_var;
    DenseMatrix<Number> Ke;
    DenseVector<Number> Ae;
    DenseVector<Number> Be;
    DenseVector<Number> Be_q;
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices(elem, A_indices_var);

        const unsigned int n_var_dofs = A_indices_var.size();
        const unsigned int n_dofs = cast_int<unsigned int>(n_var_dofs);

        fe->reinit(elem);

        Ke.resize(n_dofs, n_dofs);
        Ae.resize(n_dofs);
        Be_q.resize(n_dofs);

        for (std::size_t j = 0; j < n_dofs; j++)
            Ae(j) = system.current_solution(A_indices_var[j]);

        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            // Row is variable u1, u2, or u3, column is x, y, or z
            for (std::size_t i = 0; i < n_dofs; i++)
                for (std::size_t j = 0; j < n_dofs; j++)
                    Ke(i, j) += JxW[qp] * (phi[i][qp] * curl_phi[j][qp]);
        }

        Ke.vector_mult(Be_q, Ae);
        B_dof_map.dof_indices(elem, B_indices_var);

        // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
        // Note that this gives a discontinuous stress plot on element boundaries, which is
        // appropriate. We then also get the von Mises stress from the projected stress tensor.
        DenseMatrix<Real> Me(n_dofs, n_dofs);
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < n_dofs; i++)
                for (unsigned int j = 0; j < n_dofs; j++)
                {
                    Me(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
                }
        }

        Me.cholesky_solve(Be_q, Be);

        dof_id_type dof_index;
        // Loop over B_indices_var_values
        for (unsigned int i = 0; i < n_dofs; i++)
        {
            dof_index = B_indices_var[i];
            if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
                (dof_index < magnetic_field_system.solution->last_local_index()))
                magnetic_field_system.solution->set(dof_index, Be(i));
        }
    }
    // Should call close and update when we set vector entries directly
    magnetic_field_system.solution->close();
    magnetic_field_system.update();
}

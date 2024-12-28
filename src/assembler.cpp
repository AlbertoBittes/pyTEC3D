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
    DenseMatrix<Number> KKe;
    DenseMatrix<Number> Ke;
    DenseMatrix<Number> Ce;
    DenseMatrix<Number> Ke_aux;
    DenseMatrix<Number> Ce_aux;

    DenseVector<Number> FFe;
    DenseVector<Number> Fe;
    std::vector<Real> old_solution_values;
    DenseSubVector<Number> F(Fe);

    auto ts = eq_sys.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper");
    vtkCellData *cell_data = eq_sys.parameters.get<vtkCellData *>("cell_data");
    MaterialsMap *materials = eq_sys.parameters.get<MaterialsMap *>("materials_properties");
    const TimeStep &time_step = ts->get_time_steps()[ts->get_current_time_step()];
    ElementBuilder elem_builder(Ke, Ce, Fe, *materials, fe, cell_data, time_step, qrule);

    // Allocate memory for 3D nodal elements matrixes with transient effects
    // I can have up to 4 variables per node/edge
    // References to momentum equation right hand sides

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Start logging the shape function initialization.
        // This is done through a simple function call with
        // the name of the event to log.
        const unsigned int subdomain_id = elem->subdomain_id();
        dof_id_type elem_id = elem->id();

        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        if (ts->get_cumulative_substep() == 0) // Here add the check for nonlinear elements
        {
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
            fe->reinit(elem);

            elem_builder.compute_material_properties(elem_id, subdomain_id);
            // With one variable, we should have the same number of degrees
            // of freedom as shape functions.
            // libmesh_assert_equal_to(n_dofs, phi.size());

            Ke.resize(n_dofs, n_dofs);
            if (elem->subdomain_id() == 1)
                Ce.resize(n_dofs, n_dofs);

            elem_builder.build_element_matrix(elem, dof_map, subdomain_id, formulation, dim, fe_family);

            this->insert_ke_matrix(elem_id, Ke);
            KKe.resize(n_dofs, n_dofs);
            KKe += Ke;
            if (subdomain_id == 1)
            {
                Ce.scale(1 / time_step.get_substeps()[ts->get_current_substep()].get_dt());
                this->insert_ce_matrix(elem_id, Ce);
                if (time_step.get_transient_effect())
                    KKe += Ce;
            }
            Fe.resize(n_dofs);
            elem_builder.build_element_rhs(elem, dof_map, subdomain_id, formulation, dim, fe_family);
            this->insert_rhs_vector(elem_id, Fe);
            system.old_local_solution->get(dof_indices, old_solution_values);
        }
        // Matrix changes when we have nonlinearities or different dt
        else if (ts->get_cumulative_substep() != 0 && ts->get_prev_time_substep_dt() != ts->get_time_step_dt())
        {
            if (subdomain_id == 1)
            {
                dof_indices = this->get_dof_indices(elem_id);
                const unsigned int n_dofs = dof_indices.size();
                KKe.resize(n_dofs, n_dofs);
                KKe += this->get_ke_matrix(elem_id);
                Ce_aux = this->get_ce_matrix(elem_id);
                Ce_aux.scale(ts->get_prev_time_substep_dt() / time_step.get_substeps()[ts->get_current_substep()].get_dt());
                if (time_step.get_transient_effect())
                    KKe += Ce_aux;
                this->insert_ce_matrix(elem_id, Ce_aux);
                // this->insert_rhs_vector(elem_id, Fe);
            }
        }
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

void AssemblerLinear::insert_ke_matrix(int key, const DenseMatrix<Number> &matrix)
{
    ke_matrix_map[key] = DenseMatrix<Number>(matrix); // Insert a copy of the matrix
}

DenseMatrix<Number> &AssemblerLinear::get_ke_matrix(int key)
{
    auto it = ke_matrix_map.find(key);
    if (it != ke_matrix_map.end())
    {
        return it->second;
    }
    else
    {
        throw std::out_of_range("Key not found in ke_matrix_map");
    }
}

void AssemblerLinear::insert_ce_matrix(int key, const DenseMatrix<Number> &matrix)
{
    ce_matrix_map[key] = DenseMatrix<Number>(matrix); // Insert a copy of the matrix
}

DenseMatrix<Number> &AssemblerLinear::get_ce_matrix(int key)
{
    auto it = ce_matrix_map.find(key);
    if (it != ce_matrix_map.end())
    {
        return it->second;
    }
    else
    {
        throw std::out_of_range("Key not found in ce_matrix_map");
    }
}

void AssemblerLinear::insert_rhs_vector(int key, const DenseVector<Number> &vector)
{
    rhs_map[key] = DenseVector<Number>(vector); // Insert a copy of the vector
}

DenseVector<Number> &AssemblerLinear::get_rhs_vector(int key)
{
    auto it = rhs_map.find(key);
    if (it != rhs_map.end())
    {
        return it->second;
    }
    else
    {
        throw std::out_of_range("Key not found in rhs_map");
    }
}

void AssemblerLinear::insert_old_rhs_vector(int key, const std::vector<Real> &vector)
{
    old_rhs_map[key] = std::vector<Real>(vector); // Insert a copy of the vector
}

std::vector<Real> &AssemblerLinear::get_old_rhs_vector(int key)
{
    auto it = old_rhs_map.find(key);
    if (it != old_rhs_map.end())
    {
        return it->second;
    }
    else
    {
        throw std::out_of_range("Key not found in old_rhs_map");
    }
}

void AssemblerLinear::insert_dof_indices(int key, const std::vector<dof_id_type> &dof_indices_param)
{
    dof_indices[key] = std::vector<dof_id_type>(dof_indices_param); // Insert a copy of the vector
}

std::vector<dof_id_type> &AssemblerLinear::get_dof_indices(int key)
{
    auto it = dof_indices.find(key);
    if (it != dof_indices.end())
    {
        return it->second;
    }
    else
    {
        throw std::out_of_range("Key not found in dof_indices");
    }
}
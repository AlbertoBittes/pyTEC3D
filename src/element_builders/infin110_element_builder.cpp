#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "loads_object.h"
#include "apec_m_system_objects.h"
#include "element_types_factory.h"

void ElementBuilder::build_infin110()
{
    // Stop logging the shape function initialization.
    // If you forget to stop logging an event the PerfLog
    // object will probably catch the error and abort.
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    // Now start logging the element matrix computation
    // Numeric ids corresponding to each variable in the system
    // Now we will build the element matrix.  This involves
    // a double loop to integrate the test functions (i) against
    // the trial functions (j).
    //
    // We have split the numeric integration into two loops
    // so that we can log the matrix and right-hand-side
    // computation separately.
    //
    INFIN110Type *inf_fe_type = dynamic_cast<INFIN110Type *>(elem_data.element_type);

    FEBase *inf_fe = dynamic_cast<FEBase *>(inf_fe_type->inf_fe_array[0].get());
    // The element Jacobian * quadrature weight at each integration point.

    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = inf_fe->get_JxWxdecay_sq();
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = inf_fe->get_phi_over_decayxR();
    const std::vector<std::vector<RealGradient>> &dphi = inf_fe->get_dphi_over_decayxR();
    const std::vector<Real> &sob_w = inf_fe->get_Sobolev_weightxR_sq();
    const std::vector<RealGradient> &dsob_w = inf_fe->get_Sobolev_dweightxR_sq();

    inf_fe->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    //   const unsigned int n_sf =
    // FEInterface::n_dofs(cfe->get_fe_type(), elem);
    unsigned int n_dofs_var = inf_fe->n_shape_functions();
    unsigned int nq = inf_fe->n_quadrature_points();

    auto &Ke = this->matrices[0][0];

    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        for (unsigned int i = 0; i < n_dofs_var; ++i)
        {
            // store the test functions derivative (which contains the extra sobolev weight)
            // separately; so we don't need to recompute it for each j.
            const RealGradient dphi_i = dphi[i][qp] * sob_w[qp] + phi[i][qp] * dsob_w[qp];
            for (unsigned int j = 0; j < n_dofs_var; ++j)
            {
                Ke(i, j) += JxW[qp] * dphi_i * dphi[j][qp];
            }
        }
    }
}

void ElementBuilder::build_infin110_rhs()
{
    return;
}

void ElementBuilder::build_infin110_nodal_solution_vector()
{
    return;
}

void ElementBuilder::build_infin110_element_solution_vector()
{
    return;
}
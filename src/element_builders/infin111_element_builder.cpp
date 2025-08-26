#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>
#include <libmesh/petsc_linear_solver.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "globals.h"
#include "apec_m_system_objects.h"
#include "element_types_factory.h"

// Not implemented placeholders for INFIN111 element builder methods

void ElementBuilder::build_infin111()
{
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[4].get());
    // The element Jacobian * quadrature weight at each integration point.

    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxWxdecay_sq();
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi_over_decayxR();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi_over_decayxR();
    const std::vector<Real> &sob_w = fe->get_Sobolev_weightxR_sq();
    const std::vector<RealGradient> &dsob_w = fe->get_Sobolev_dweightxR_sq();

    fe->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    //   const unsigned int n_sf =
    // FEInterface::n_dofs(cfe->get_fe_type(), elem);
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    auto &Ke = this->matrices[0][0];
    // unsigned int n_nodes = elem_data.elem->n_nodes();
    // unsigned int n_dof_check = FEInterface::n_dofs(fe->get_fe_type(), elem_data.elem);

    // n_dof_check = FEInterface::n_dofs(fe->get_fe_type(), 0, elem_data.elem);

    // n_dof_check = FEInterface::n_dofs_at_node(fe->get_fe_type(), elem_data.elem, 0);

    // n_dof_check = FEInterface::n_dofs_per_elem(fe->get_fe_type(), elem_data.elem);

    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        for (unsigned int i = 0; i < n_dofs_var; ++i)
        {
            // store the test functions derivative (which contains the extra sobolev weight)
            // separately; so we don't need to recompute it for each j.
            const RealGradient dphi_i = dphi[i][qp] * sob_w[qp] + phi[i][qp] * dsob_w[qp];
            for (unsigned int j = 0; j < n_dofs_var; ++j)
            {
                for (unsigned int k = 0; k < 3; ++k)
                    Ke(i + k * n_dofs_var, j + k * n_dofs_var) += JxW[qp] * dphi_i * dphi[j][qp];
            }
        }
    }
    Ke.scale(1 / Globals::MU_0);
    return;
}

void ElementBuilder::build_infin111_rhs()
{
    return;
}

void ElementBuilder::build_infin111_solution_vectors()
{
    return;
}
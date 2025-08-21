#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/boundary_info.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "apec_m_system_objects.h"
#include "element_types_factory.h"

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236()
{
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
    // Now start logging the element matrix computation
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[0].get());
    FEVectorBase *fe_vec = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[1].get());
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_vec->get_phi();

    const std::vector<std::vector<RealGradient>> &curl_phi_vec = fe_vec->get_curl_phi();

    // const std::vector<std::vector<RealVectorValue>> &div_phi_vec = fe_vec->get_div_phi();
    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = fe->n_quadrature_points();

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = materials_object->get_mu_matrix();

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][SOLID97DOFMatrices::Ke_SOLID97];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][SOLID97DOFMatrices::Ce_SOLID97];
    // auto &Me = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::Me_SOLID97];
    // auto &AAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::AAe_SOLID97];
    // auto &JAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JAe_SOLID97];
    // auto &JVe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JVe_SOLID97];
    // auto &BAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::BAe_SOLID97];
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
            {
                Ke(dof_i, dof_j) += JxW[qp] / mu_matrix[0] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
            }

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        Real CAV_value;
        materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = materials_object->get_rs_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {
                    JVe(qp + k * n_dofs_var, dof_i) += dphi[dof_i][qp](k);
                    JAe(qp + k * n_dofs_var, dof_i + k * n_dofs_var) += phi[dof_i][qp];
                }
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                    {
                        Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / rs_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        CAV_value = JxW[qp] / rs_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
                        Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) += CAV_value;
                        Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) += CAV_value;
                        Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) += JxW[qp] * dphi[dof_i][qp](l) / rs_matrix[l] * dphi[dof_j][qp](l);
                    }
            }
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i + 3 * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
    // The penalty values, for solution boundary trace and flux.

    // Here use the formula for quadrature
    std::vector<boundary_id_type> side_boundary_ids;
    FEBase *fe_face = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[2].get());

    // The following loops over the sides of the element.
    // If the element has no neighbor on a side then that
    // side MUST live on a boundary of the domain.
    Real penalty_base;
    char cont = 0;
    auto boundary_info = &this->get_parent_object()->get_mesh().get_boundary_info();

    for (auto s : elem_data.elem->side_index_range())
        if (elem_data.elem->neighbor_ptr(s) == nullptr)
        {
            boundary_info->boundary_ids(elem_data.elem, s, side_boundary_ids);
            if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), Globals::FLUX_PARALLEL_ID_BC) != side_boundary_ids.end())
            {
                if (cont == 0)
                {
                    penalty_base = 1e7 * Ke.max() + 1e7 * Ce.max();
                    cont++;
                }
                // The value of the shape functions at the quadrature
                // points.
                const std::vector<std::vector<Real>> &phi_face =
                    fe_face->get_phi();

                // The value of the shape function derivatives at the
                // quadrature points.
                const std::vector<std::vector<RealGradient>> &dphi_face =
                    fe_face->get_dphi();

                // The Jacobian * Quadrature Weight at the quadrature
                // points on the face.
                const std::vector<Real> &JxW_face = fe_face->get_JxW();

                // The XYZ locations (in physical space) of the
                // quadrature points on the face.  This is where
                // we will interpolate the boundary value function.
                const std::vector<Point> &qface_point = fe_face->get_xyz();
                const std::vector<std::vector<Point>> &face_tangents = fe_face->get_tangents();

                // Compute the shape function values on the element
                // face.
                fe_face->reinit(elem_data.elem, s);

                // Loop over the face quadrature points for integration.
                for (unsigned int qp = 0; qp < fe_face->n_quadrature_points(); qp++)
                {
                    // Matrix contribution of the L2 projection.
                    // Note that the basis function values are
                    // integrated against test function values while
                    // basis fluxes are integrated against test function
                    // fluxes.
                    for (unsigned int i = 0; i < n_dofs_var; i++)
                        for (unsigned int j = 0; j < n_dofs_var; j++)
                            for (unsigned int l = 0; l < 3; l++)
                                for (unsigned int k = 0; k < 3; k++)
                                    for (unsigned int h = 0; h < face_tangents[0].size(); h++)
                                        Ke(i + k * n_dofs_var, j + l * n_dofs_var) += JxW_face[qp] *
                                                                                      (penalty_base * phi_face[i][qp] *
                                                                                       phi_face[j][qp] * face_tangents[qp][h](k) *
                                                                                       face_tangents[qp][h](l));
                }
            }
            if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), Globals::FLUX_NORMAL_ID_BC) != side_boundary_ids.end())
            {
                // The value of the shape functions at the quadrature
                // points.
                const std::vector<std::vector<Real>> &phi_face =
                    fe_face->get_phi();

                // The value of the shape function derivatives at the
                // quadrature points.
                const std::vector<std::vector<RealGradient>> &dphi_face =
                    fe_face->get_dphi();

                // The Jacobian * Quadrature Weight at the quadrature
                // points on the face.
                const std::vector<Real> &JxW_face = fe_face->get_JxW();

                // The XYZ locations (in physical space) of the
                // quadrature points on the face.  This is where
                // we will interpolate the boundary value function.
                const std::vector<Point> &qface_point = fe_face->get_xyz();
                const std::vector<Point> &face_normals = fe_face->get_normals();

                // Compute the shape function values on the element
                // face.
                fe_face->reinit(elem_data.elem, s);

                // Loop over the face quadrature points for integration.
                for (unsigned int qp = 0; qp < fe_face->n_quadrature_points(); qp++)
                {
                    // Matrix contribution of the L2 projection.
                    // Note that the basis function values are
                    // integrated against test function values while
                    // basis fluxes are integrated against test function
                    // fluxes.
                    // Compute the shape function values on the element
                    // face.
                    fe_face->reinit(elem_data.elem, s);

                    // Loop over the face quadrature points for integration.
                    for (unsigned int qp = 0; qp < fe_face->n_quadrature_points(); qp++)
                    {
                        // Matrix contribution of the L2 projection.
                        // Note that the basis function values are
                        // integrated against test function values while
                        // basis fluxes are integrated against test function
                        // fluxes.
                        for (unsigned int i = 0; i < n_dofs_var; i++)
                            for (unsigned int j = 0; j < n_dofs_var; j++)
                                for (unsigned int l = 0; l < 3; l++)
                                    for (unsigned int k = 0; k < 3; k++)
                                        Ke(i + k * n_dofs_var, j + l * n_dofs_var) += JxW_face[qp] *
                                                                                      (penalty_base * phi_face[i][qp] *
                                                                                       phi_face[j][qp] * face_normals[qp](k) *
                                                                                       face_normals[qp](l));
                    }
                }
            }
        }
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236_rhs()
{
    return;
}

void ElementBuilder::build_solid236_solution_vectors()
{
    return;
}

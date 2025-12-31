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
#include "apec_m_system.h"
#include <libmesh/boundary_info.h>
#include "element_types_factory.h"

#include "variables_map.h"

// {0: "Ke", 1: "Ce", 2: "Me", 3: "AAe", 4: "JAe", 5: "JVe"}
void ElementBuilder::build_solid97()
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
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[2].get());
    FEVectorBase *fe_vec = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[0].get());
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

    const std::vector<std::vector<Real>> &div_phi_vec = fe_vec->get_div_phi();

    const std::vector<std::vector<RealTensor>> &dphi_vec = fe_vec->get_dphi();
    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = fe->n_quadrature_points();

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_nu_matrix(this->elem_data.material_num);
    auto &nu_matrix = materials_object->get_nu_matrix();

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ke];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ce];
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::Me];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::AAe];
    auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];
    auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JVe];
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::BAe];
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                for (unsigned int l = 0; l < 3; l++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) +=
                            JxW[qp] * nu_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int l = 0; l < 3; l++)
                    BAe(qp + k * nq, dof_i + l * n_dofs_var) += curl_phi_vec[dof_i * 3 + l][qp](k);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int l = 0; l < 3; l++)
                        AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                            JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int l = 0; l < 3; l++)
                        Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                            JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                JAe(qp + k * nq, dof_i + k * n_dofs_var) += phi[dof_i][qp];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                JVe(qp + k * nq, dof_i) += dphi[dof_i][qp](k);

    // for (unsigned int qp = 0; qp < nq; qp++)
    //     for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
    //         for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
    //             for (unsigned int k = 0; k < 3; k++)
    //                 Ke(dof_i, dof_j) += JxW[qp] * nu_matrix[k] * (curl_phi_vec[dof_i / 3 + (dof_i % 3) * n_dofs_var][qp](k) * curl_phi_vec[dof_j / 3 + (dof_j % 3) * n_dofs_var][qp](k));

    // Real avg_nu = (nu_matrix[0] + nu_matrix[1] + nu_matrix[2]) / 3.0;
    // for (unsigned int qp = 0; qp < nq; qp++)
    //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
    //         for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
    //             for (unsigned int k = 0; k < 3; k++)
    //                 for (unsigned int l = 0; l < 3; l++)
    //                 {
    //                     Ke(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * nu_matrix[k] * (curl_phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
    //                     // for (unsigned int m = 0; m < 3; m++)
    //                     //     Ke(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * nu_matrix[m] * dphi_vec[dof_i * 3 + k][qp](m, m) * dphi_vec[dof_j * 3 + l][qp](m, m);
    //                     Ke(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * avg_nu * (div_phi_vec[dof_i * 3 + k][qp] * div_phi_vec[dof_j * 3 + l][qp]);
    //                 }

    // switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
    // {
    // case PhysicalRegionEnum::NON_MAGNETIC:
    // case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    // {
    //     Real avg_nu = (nu_matrix[0] + nu_matrix[1] + nu_matrix[2]) / 3.0;
    //     for (unsigned int qp = 0; qp < nq; qp++)
    //         for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
    //             for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
    //                 Ke(dof_i, dof_j) += JxW[qp] * avg_nu * (div_phi_vec[dof_i / 3 + (dof_i % 3) * n_dofs_var][qp] * div_phi_vec[dof_j / 3 + (dof_j % 3) * n_dofs_var][qp]);

    //     break;
    // }
    // }

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
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
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        auto real_constants_object = this->context<RealConstantsObject>();
        if (!real_constants_object || !real_constants_object->has_constant(this->elem_data.real_constant_num))
            libmesh_error_msg("Missing real constants for STRANDED_COIL elements.");

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real vc = rc[3];
        const Real tx = rc[5];
        const Real ty = rc[6];
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL elements.");

        const unsigned int n_dofs_A = n_dofs_var * 3;
        if (n_dofs >= n_dofs_A + 2 * n_dofs_var)
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        const Real kav_x = -(nc / (sc * r_coil)) * tx * nn;
                        const Real kav_y = -(nc / (sc * r_coil)) * ty * nn;
                        const Real kav_z = -(nc / (sc * r_coil)) * tz * nn;
                        const Real kvv = (1.0 / (r_coil * vc)) * nn;

                        Ce(dof_i + 0 * n_dofs_var, dof_j + n_dofs_A) += kav_x;
                        Ce(dof_i + 1 * n_dofs_var, dof_j + n_dofs_A) += kav_y;
                        Ce(dof_i + 2 * n_dofs_var, dof_j + n_dofs_A) += kav_z;
                        Ce(dof_i + 0 * n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kav_x;
                        Ce(dof_i + 1 * n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kav_y;
                        Ce(dof_i + 2 * n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kav_z;

                        Ce(dof_i + n_dofs_A, dof_j + 0 * n_dofs_var) += sym * kav_x;
                        Ce(dof_i + n_dofs_A, dof_j + 1 * n_dofs_var) += sym * kav_y;
                        Ce(dof_i + n_dofs_A, dof_j + 2 * n_dofs_var) += sym * kav_z;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + 0 * n_dofs_var) -= kav_x;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + 1 * n_dofs_var) -= kav_y;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + 2 * n_dofs_var) -= kav_z;

                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kvv;
                    }
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        auto real_constants_object = this->context<RealConstantsObject>();
        if (!real_constants_object || !real_constants_object->has_constant(this->elem_data.real_constant_num))
            libmesh_error_msg("Missing real constants for STRANDED_COIL_SOLENOIDAL elements.");

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real vc = rc[3];
        const Real tx = rc[5];
        const Real ty = rc[6];
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL_SOLENOIDAL elements.");

        const unsigned int n_dofs_A = n_dofs_var * 3;
        if (n_dofs >= n_dofs_A + 2 * n_dofs_var)
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        const Real kav_x = -(nc / (sc * r_coil)) * tx * nn;
                        const Real kav_y = -(nc / (sc * r_coil)) * ty * nn;
                        const Real kav_z = -(nc / (sc * r_coil)) * tz * nn;
                        const Real kvv = (1.0 / (r_coil * vc)) * nn;

                        Ke(dof_i + 0 * n_dofs_var, dof_j + n_dofs_A) += kav_x;
                        Ke(dof_i + 1 * n_dofs_var, dof_j + n_dofs_A) += kav_y;
                        Ke(dof_i + 2 * n_dofs_var, dof_j + n_dofs_A) += kav_z;
                        Ke(dof_i + 0 * n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kav_x;
                        Ke(dof_i + 1 * n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kav_y;
                        Ke(dof_i + 2 * n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kav_z;

                        Ce(dof_i + n_dofs_A, dof_j + 0 * n_dofs_var) += sym * kav_x;
                        Ce(dof_i + n_dofs_A, dof_j + 1 * n_dofs_var) += sym * kav_y;
                        Ce(dof_i + n_dofs_A, dof_j + 2 * n_dofs_var) += sym * kav_z;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + 0 * n_dofs_var) -= kav_x;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + 1 * n_dofs_var) -= kav_y;
                        Ce(dof_i + n_dofs_A + n_dofs_var, dof_j + 2 * n_dofs_var) -= kav_z;

                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A + n_dofs_var) -= kvv;
                        Ke(dof_i + n_dofs_A + n_dofs_var, dof_j + n_dofs_A) -= kvv;
                        Ke(dof_i + n_dofs_A + n_dofs_var, dof_j + n_dofs_A + n_dofs_var) += kvv;
                    }
        }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) +=
                            JxW[qp] * sigma_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) +=
                            JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) +=
                            JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) +=
                            JxW[qp] * dphi[dof_i][qp](l) * sigma_matrix[l] * dphi[dof_j][qp](l);
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i + 3 * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    JVe(qp + k * nq, dof_i) += dphi[dof_i][qp](k);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    JAe(qp + k * nq, dof_i + k * n_dofs_var) += phi[dof_i][qp];

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) +=
                            JxW[qp] * sigma_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ke(dof_i + l * n_dofs_var, n_dofs_A + dof_j) +=
                            JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) +=
                            JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Ke(n_dofs_A + dof_i, n_dofs_A + dof_j) +=
                            JxW[qp] * dphi[dof_i][qp](l) * sigma_matrix[l] * dphi[dof_j][qp](l);
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
    // The penalty values, for solution boundary trace and flux.

    // Here use the formula for quadrature
    std::vector<boundary_id_type> side_boundary_ids;
    FEBase *fe_face = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[1].get());

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
            if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), 32760) != side_boundary_ids.end())
            {
                if (cont == 0)
                {
                    penalty_base = 1e20 * Ke.max() + 1e20 * Ce.max();
                    cont++;
                }
                // The value of the shape functions at the quadrature
                // points.
                const std::vector<std::vector<Real>> &phi_face =
                    fe_face->get_phi();

                // The value of the shape function derivatives at the
                // quadrature points.
                // const std::vector<std::vector<RealGradient>> &dphi_face =
                //     fe_face->get_dphi();

                // The Jacobian * Quadrature Weight at the quadrature
                // points on the face.
                const std::vector<Real> &JxW_face = fe_face->get_JxW();

                // The XYZ locations (in physical space) of the
                // quadrature points on the face.  This is where
                // we will interpolate the boundary value function.
                // const std::vector<Point> &qface_point = fe_face->get_xyz();
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
            if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), 32761) != side_boundary_ids.end())
            {
                if (cont == 0)
                {
                    penalty_base = 1e20 * Ke.max() + 1e20 * Ce.max();
                    cont++;
                }
                // The value of the shape functions at the quadrature
                // points.
                const std::vector<std::vector<Real>> &phi_face =
                    fe_face->get_phi();

                // The value of the shape function derivatives at the
                // quadrature points.
                // const std::vector<std::vector<RealGradient>> &dphi_face =
                //     fe_face->get_dphi();

                // The Jacobian * Quadrature Weight at the quadrature
                // points on the face.
                const std::vector<Real> &JxW_face = fe_face->get_JxW();

                // The XYZ locations (in physical space) of the
                // quadrature points on the face.  This is where
                // we will interpolate the boundary value function.
                // const std::vector<Point> &qface_point = fe_face->get_xyz();
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

void ElementBuilder::build_solid97_rhs()
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
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    unsigned int n_dofs_var = elem_data.elem->n_nodes();

    std::vector<dof_id_type> dof_indices_loads;
    std::vector<Real> loads_vec;

    System *main_system = &this->get_parent_object()->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    System *nodal_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::NODAL_LOADS);
    System *element_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_LOADS);

    DofMap *dof_map_nodal_loads = &nodal_loads_system->get_dof_map();
    DofMap *dof_map_element_loads = &element_loads_system->get_dof_map();

    NumericVector<Number> *nodal_loads_vector = nodal_loads_system->request_vector(2);
    NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        // Pay attention here if you add other variables to elemental loads system
        dof_map_element_loads->dof_indices(elem_data.elem, dof_indices_loads);
        element_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i + k * n_dofs_var) *= loads_vec[k];

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    case PhysicalRegionEnum::STRANDED_COIL:
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        std::string volt_var_name;
        if (PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL == static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]) ||
            PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL == static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
            volt_var_name = "V";
        else
            volt_var_name = "VOLT";
        if (nodal_loads_system->has_variable("AMPS") && main_system->has_variable(volt_var_name))
        {
            unsigned int amps_var_num = nodal_loads_system->variable_number("AMPS");
            unsigned int volt_var_num = main_system->variable_number(volt_var_name);

            unsigned int proc_id = nodal_loads_system->comm().rank();
            dof_indices_loads.clear();

            for (const auto &node : elem_data.elem->node_ref_range())
                if (node.processor_id() == proc_id)
                    if (node.n_dofs(ApecMSystemBase::NODAL_LOADS, amps_var_num))
                        dof_indices_loads.push_back(node.dof_number(ApecMSystemBase::NODAL_LOADS, amps_var_num, 0));

            // Apply nodal AMPS loads
            nodal_loads_vector->get(dof_indices_loads, loads_vec);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                const dof_id_type dof_id = dof_indices_loads[dof_i];
                if (!nodal_loads_dof_cache.insert(dof_id).second)
                    continue;
                Fe(dof_i + volt_var_num * n_dofs_var) = loads_vec[dof_i];
            }
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    return;
}

void ElementBuilder::build_solid97_solution_vectors()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    std::vector<dof_id_type> dof_indices;

    auto materials_object = this->context<MaterialsObject>();
    materials_object->fill_nu_matrix(this->elem_data.material_num);
    auto &nu_matrix = materials_object->get_nu_matrix();

    auto real_constants_object = this->context<RealConstantsObject>();

    ApecMLinearSystem *main_sys = dynamic_cast<ApecMLinearSystem *>(&this->get_parent_object()->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    System *nodal_solution_system = &this->get_parent_object()->get_system(ApecMSystemBase::NODAL_SOLUTION);
    System *element_solution_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_SOLUTION);

    this->fill_solution_vector_support(main_sys->current_local_solution.get(), ApecMSystemBase::DOF_SOLUTION_SYSTEM);

    // Get DOF Solutions in the nodes
    DenseVector<Number> An;
    for (unsigned int k = 0; k < 3; k++)
    {
        An.append(solution_vector_support.sol_vectors[k]);
    }

    // Nodal solution
    if (this->compute_nodal_solution_vectors)
    {
        unsigned int n_dofs_var = elem_data.elem->n_nodes();

        auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::Me];
        auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::AAe];
        this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS].resize(NodalSolutionVectors::NUM_NODAL_SOLUTION_VECTORS);
        auto &Bn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::B];
        auto &Hn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::H];
        // Compute B in nodes
        DenseVector<Number> Bn_x;
        AAe.vector_mult(Bn_x, An);
        Me.cholesky_solve(Bn_x, Bn);

        // Compute H in nodes
        Hn.resize(Bn.size());
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Hn(dof_i + k * n_dofs_var) += Bn(dof_i + k * n_dofs_var) * nu_matrix[k];
    }

    // Elemental solution
    this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS].resize(ElementSolutionVectors::NUM_ELEMENT_SOLUTION_VECTORS);
    auto &Be = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::Be];
    auto &He = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::He];
    auto &Je = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::Je];
    Be.resize(3);
    He.resize(3);
    Je.resize(3);

    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::BAe];
    unsigned int n_qp = BAe.m() / 3;

    // Compute B in elements
    DenseVector<Number> Bq;
    BAe.vector_mult(Bq, An);
    for (unsigned int k = 0; k < 3; k++)
    {
        for (unsigned int qp = 0; qp < n_qp; qp++)
            Be(k) += Bq(qp + k * n_qp);
        Be(k) = Be(k) / n_qp;
    }

    // Compute H in elements
    for (unsigned int k = 0; k < 3; k++)
        He(k) = Be(k) * nu_matrix[k];

    // Compute J in elements
    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]);
    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        System *element_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_LOADS);
        DofMap *dof_map_element_loads = &element_loads_system->get_dof_map();
        dof_map_element_loads->dof_indices(elem_data.elem, dof_indices);
        NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);
        element_loads_vector->get(dof_indices, Je.get_values());
        // DenseVector<Number> Jq;
        // BAe.vector_mult(Jq, Hn);
        // for (unsigned int k = 0; k < 3; k++)
        // {
        //     for (unsigned int qp = 0; qp < n_qp; qp++)
        //         Je(k) += Jq(qp + k * n_qp);
        //     Je(k) = Je(k) / n_qp;
        // }

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];
        auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JVe];

        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        DenseVector<Number> Vn;

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
            Vn = solution_vector_support.sol_vectors[3];

        fill_solution_vector_support(main_sys->request_vector("solution_time_derivative"), 0); // solution_derivative

        // Get derivative of A DOF Solutions (dA/dt, dV/dt)
        DenseVector<Number> An_dt;
        for (unsigned int k = 0; k < 3; k++)
        {
            An_dt.append(solution_vector_support.sol_vectors[k]);
        }

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
            Vn = solution_vector_support.sol_vectors[3];

        // Compute dA/dt in quadrature points
        DenseVector<Number> Jq_A;
        JAe.vector_mult(Jq_A, An_dt);

        // Compute dV/dt in quadrature points
        DenseVector<Number> Jq_V;
        JVe.vector_mult(Jq_V, Vn);

        // sum the two contributions
        Jq_A += Jq_V;

        // average on quadrature points and scale by conductivity
        for (unsigned int k = 0; k < 3; k++)
        {
            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(k) -= Jq_A(qp + k * n_qp);
            Je(k) = Je(k) * sigma_matrix[k] / n_qp;
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        Real Vn = 0.0;
        Real En = 0.0;

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
        {
            Vn = solution_vector_support.sol_vectors[3][0];
            En = solution_vector_support.sol_vectors[4][0];
        }

        fill_solution_vector_support(main_sys->request_vector("solution_time_derivative"), ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
        {
            Vn = solution_vector_support.sol_vectors[3][0];
            En = solution_vector_support.sol_vectors[4][0];
        }

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];

        const VectorValue t_vec(rc[5], rc[6], rc[7]);
        const Real r_coil = rc[8];

        const Real coeff = nc / (sc * r_coil) * (Vn - En);
        if (sc == 0.0 || r_coil == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R) for STRANDED_COIL_SOLENOIDAL elements.");

        for (unsigned int k = 0; k < 3; k++)
            Je(k) = coeff * t_vec(k);
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    if (this->compute_nodal_solution_vectors)
    {
        // Set B and H in nodes
        DofMap *dof_map_nodal_solution = &nodal_solution_system->get_dof_map();
        dof_map_nodal_solution->dof_indices(elem_data.elem, dof_indices);
        auto &Bn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::B];
        auto &Hn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::H];
        Bn.append(Hn);
        for (unsigned int i = 0; i < dof_indices.size(); i++)
        {
            dof_id_type dof_index = dof_indices[i];
            if ((nodal_solution_system->solution->first_local_index() <= dof_index) &&
                (dof_index < nodal_solution_system->solution->last_local_index()))
                nodal_solution_system->solution->set(dof_index, Bn(i));
        }
    }
    DofMap *dof_map_element_solution = &element_solution_system->get_dof_map();
    dof_map_element_solution->dof_indices(elem_data.elem, dof_indices);
    Je.append(Be);
    Je.append(He);
    for (unsigned int i = 0; i < dof_indices.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((element_solution_system->solution->first_local_index() <= dof_index) &&
            (dof_index < element_solution_system->solution->last_local_index()))
            element_solution_system->solution->set(dof_index, Je(i));
    }
}

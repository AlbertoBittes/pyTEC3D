#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

// Define Gauss quadrature rules.
#include "element_types.h" // Include the header defining ElementType
#include "apec_m_system.h"
#include "element_types_factory.h"
#include "variables_map.h"

namespace
{
    constexpr unsigned int kNodalMe = 0;
    constexpr unsigned int kNodalAAe = 1;
    constexpr unsigned int kElementJAe = 0;
    constexpr unsigned int kElementJVe = 1;
    constexpr unsigned int kElementBAe = 2;
}

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
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[1].get());
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

    // const std::vector<std::vector<Real>> &div_phi_vec = fe_vec->get_div_phi();

    const std::vector<std::vector<RealTensor>> &dphi_vec = fe_vec->get_dphi();
    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = fe->n_quadrature_points();

    const auto nu_matrix = this->get_mu_matrix(this->elem_data.material_num);

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ke];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ce];
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
    auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    auto build_ke = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) +=
                                JxW[qp] * nu_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
    };

    auto build_bae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int l = 0; l < 3; l++)
                        BAe(qp + k * nq, dof_i + l * n_dofs_var) += curl_phi_vec[dof_i * 3 + l][qp](k);
    };

    auto build_aae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        for (unsigned int l = 0; l < 3; l++)
                            AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                                JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
    };

    auto build_me = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        for (unsigned int l = 0; l < 3; l++)
                            Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                                JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);
    };

    auto build_jae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    JAe(qp + k * nq, dof_i + k * n_dofs_var) += phi[dof_i][qp];
    };

    auto build_jve = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    JVe(qp + k * nq, dof_i) += dphi[dof_i][qp](k);
    };

    build_ke();
    build_bae();
    build_aae();
    build_me();
    build_jae();
    build_jve();

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

    // switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]))
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

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto build_fa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int k = 0; k < 3; k++)
                        Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        };
        build_fa();
        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        const auto sigma_matrix = this->get_sigma_matrix(this->elem_data.material_num);

        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) +=
                                JxW[qp] * sigma_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);
        };

        auto build_cav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) +=
                                JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) +=
                                JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
        };

        auto build_cvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) +=
                                JxW[qp] * dphi[dof_i][qp](l) * sigma_matrix[l] * dphi[dof_j][qp](l);
        };

        build_caa();
        build_cav();
        build_cva();
        build_cvv();
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i + 3 * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        const auto sigma_matrix = this->get_sigma_matrix(this->elem_data.material_num);

        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) +=
                                JxW[qp] * sigma_matrix[l] * (phi[dof_i][qp] * phi[dof_j][qp]);
        };

        auto build_kav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ke(dof_i + l * n_dofs_var, n_dofs_A + dof_j) +=
                                JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) +=
                                JxW[qp] * sigma_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
        };

        auto build_kvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Ke(n_dofs_A + dof_i, n_dofs_A + dof_j) +=
                                JxW[qp] * dphi[dof_i][qp](l) * sigma_matrix[l] * dphi[dof_j][qp](l);
        };

        build_caa();
        build_kav();
        build_cva();
        build_kvv();
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        const auto rc = this->get_real_constants(this->elem_data.real_constant_num);
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

        auto build_coil_coupling = [&]()
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
        };
        build_coil_coupling();
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        const auto rc = this->get_real_constants(this->elem_data.real_constant_num);
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

        auto build_coil_coupling = [&]()
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
                        Ke(dof_i + n_dofs_A + n_dofs_var, dof_j + n_dofs_A + n_dofs_var) -= kvv;
                    }
        };
        build_coil_coupling();
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
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

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
        const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
        Real js[3] = {
            js_array->GetComponent(elem_id, 0),
            js_array->GetComponent(elem_id, 1),
            js_array->GetComponent(elem_id, 2)};

        auto build_fa = [&]()
        {
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    Fe(dof_i + k * n_dofs_var) *= js[k];
        };
        build_fa();

        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    case PhysicalRegionEnum::STRANDED_COIL:
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        std::string volt_var_name;
        if (PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL == static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]) ||
            PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL == static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]))
            volt_var_name = "V";
        else
            volt_var_name = "VOLT";
        if (dof_solution_system_->has_variable(volt_var_name))
        {
            auto build_fv = [&]()
            {
                const unsigned int volt_var_num = dof_solution_system_->variable_number(volt_var_name);
                vtkDataArray *amps_array = vtk_grid_->GetPointData()->GetArray("AMPS");

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                {
                    const Node &node = elem_data.elem->node_ref(dof_i);
                    const vtkIdType node_id = static_cast<vtkIdType>(node.id());
                    const Real amps_value = amps_array->GetComponent(node_id, 0);
                    const dof_id_type dof_id = node.dof_number(ApecMSystemBase::DOF_SOLUTION_SYSTEM, volt_var_num, 0);
                    if (!nodal_loads_dof_cache.insert(dof_id).second)
                        continue;
                    Fe(dof_i + volt_var_num * n_dofs_var) = amps_value;
                }
            };
            build_fv();
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
    this->build_solid97_nodal_solution();
    this->build_solid97_element_solution();
}

void ElementBuilder::build_solid97_nodal_solution()
{
    if (!this->compute_nodal_solution_vectors)
        return;

    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::NODAL_SOLUTION);
    if (api_vars.empty())
        return;

    const auto nu_matrix = this->get_mu_matrix(this->elem_data.material_num);

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     ApecMSystemBase::DOF_SOLUTION_SYSTEM);

    DenseVector<Number> An;
    for (unsigned int k = 0; k < 3; k++)
        An.append(sol_vectors[k]);

    unsigned int n_dofs_var = elem_data.elem->n_nodes();

    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    auto &nodal_vectors = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS];
    nodal_vectors.clear();
    nodal_vectors.resize(api_vars.size());
    DenseVector<Number> Bn;
    DenseVector<Number> Hn;

    bool have_b = false;
    bool have_h = false;

    auto build_b = [&]()
    {
        DenseVector<Number> Bn_x;
        AAe.vector_mult(Bn_x, An);
        Me.cholesky_solve(Bn_x, Bn);
        have_b = true;
    };

    auto build_h = [&]()
    {
        if (!have_b)
            build_b();
        Hn.resize(Bn.size());
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Hn(dof_i + k * n_dofs_var) += Bn(dof_i + k * n_dofs_var) * nu_matrix[k];
        have_h = true;
    };

    this->fill_vectors_by_api(api_vars, nodal_vectors, [&](APIVariableName api_var) -> const DenseVector<Number> *
                              {
                                  switch (api_var)
                                  {
                                  case APIVariableName::B:
                                      build_b();
                                      return have_b ? &Bn : nullptr;
                                  case APIVariableName::H:
                                      build_h();
                                      return have_h ? &Hn : nullptr;
                                  default:
                                      return nullptr;
                                  } });

    this->append_and_write_solution(*nodal_solution_system_, nodal_vectors);
}

void ElementBuilder::build_solid97_element_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::ELEMENT_SOLUTION);
    if (api_vars.empty())
        return;

    std::vector<dof_id_type> dof_indices;

    const auto nu_matrix = this->get_mu_matrix(this->elem_data.material_num);

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     ApecMSystemBase::DOF_SOLUTION_SYSTEM);

    DenseVector<Number> An;
    for (unsigned int k = 0; k < 3; k++)
        An.append(sol_vectors[k]);

    auto &element_vectors = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS];
    element_vectors.clear();
    element_vectors.resize(api_vars.size());
    DenseVector<Number> Be;
    DenseVector<Number> He;
    DenseVector<Number> Je;

    bool have_be = false;
    bool have_he = false;
    bool have_je = false;

    auto build_be = [&]()
    {
        auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];
        unsigned int n_qp = BAe.m() / 3;
        Be.resize(3);

        DenseVector<Number> Bq;
        BAe.vector_mult(Bq, An);
        for (unsigned int k = 0; k < 3; k++)
        {
            for (unsigned int qp = 0; qp < n_qp; qp++)
                Be(k) += Bq(qp + k * n_qp);
            Be(k) = Be(k) / n_qp;
        }
        have_be = true;
    };

    auto build_he = [&]()
    {
        if (!have_be)
            build_be();
        He.resize(3);
        for (unsigned int k = 0; k < 3; k++)
            He(k) = Be(k) * nu_matrix[k];
        have_he = true;
    };

    auto build_je = [&]()
    {
        auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];
        unsigned int n_qp = BAe.m() / 3;
        Je.resize(3);

        PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);
        switch (physical_region)
        {
        case PhysicalRegionEnum::NON_MAGNETIC:
        case PhysicalRegionEnum::SOFT_MAGNETIC:
        case PhysicalRegionEnum::HARD_MAGNETIC:
            break;
        case PhysicalRegionEnum::STRANDED_CONDUCTOR:
        {
            vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
            const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
            for (unsigned int k = 0; k < 3; k++)
                Je(k) = js_array->GetComponent(elem_id, k);
            break;
        }
        case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        case PhysicalRegionEnum::SOLID_CONDUCTOR:
        case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
        {
            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
            auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];

            const auto sigma_matrix = this->get_sigma_matrix(this->elem_data.material_num);

            DenseVector<Number> Vn;

            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
                Vn = sol_vectors[3];

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"), 0); // solution_derivative

            DenseVector<Number> An_dt;
            for (unsigned int k = 0; k < 3; k++)
                An_dt.append(sol_vectors_dt[k]);

            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR || physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR)
                Vn = sol_vectors_dt[3];

            DenseVector<Number> Jq_A;
            JAe.vector_mult(Jq_A, An_dt);

            DenseVector<Number> Jq_V;
            JVe.vector_mult(Jq_V, Vn);

            Jq_A += Jq_V;

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

            if (physical_region == PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL)
            {
                Vn = sol_vectors[3][0];
                En = sol_vectors[4][0];
            }

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            if (physical_region == PhysicalRegionEnum::STRANDED_COIL)
            {
                Vn = sol_vectors_dt[3][0];
                En = sol_vectors_dt[4][0];
            }

            const auto rc = this->get_real_constants(this->elem_data.real_constant_num);
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

        have_je = true;
    };

    this->fill_vectors_by_api(api_vars, element_vectors, [&](APIVariableName api_var) -> const DenseVector<Number> *
                              {
                                  switch (api_var)
                                  {
                                  case APIVariableName::Be:
                                      build_be();
                                      return have_be ? &Be : nullptr;
                                  case APIVariableName::He:
                                      build_he();
                                      return have_he ? &He : nullptr;
                                  case APIVariableName::Jt:
                                      build_je();
                                      return have_je ? &Je : nullptr;
                                  default:
                                      return nullptr;
                                  } });

    this->append_and_write_solution(*element_solution_system_, element_vectors);
}

/**
 * @file av3dnodal_element_builder.cpp
 * @brief ElementBuilder::build_av3dnodal() — element matrices/vectors for the AV3D nodal formulation,
 *        plus matching nodal-solution and element-solution post-processing methods.
 */
#include "element_builder.h"
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>
#include <libmesh/petsc_linear_solver.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <cmath>
#include <stdexcept>

#include "element_types.h"
#include "pytec3d_system.h"
#include "element_types_factory.h"
#include "variables_map.h"
#include "logging.h"

using libMesh::DenseVector;
using libMesh::dof_id_type;
using libMesh::FEBase;
using libMesh::FEInterface;
using libMesh::FEVectorBase;
using libMesh::Node;
using libMesh::Number;
using libMesh::Point;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::RealVectorValue;
using libMesh::TensorValue;

// Matrix/vector index legend: {0: Ke, 1: Ce, 2: Me, 3: AAe, 4: JAe, 5: JVe}
void ElementBuilder::build_av3dnodal()
{
    FEBase *fe = &elem_data.element_type->scalar_fe(1);
    FEVectorBase *fe_vec = &elem_data.element_type->vector_fe(0);
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
    const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_vec->get_phi();
    const std::vector<std::vector<RealGradient>> &curl_phi_vec = fe_vec->get_curl_phi();
    const std::vector<std::vector<Real>> &div_phi_vec = fe_vec->get_div_phi();

    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int n_dofs_A = fe_vec->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
    const Real nu_x = 1 / (kMu0 * this->material_value(material_id, MAT_MUR_X));
    const Real nu_y = 1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y));
    const Real nu_z = 1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z));
    TensorValue<Real> nu_tensor(
        nu_x, 0.0, 0.0,
        0.0, nu_y, 0.0,
        0.0, 0.0, nu_z);
    this->rotate_diag_tensor_for_elem(nu_tensor);
    const Real nu_penalty = std::cbrt(std::abs(nu_x * nu_y * nu_z));

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ke];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ce];
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
    auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];

    auto build_ke = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                {
                    const unsigned int dof_idx_i = 3 * (dof_i % n_dofs_var) + (dof_i / n_dofs_var);
                    const unsigned int dof_idx_j = 3 * (dof_j % n_dofs_var) + (dof_j / n_dofs_var);

                    const auto &curl_i = curl_phi_vec[dof_idx_i][qp];
                    const auto &curl_j = curl_phi_vec[dof_idx_j][qp];
                    Real curl_term = 0.0;
                    for (unsigned int k = 0; k < 3; ++k)
                        for (unsigned int l = 0; l < 3; ++l)
                            curl_term += curl_i(k) * nu_tensor(k, l) * curl_j(l);

                    const Real div_i = div_phi_vec[dof_idx_i][qp];
                    const Real div_j = div_phi_vec[dof_idx_j][qp];
                    const Real gauge_term = nu_penalty * div_i * div_j;

                    Ke(dof_i, dof_j) += w * (curl_term + gauge_term);
                }
        }
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
        {
            const Real w = JxW[qp];
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        for (unsigned int l = 0; l < 3; l++)
                            AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                                w * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
        }
    };

    auto build_me = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        for (unsigned int l = 0; l < 3; l++)
                            Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                                w * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);
        }
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

    if (!this->nonlinear_element())
        build_ke();
    build_bae();
    build_aae();
    build_me();
    build_jae();
    build_jve();

    auto physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
        break;
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        TensorValue<Real> sigma_tensor(
            1.0 / this->material_value(material_id, MAT_RHO_X), 0.0, 0.0,
            0.0, 1.0 / this->material_value(material_id, MAT_RHO_Y), 0.0,
            0.0, 0.0, 1.0 / this->material_value(material_id, MAT_RHO_Z));
        this->rotate_diag_tensor_for_elem(sigma_tensor);

        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = phi[dof_i][qp] * phi[dof_j][qp];
                        for (unsigned int k = 0; k < 3; ++k)
                            for (unsigned int l = 0; l < 3; ++l)
                                Ce(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) +=
                                    w * sigma_tensor(k, l) * nn;
                    }
            }
        };

        auto build_cav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real phi_i = phi[dof_i][qp];
                        for (unsigned int k = 0; k < 3; ++k)
                            for (unsigned int l = 0; l < 3; ++l)
                                Ce(dof_i + k * n_dofs_var, n_dofs_A + dof_j) +=
                                    w * sigma_tensor(k, l) * phi_i * dphi[dof_j][qp](l);
                    }
            }
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real phi_i = phi[dof_i][qp];
                        for (unsigned int k = 0; k < 3; ++k)
                            for (unsigned int l = 0; l < 3; ++l)
                                Ce(n_dofs_A + dof_j, dof_i + k * n_dofs_var) +=
                                    w * sigma_tensor(l, k) * phi_i * dphi[dof_j][qp](l);
                    }
            }
        };

        auto build_cvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        for (unsigned int k = 0; k < 3; ++k)
                            for (unsigned int l = 0; l < 3; ++l)
                                Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) +=
                                    w * dphi[dof_i][qp](k) * sigma_tensor(k, l) * dphi[dof_j][qp](l);
                    }
            }
        };

        auto build_kav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real phi_i = phi[dof_i][qp];
                        for (unsigned int k = 0; k < 3; ++k)
                            for (unsigned int l = 0; l < 3; ++l)
                                Ke(dof_i + k * n_dofs_var, n_dofs_A + dof_j) +=
                                    w * sigma_tensor(k, l) * phi_i * dphi[dof_j][qp](l);
                    }
            }
        };

        auto build_kvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        for (unsigned int k = 0; k < 3; ++k)
                            for (unsigned int l = 0; l < 3; ++l)
                                Ke(n_dofs_A + dof_i, n_dofs_A + dof_j) +=
                                    w * dphi[dof_i][qp](k) * sigma_tensor(k, l) * dphi[dof_j][qp](l);
                    }
            }
        };

        build_caa();
        if (elem_data.element_type->get_OPTIONs()[4] == 0 &&
            physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR)
        {
            build_kav();
            build_cva();
            build_kvv();
        }
        else
        {
            build_cav();
            build_cva();
            build_cvv();
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        const Real sc = this->constant_value(constant_id, CONST_CS);
        const Real nc = this->constant_value(constant_id, CONST_TURNS);
        const Real vc = this->constant_value(constant_id, CONST_VOL);
        RealVectorValue t_global(
            this->constant_value(constant_id, CONST_T_X),
            this->constant_value(constant_id, CONST_T_Y),
            this->constant_value(constant_id, CONST_T_Z));
        this->rotate_vector_for_elem(t_global);
        const Real tx = t_global(0);
        const Real ty = t_global(1);
        const Real tz = t_global(2);
        const Real r_coil = this->constant_value(constant_id, CONST_R);
        const Real sym = this->constant_value(constant_id, CONST_S);

        auto build_coil_coupling_timeintegrated = [&]()
        {
            const Real kav_coeff = -(nc / (sc * r_coil));
            const Real kvv_coeff = (1.0 / (r_coil * vc));
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        const Real kav_x = kav_coeff * tx * nn;
                        const Real kav_y = kav_coeff * ty * nn;
                        const Real kav_z = kav_coeff * tz * nn;
                        const Real kvv = kvv_coeff * nn;

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

        auto build_coil_coupling_notimeintegrated = [&]()
        {
            const Real kav_coeff = -(nc / (sc * r_coil));
            const Real kvv_coeff = (1.0 / (r_coil * vc));
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        const Real kav_x = kav_coeff * tx * nn;
                        const Real kav_y = kav_coeff * ty * nn;
                        const Real kav_z = kav_coeff * tz * nn;
                        const Real kvv = kvv_coeff * nn;

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

        if (elem_data.element_type->get_OPTIONs()[4] == 0)
            build_coil_coupling_notimeintegrated();
        else
            build_coil_coupling_timeintegrated();
        break;
    }
    default:
        pytec3d::logging::error_and_abort("Physical region not implemented");
    }
}

void ElementBuilder::build_av3dnodal_rhs()
{
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};
    this->build_vectors(&vector_sizes, VectorsSet::RHS_VECTORS);
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    FEVectorBase *fe_vec = &elem_data.element_type->vector_fe(0);
    FEBase *fe = &elem_data.element_type->scalar_fe(1);
    unsigned int n_dofs_A = FEInterface::n_shape_functions(fe_vec->get_fe_type(),
                                                           elem_data.elem);
    unsigned int n_dofs_var = FEInterface::n_shape_functions(fe->get_fe_type(),
                                                             elem_data.elem);

    const std::vector<Real> &JxW = fe_vec->get_JxW();
    const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_vec->get_phi();
    const std::vector<std::vector<RealGradient>> &curl_phi_vec = fe_vec->get_curl_phi();

    unsigned int nq;

    RealVectorValue mg(0.0, 0.0, 0.0);

    auto physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        fe_vec->reinit(elem_data.elem);
        nq = fe_vec->n_quadrature_points();
    }
    default:
    {
        const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);

        RealVectorValue mg(
            this->material_value(material_id, MAT_M_X),
            this->material_value(material_id, MAT_M_Y),
            this->material_value(material_id, MAT_M_Z));
        this->rotate_vector_for_elem(mg);
        if (mg(0) != 0.0 || mg(1) != 0.0 || mg(2) != 0.0)
        {
            fe_vec->reinit(elem_data.elem);
            nq = fe_vec->n_quadrature_points();
        }
    }
    }

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
        RealVectorValue js(
            js_array->GetComponent(elem_id, 0),
            js_array->GetComponent(elem_id, 1),
            js_array->GetComponent(elem_id, 2));
        this->rotate_vector_for_elem(js);

        auto build_fa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                        Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi_vec[dof_i * 3 + k][qp](k) * js(k);
        };
        build_fa();

        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        std::string volt_var_name = "V";
        if (physical_region == PhysicalRegionEnum::STRANDED_COIL)
            volt_var_name = "VOLT";
        if (dof_solution_system_->has_variable(volt_var_name))
        {
            auto build_fv = [&]()
            {
                const unsigned int volt_var_num = dof_solution_system_->variable_number(volt_var_name);
                vtkDataArray *amps_curr_array = nullptr;
                if (volt_var_name == "V")
                    amps_curr_array = vtk_grid_->GetPointData()->GetArray("AMPS");
                else
                    amps_curr_array = vtk_grid_->GetPointData()->GetArray("CURR");

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                {
                    const Node &node = elem_data.elem->node_ref(dof_i);
                    const vtkIdType node_id = static_cast<vtkIdType>(node.id());
                    const Real amps_value = amps_curr_array->GetComponent(node_id, 0);
                    const dof_id_type dof_id = node.dof_number(dof_solution_system_->number(), volt_var_num, 0);
                    if (!this->should_inject_nodal_load_dof(volt_var_name, dof_id))
                        continue;
                    Fe(dof_i + n_dofs_A) = amps_value;
                }
            };
            build_fv();
        }
        break;
    }
    default:
        pytec3d::logging::error_and_abort("Physical region not implemented");
    }

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
        break;
    default:
    {
        if (mg(0) != 0.0 || mg(1) != 0.0 || mg(2) != 0.0)
        {
            auto build_mg = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                    for (unsigned int k = 0; k < 3; k++)
                        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                        {
                            const auto &curl_w = curl_phi_vec[dof_i * 3 + k][qp];
                            const Real curl_dot_m =
                                curl_w(0) * mg(0) +
                                curl_w(1) * mg(1) +
                                curl_w(2) * mg(2);
                            Fe(dof_i + k * n_dofs_var) += JxW[qp] * curl_dot_m;
                        }
            };
            build_mg();
        }
    }
    }

    return;
}

void ElementBuilder::build_av3dnodal_nodal_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, pyTEC3DSystemBase::NODAL_SOLUTION);
    if (api_vars.empty())
        return;

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    TensorValue<Real> nu_tensor(
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)), 0.0, 0.0,
        0.0, 1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)), 0.0,
        0.0, 0.0, 1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z)));
    this->rotate_diag_tensor_for_elem(nu_tensor);

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);

    DenseVector<Number> An;
    for (unsigned int k = 0; k < 3; k++)
        An.append(sol_vectors[k]);

    FEBase *fe = &elem_data.element_type->scalar_fe(1);
    unsigned int n_dofs_var = FEInterface::n_shape_functions(fe->get_fe_type(),
                                                             elem_data.elem);

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
        if (this->nonlinear_element())
        {
            const auto &material_bh = bh_curves_->find(material_id)->second;
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; ++dof_i)
            {
                Real B[3] = {
                    static_cast<Real>(Bn(dof_i + 0 * n_dofs_var)),
                    static_cast<Real>(Bn(dof_i + 1 * n_dofs_var)),
                    static_cast<Real>(Bn(dof_i + 2 * n_dofs_var))};
                Real H[3];
                Real jac[3][3];
                this->evaluate_isotropic_bh(material_bh, B, 3, H, jac);
                Hn(dof_i + 0 * n_dofs_var) = H[0];
                Hn(dof_i + 1 * n_dofs_var) = H[1];
                Hn(dof_i + 2 * n_dofs_var) = H[2];
            }
        }
        else
        {
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; ++dof_i)
                for (unsigned int k = 0; k < 3; ++k)
                    for (unsigned int l = 0; l < 3; ++l)
                        Hn(dof_i + k * n_dofs_var) += nu_tensor(k, l) * Bn(dof_i + l * n_dofs_var);
        }
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

void ElementBuilder::build_av3dnodal_element_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, pyTEC3DSystemBase::ELEMENT_SOLUTION);
    if (api_vars.empty())
        return;

    std::vector<dof_id_type> dof_indices;

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
    TensorValue<Real> nu_tensor(
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)), 0.0, 0.0,
        0.0, 1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)), 0.0,
        0.0, 0.0, 1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z)));
    this->rotate_diag_tensor_for_elem(nu_tensor);

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);

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
        if (this->nonlinear_element())
        {
            const auto &material_bh = bh_curves_->find(material_id)->second;
            Real B[3] = {
                static_cast<Real>(Be(0)),
                static_cast<Real>(Be(1)),
                static_cast<Real>(Be(2))};
            Real H[3];
            Real jac[3][3];
            this->evaluate_isotropic_bh(material_bh, B, 3, H, jac);
            He(0) = H[0];
            He(1) = H[1];
            He(2) = H[2];
        }
        else
        {
            for (unsigned int k = 0; k < 3; ++k)
                for (unsigned int l = 0; l < 3; ++l)
                    He(k) += nu_tensor(k, l) * Be(l);
        }
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
            RealVectorValue js(
                js_array->GetComponent(elem_id, 0),
                js_array->GetComponent(elem_id, 1),
                js_array->GetComponent(elem_id, 2));
            this->rotate_vector_for_elem(js);
            for (unsigned int k = 0; k < 3; k++)
                Je(k) = js(k);
            break;
        }
        case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        case PhysicalRegionEnum::SOLID_CONDUCTOR:
        {
            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
            auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];

            TensorValue<Real> sigma_tensor(
                1.0 / this->material_value(material_id, MAT_RHO_X), 0.0, 0.0,
                0.0, 1.0 / this->material_value(material_id, MAT_RHO_Y), 0.0,
                0.0, 0.0, 1.0 / this->material_value(material_id, MAT_RHO_Z));
            this->rotate_diag_tensor_for_elem(sigma_tensor);

            DenseVector<Number> Vn;
            const unsigned int v_var_num = dof_solution_system_->variable_number("V");

            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR &&
                elem_data.element_type->get_OPTIONs()[4] == 0)
                Vn = sol_vectors[v_var_num];

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            DenseVector<Number> An_dt;
            for (unsigned int k = 0; k < 3; k++)
                An_dt.append(sol_vectors_dt[k]);

            if (physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR ||
                elem_data.element_type->get_OPTIONs()[4] != 0)
                Vn = sol_vectors_dt[v_var_num];

            DenseVector<Number> Jq_A;
            JAe.vector_mult(Jq_A, An_dt);

            DenseVector<Number> Jq_V;
            JVe.vector_mult(Jq_V, Vn);

            Jq_A += Jq_V;

            RealVectorValue j_avg(0.0, 0.0, 0.0);
            for (unsigned int l = 0; l < 3; ++l)
            {
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                    j_avg(l) -= Jq_A(qp + l * n_qp);
                j_avg(l) /= n_qp;
            }
            for (unsigned int k = 0; k < 3; ++k)
                for (unsigned int l = 0; l < 3; ++l)
                    Je(k) += sigma_tensor(k, l) * j_avg(l);
            break;
        }
        case PhysicalRegionEnum::STRANDED_COIL:
        {
            Real Vn = 0.0;
            Real En = 0.0;
            const unsigned int volt_var_num = dof_solution_system_->variable_number("VOLT");
            const unsigned int emf_var_num = dof_solution_system_->variable_number("EMF");

            if (elem_data.element_type->get_OPTIONs()[4] == 0)
            {
                Vn = sol_vectors[volt_var_num][0];
                En = sol_vectors[emf_var_num][0];
            }

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            if (elem_data.element_type->get_OPTIONs()[4] != 0)
            {
                Vn = sol_vectors_dt[volt_var_num][0];
                En = sol_vectors_dt[emf_var_num][0];
            }

            const Real sc = this->constant_value(constant_id, CONST_CS);
            const Real nc = this->constant_value(constant_id, CONST_TURNS);
            RealVectorValue t_global(
                this->constant_value(constant_id, CONST_T_X),
                this->constant_value(constant_id, CONST_T_Y),
                this->constant_value(constant_id, CONST_T_Z));
            this->rotate_vector_for_elem(t_global);
            const Real r_coil = this->constant_value(constant_id, CONST_R);

            const Real coeff = nc / (sc * r_coil) * (Vn - En);

            for (unsigned int k = 0; k < 3; k++)
                Je(k) = coeff * t_global(k);
            break;
        }
        default:
            pytec3d::logging::error_and_abort("Physical region not implemented");
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

void ElementBuilder::build_av3dnodal_tangent_stiffness_residual()
{
    FEBase *fe = &elem_data.element_type->scalar_fe(1);
    FEVectorBase *fe_vec = &elem_data.element_type->vector_fe(0);
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<RealVectorValue>> &curl_phi_vec = fe_vec->get_curl_phi();
    const std::vector<std::vector<Real>> &div_phi_vec = fe_vec->get_div_phi();

    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs}};
    this->build_matrices(&matrix_sizes, MatricesSet::TANGENT_STIFFNESS_MATRICES);
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};
    this->build_vectors(&vector_sizes, VectorsSet::RESIDUAL_VECTORS);

    const unsigned int n_dofs_var = fe->n_shape_functions();
    const unsigned int n_dofs_A = fe_vec->n_shape_functions();
    const unsigned int nq = fe->n_quadrature_points();

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    if (!this->nonlinear_element())
        return;

    const auto &material_bh = bh_curves_->find(material_id)->second;

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);
    DenseVector<Number> An;
    for (unsigned int k = 0; k < 3; ++k)
        An.append(sol_vectors[k]);

    auto &Ke = this->matrices[MatricesSet::TANGENT_STIFFNESS_MATRICES][0];
    auto &Re = this->vectors[VectorsSet::RESIDUAL_VECTORS][0];

    libMesh::DenseMatrix<libMesh::Number> BAe_qp;
    BAe_qp.resize(3 * nq, n_dofs_A);
    for (unsigned int qp = 0; qp < nq; ++qp)
        for (unsigned int k = 0; k < 3; ++k)
            for (unsigned int dof_j = 0; dof_j < n_dofs_A; ++dof_j)
            {
                const unsigned int dof_idx_j = 3 * (dof_j % n_dofs_var) + (dof_j / n_dofs_var);
                BAe_qp(qp + k * nq, dof_j) += curl_phi_vec[dof_idx_j][qp](k);
            }

    DenseVector<Number> Bq;
    BAe_qp.vector_mult(Bq, An);

    libMesh::DenseMatrix<libMesh::Number> DAe_qp;
    DAe_qp.resize(nq, n_dofs_A);
    for (unsigned int qp = 0; qp < nq; ++qp)
        for (unsigned int dof_j = 0; dof_j < n_dofs_A; ++dof_j)
        {
            const unsigned int dof_idx_j = 3 * (dof_j % n_dofs_var) + (dof_j / n_dofs_var);
            DAe_qp(qp, dof_j) += div_phi_vec[dof_idx_j][qp];
        }

    DenseVector<Number> divA_qp;
    DAe_qp.vector_mult(divA_qp, An);

    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        Real B[3] = {
            static_cast<Real>(Bq(qp + 0 * nq)),
            static_cast<Real>(Bq(qp + 1 * nq)),
            static_cast<Real>(Bq(qp + 2 * nq))};
        Real H[3];
        Real jac[3][3];
        Real nu_penalty;
        this->evaluate_isotropic_bh(material_bh, B, 3, H, jac, &nu_penalty);

        const Real div_a_qp = static_cast<Real>(divA_qp(qp));

        for (unsigned int dof_i = 0; dof_i < n_dofs_A; ++dof_i)
        {
            const unsigned int dof_idx_i = 3 * (dof_i % n_dofs_var) + (dof_i / n_dofs_var);
            Real residual_term = 0.0;
            for (unsigned int a = 0; a < 3; ++a)
                residual_term += curl_phi_vec[dof_idx_i][qp](a) * H[a];
            residual_term += nu_penalty * div_phi_vec[dof_idx_i][qp] * div_a_qp;
            Re(dof_i) += JxW[qp] * residual_term;

            for (unsigned int dof_j = 0; dof_j < n_dofs_A; ++dof_j)
            {
                const unsigned int dof_idx_j = 3 * (dof_j % n_dofs_var) + (dof_j / n_dofs_var);
                Real tangential_term = 0.0;
                for (unsigned int a = 0; a < 3; ++a)
                    for (unsigned int b = 0; b < 3; ++b)
                        tangential_term += curl_phi_vec[dof_idx_i][qp](a) * jac[a][b] * curl_phi_vec[dof_idx_j][qp](b);
                const Real gauge_term = nu_penalty * div_phi_vec[dof_idx_i][qp] * div_phi_vec[dof_idx_j][qp];
                Ke(dof_i, dof_j) += JxW[qp] * (tangential_term + gauge_term);
            }
        }
    }
}

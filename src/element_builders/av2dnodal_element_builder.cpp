/**
 * @file av2dnodal_element_builder.cpp
 * @brief ElementBuilder::build_av2dnodal() — element matrices/vectors for the AV2D nodal
 *        formulation, plus matching nodal-solution and element-solution post-processing methods.
 */
#include "element_builder.h"
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

#include <stdexcept>

#include "element_types.h"
#include "apec_m_system.h"
#include "coordinate_system.h"
#include "element_types_factory.h"
#include "variables_map.h"
#include "apec_m_logging.h"

using libMesh::DenseVector;
using libMesh::dof_id_type;
using libMesh::Elem;
using libMesh::FEBase;
using libMesh::Node;
using libMesh::Number;
using libMesh::Point;
using libMesh::Real;
using libMesh::RealGradient;

void ElementBuilder::build_av2dnodal()
{
    FEBase *fe = &elem_data.element_type->scalar_fe(0);
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
    const std::vector<Point> &xyz = fe->get_xyz();

    fe->reinit(elem_data.elem);
    this->prepare_element_matrix_building();

    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
    const Real nu_local[3] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z))};
    Real nu_tensor[3][3] = {{0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0}};
    this->rotate_diag_tensor_for_elem(nu_local, nu_tensor);

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ke];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ce];
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
    auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];

    unsigned int axisymm = elem_data.element_type->get_OPTIONs()[3];

    const Real r_sign = axisymm ? -1.0 : 1.0;
    const Real two_pi = 2.0 * libMesh::pi;

    Real thk = this->constant_value(constant_id, CONST_THK);

    std::vector<Real> r_qp(nq, 0.0);
    std::vector<Real> r_mult_qp(nq, 1.0);
    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        const Real r = std::sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        r_qp[qp] = r;
        r_mult_qp[qp] = axisymm ? (two_pi * r) : 1.0;
    }

    auto build_ke = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            const Real r = r_qp[qp];
            const Real r_mult = r_mult_qp[qp];

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                {
                    const Real bx_i = r_sign * dphi[dof_i][qp](1);
                    const Real by_i = -r_sign * (dphi[dof_i][qp](0) + axisymm * phi[dof_i][qp] / r);
                    const Real bx_j = r_sign * dphi[dof_j][qp](1);
                    const Real by_j = -r_sign * (dphi[dof_j][qp](0) + axisymm * phi[dof_j][qp] / r);
                    const Real hx_j = nu_tensor[0][0] * bx_j + nu_tensor[0][1] * by_j;
                    const Real hy_j = nu_tensor[1][0] * bx_j + nu_tensor[1][1] * by_j;
                    Ke(dof_i, dof_j) += w * (bx_i * hx_j + by_i * hy_j) * r_mult * thk;
                }
        }
    };

    auto build_bae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            // Here maybe is more suited the centroid radius
            const Real r = r_qp[qp];

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                BAe(qp, dof_i) += r_sign * dphi[dof_i][qp](1);
                BAe(qp + nq, dof_i) -= r_sign * (dphi[dof_i][qp](0) + axisymm * phi[dof_i][qp] / r);
            }
        }
    };

    auto build_jae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                JAe(qp, dof_i) += phi[dof_i][qp];
    };

    auto build_aae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            const Real r_mult = r_mult_qp[qp];

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    AAe(dof_i, dof_j) += w * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](1)) * r_mult * thk;
        }

        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            const Real r = r_qp[qp];
            const Real r_mult = r_mult_qp[qp];

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    AAe(dof_i + n_dofs_var, dof_j) -=
                        w * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](0) - axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r) * r_mult * thk;
        }
    };

    auto build_me = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            const Real r_mult = r_mult_qp[qp];

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Me(dof_i, dof_j) += w * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
        }

        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            const Real r_mult = r_mult_qp[qp];

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Me(dof_i + n_dofs_var, dof_j + n_dofs_var) += w * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
        }
    };

    if (!this->nonlinear_element())
        build_ke();
    build_bae();
    build_jae();
    build_aae();
    build_me();

    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

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
        const Real sigma_local[3] = {
            1.0 / this->material_value(material_id, MAT_RHO_X),
            1.0 / this->material_value(material_id, MAT_RHO_Y),
            1.0 / this->material_value(material_id, MAT_RHO_Z)};
        Real sigma_tensor[3][3] = {{0.0, 0.0, 0.0},
                                   {0.0, 0.0, 0.0},
                                   {0.0, 0.0, 0.0}};
        this->rotate_diag_tensor_for_elem(sigma_local, sigma_tensor);
        const Real sigma = sigma_tensor[2][2];

        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                const Real r_mult = r_mult_qp[qp];

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        Ce(dof_i, dof_j) += w * sigma * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
            }
        };

        build_caa();

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR)
        {
            auto build_coupling_timeintegrated = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                {
                    const Real w = JxW[qp];
                    const Real r_mult = r_mult_qp[qp];
                    const Real inv_rmult_thk = 1.0 / (r_mult * thk);

                    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                        for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        {
                            const Real w_sigma_nn = w * sigma * (phi[dof_i][qp] * phi[dof_j][qp]);
                            Ce(dof_i, dof_j + n_dofs_var) -=
                                w_sigma_nn;
                            Ce(dof_i + n_dofs_var, dof_j) -=
                                w_sigma_nn;
                            Ce(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                                w_sigma_nn * inv_rmult_thk;
                        }
                }
            };
            auto build_coupling_notimeintegrated = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                {
                    const Real w = JxW[qp];
                    const Real r_mult = r_mult_qp[qp];
                    const Real inv_rmult_thk = 1.0 / (r_mult * thk);

                    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                        for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        {
                            const Real w_sigma_nn = w * sigma * (phi[dof_i][qp] * phi[dof_j][qp]);
                            Ke(dof_i, dof_j + n_dofs_var) -=
                                w_sigma_nn;
                            Ce(dof_i + n_dofs_var, dof_j) -=
                                w_sigma_nn;
                            Ke(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                                w_sigma_nn * inv_rmult_thk;
                        }
                }
            };
            if (elem_data.element_type->get_OPTIONs()[4] == 0)
                build_coupling_notimeintegrated();
            else
                build_coupling_timeintegrated();
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        const Real sc = this->constant_value(constant_id, CONST_CS);
        const Real nc = this->constant_value(constant_id, CONST_TURNS);
        const Real rad = this->constant_value(constant_id, CONST_RAD);
        Real rad_factor = 1.0;
        if (axisymm)
            rad_factor = 2 * libMesh::pi * rad;
        const Real vc = rad_factor * sc * thk;
        const Real t_local[3] = {
            this->constant_value(constant_id, CONST_T_X),
            this->constant_value(constant_id, CONST_T_Y),
            this->constant_value(constant_id, CONST_T_Z)};
        Real t_global[3] = {0.0, 0.0, 0.0};
        this->rotate_vector_for_elem(t_local, t_global);
        const Real tz = t_global[2];
        const Real r_coil = this->constant_value(constant_id, CONST_R);
        const Real sym = this->constant_value(constant_id, CONST_S);

        const Real kav_scale = -(nc / (sc * r_coil)) * tz;
        const Real kvv_scale = 1.0 / (r_coil * vc);

        auto build_coil_coupling_timeintegrated = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                const Real r_mult = r_mult_qp[qp];

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = w * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
                        const Real kav = kav_scale * nn;
                        const Real kvv = kvv_scale * nn;

                        // Time-integrated A-VOLT-EMF formulation (OPTION(2)=0/1).
                        Ce(dof_i, dof_j + n_dofs_var) += kav;
                        Ce(dof_i, dof_j + 2 * n_dofs_var) -= kav;
                        Ce(dof_i + n_dofs_var, dof_j + n_dofs_var) += sym * kvv;
                        // Ke(dof_i + n_dofs_var, dof_j + 2 * n_dofs_var) -= kvv;
                        // Ke(dof_i + 2 * n_dofs_var, dof_j + n_dofs_var) -= kvv;
                        Ce(dof_i + 2 * n_dofs_var, dof_j + 2 * n_dofs_var) -= kvv;

                        Ce(dof_i + n_dofs_var, dof_j) += sym * kav;
                        Ce(dof_i + 2 * n_dofs_var, dof_j) -= kav;
                    }
            }
        };
        auto build_coil_coupling_notimeintegrated = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                const Real r_mult = r_mult_qp[qp];

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = w * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
                        const Real kav = kav_scale * nn;
                        const Real kvv = kvv_scale * nn;

                        // No-time-integrated A-VOLT-EMF formulation (OPTION(2)=0/1).
                        Ke(dof_i, dof_j + n_dofs_var) += kav;
                        Ke(dof_i, dof_j + 2 * n_dofs_var) -= kav;
                        Ke(dof_i + n_dofs_var, dof_j + n_dofs_var) += sym * kvv;
                        // Ke(dof_i + n_dofs_var, dof_j + 2 * n_dofs_var) -= kvv;
                        // Ke(dof_i + 2 * n_dofs_var, dof_j + n_dofs_var) -= kvv;
                        Ke(dof_i + 2 * n_dofs_var, dof_j + 2 * n_dofs_var) -= kvv;

                        Ce(dof_i + n_dofs_var, dof_j) += sym * kav;
                        Ce(dof_i + 2 * n_dofs_var, dof_j) -= kav;
                    }
            }
        };
        if (elem_data.element_type->get_OPTIONs()[4] == 0)
            build_coil_coupling_notimeintegrated();
        else
            build_coil_coupling_timeintegrated();
        break;
    }
    default:
        apec::logging::error_and_abort("Physical region not implemented");
    }
}

void ElementBuilder::build_av2dnodal_rhs()
{
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};
    this->build_vectors(&vector_sizes, VectorsSet::RHS_VECTORS);
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    FEBase *fe = &elem_data.element_type->scalar_fe(0);
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
    const std::vector<Point> &xyz = fe->get_xyz();

    unsigned int n_dofs_var = elem_data.elem->n_nodes();
    unsigned int nq;

    Real mg[2] = {0.0, 0.0};
    Real thk;

    unsigned int axisymm = elem_data.element_type->get_OPTIONs()[3];

    auto physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        fe->reinit(elem_data.elem);
        nq = fe->n_quadrature_points();
        const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
        thk = this->constant_value(constant_id, CONST_THK);
    }
    default:
    {
        const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);

        const Real mg_local[3] = {
            this->material_value(material_id, MAT_M_X),
            this->material_value(material_id, MAT_M_Y),
            this->material_value(material_id, MAT_M_Z)};
        Real mg_global[3] = {0.0, 0.0, 0.0};
        this->rotate_vector_for_elem(mg_local, mg_global);
        mg[0] = mg_global[0];
        mg[1] = mg_global[1];

        if (mg[0] != 0.0 || mg[1] != 0.0)
        {
            fe->reinit(elem_data.elem);
            nq = fe->n_quadrature_points();
            const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
            thk = this->constant_value(constant_id, CONST_THK);
        }
    }
    }

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
        const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
        const Real js_local[3] = {
            js_array->GetComponent(elem_id, 0),
            js_array->GetComponent(elem_id, 1),
            js_array->GetComponent(elem_id, 2)};
        Real js_global[3] = {0.0, 0.0, 0.0};
        this->rotate_vector_for_elem(js_local, js_global);
        const Real js_z = js_global[2];

        auto build_fa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    Fe(dof_i) += JxW[qp] * phi[dof_i][qp] * r_mult * thk * js_z;
            }
        };
        build_fa();
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        if (dof_solution_system_->has_variable("VOLT") && dof_solution_system_->has_variable("A_z"))
        {
            auto build_fv = [&]()
            {
                const unsigned int volt_var_num = dof_solution_system_->variable_number("VOLT");
                vtkDataArray *amps_array = vtk_grid_->GetPointData()->GetArray("AMPS");

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                {
                    const Node &node = elem_data.elem->node_ref(dof_i);
                    const vtkIdType node_id = static_cast<vtkIdType>(node.id());
                    const Real amps_value = amps_array->GetComponent(node_id, 0);
                    const dof_id_type dof_id = node.dof_number(ApecMSystemBase::DOF_SOLUTION_SYSTEM, volt_var_num, 0);
                    if (!nodal_loads_dof_cache.insert(dof_id).second)
                        continue;
                    Fe(dof_i + n_dofs_var) = amps_value;
                }
            };
            build_fv();
        }
        break;
    }
    default:
        apec::logging::error_and_abort("Physical region not implemented");
    }

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
        break;
    default:
    {
        if (mg[0] != 0.0 || mg[1] != 0.0)
        {
            // Same curl operator (bx_i, by_i) as build_ke()/build_bae(), so that the magnetization
            // source term stays consistent with B = curl(A) under the axisymmetric r_sign/curvature terms.
            const Real r_sign = axisymm ? -1.0 : 1.0;
            auto build_mg = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                {
                    Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                    Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    {
                        const Real bx_i = r_sign * dphi[dof_i][qp](1);
                        const Real by_i = -r_sign * (dphi[dof_i][qp](0) + axisymm * phi[dof_i][qp] / r);
                        Fe(dof_i) += JxW[qp] * (bx_i * mg[0] + by_i * mg[1]) * thk * r_mult;
                    }
                }
            };
            build_mg();
        }
    }
    }
    return;
}

void ElementBuilder::build_av2dnodal_nodal_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::NODAL_SOLUTION);
    if (api_vars.empty())
        return;

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const Real nu_local[3] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z))};
    Real nu_tensor[3][3] = {{0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0}};
    this->rotate_diag_tensor_for_elem(nu_local, nu_tensor);

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    const auto &An = sol_vectors[0];

    auto &nodal_vectors = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS];
    nodal_vectors.clear();
    nodal_vectors.resize(api_vars.size());
    DenseVector<Number> Bn;
    DenseVector<Number> Hn;
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    unsigned int n_dofs_var = elem_data.elem->n_nodes();

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
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                const Real bx = Bn(dof_i);
                const Real by = Bn(dof_i + n_dofs_var);
                Real B[3] = {bx, by, 0.0};
                Real H[3];
                Real jac[3][3];
                this->evaluate_isotropic_bh(material_bh, B, 2, H, jac);
                Hn(dof_i) = H[0];
                Hn(dof_i + n_dofs_var) = H[1];
            }
        }
        else
        {
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                const Real bx = Bn(dof_i);
                const Real by = Bn(dof_i + n_dofs_var);
                Hn(dof_i) = nu_tensor[0][0] * bx + nu_tensor[0][1] * by;
                Hn(dof_i + n_dofs_var) = nu_tensor[1][0] * bx + nu_tensor[1][1] * by;
            }
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

void ElementBuilder::build_av2dnodal_element_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::ELEMENT_SOLUTION);
    if (api_vars.empty())
        return;

    std::vector<dof_id_type> dof_indices;

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
    const Real nu_local[3] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z))};
    Real nu_tensor[3][3] = {{0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0}};
    this->rotate_diag_tensor_for_elem(nu_local, nu_tensor);

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    const auto &An = sol_vectors[0];

    bool axisymm = this->elem_data.element_type->get_OPTIONs()[3];
    Point centroid = this->elem_data.elem->vertex_average();

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
        unsigned int n_qp = BAe.m() / 2;
        Be.resize(2);

        DenseVector<Number> Bq(2 * n_qp);
        BAe.vector_mult(Bq, An);
        for (unsigned int k = 0; k < 2; k++)
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
        He.resize(2);
        if (this->nonlinear_element())
        {
            const auto &material_bh = bh_curves_->find(material_id)->second;
            Real B[3] = {
                static_cast<Real>(Be(0)),
                static_cast<Real>(Be(1)),
                0.0};
            Real H[3];
            Real jac[3][3];
            this->evaluate_isotropic_bh(material_bh, B, 2, H, jac);
            He(0) = H[0];
            He(1) = H[1];
        }
        else
        {
            He(0) = nu_tensor[0][0] * Be(0) + nu_tensor[0][1] * Be(1);
            He(1) = nu_tensor[1][0] * Be(0) + nu_tensor[1][1] * Be(1);
        }
        have_he = true;
    };

    auto build_je = [&]()
    {
        auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];
        unsigned int n_qp = BAe.m() / 2;
        Je.resize(1);

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
            const Real js_local[3] = {
                js_array->GetComponent(elem_id, 0),
                js_array->GetComponent(elem_id, 1),
                js_array->GetComponent(elem_id, 2)};
            Real js_global[3] = {0.0, 0.0, 0.0};
            this->rotate_vector_for_elem(js_local, js_global);
            Je(0) = js_global[2];
            break;
        }
        case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        case PhysicalRegionEnum::SOLID_CONDUCTOR:
        {
            Real Vn = 0.0;
            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR &&
                elem_data.element_type->get_OPTIONs()[4] == 0)
                Vn = sol_vectors[1][0];

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR &&
                elem_data.element_type->get_OPTIONs()[4] != 0)
                Vn = sol_vectors_dt[1][0];

            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];

            const Real sigma_local[3] = {
                1.0 / this->material_value(material_id, MAT_RHO_X),
                1.0 / this->material_value(material_id, MAT_RHO_Y),
                1.0 / this->material_value(material_id, MAT_RHO_Z)};
            Real sigma_tensor[3][3] = {{0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0}};
            this->rotate_diag_tensor_for_elem(sigma_local, sigma_tensor);
            const Real sigma = sigma_tensor[2][2];

            const auto &An_dt = sol_vectors_dt[0];

            DenseVector<Number> Jq_A(n_qp);
            JAe.vector_mult(Jq_A, An_dt);

            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(0) -= Jq_A(qp);
            Je(0) = Je(0) * sigma / n_qp;

            if (physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR)
                break;

            Real thk = this->constant_value(constant_id, CONST_THK);

            Real V_factor = sigma * Vn / thk;
            if (axisymm)
            {
                FEBase *fe = &elem_data.element_type->scalar_fe(1);
                const std::vector<Point> &xyz = fe->get_xyz();
                fe->reinit(elem_data.elem);
                Real Jv = 0.0;
                for (unsigned int qp = 0; qp < JAe.m(); qp++)
                    for (unsigned int dof_i = 0; dof_i < JAe.n(); dof_i++)
                        Jv += JAe(qp, dof_i) / (2 * libMesh::pi * xyz[qp](0));
                V_factor *= Jv / n_qp;
            }
            Je(0) += V_factor;
            break;
        }
        case PhysicalRegionEnum::STRANDED_COIL:
        {
            Real Vn = 0.0;
            Real En = 0.0;

            if (elem_data.element_type->get_OPTIONs()[4] == 0)
            {
                Vn = sol_vectors[1][0];
                En = sol_vectors[2][0];
            }

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            if (elem_data.element_type->get_OPTIONs()[4] != 0)
            {
                Vn = sol_vectors_dt[1][0];
                En = sol_vectors_dt[2][0];
            }

            const Real sc = this->constant_value(constant_id, CONST_CS);
            const Real nc = this->constant_value(constant_id, CONST_TURNS);
            const Real t_local[3] = {
                this->constant_value(constant_id, CONST_T_X),
                this->constant_value(constant_id, CONST_T_Y),
                this->constant_value(constant_id, CONST_T_Z)};
            Real t_global[3] = {0.0, 0.0, 0.0};
            this->rotate_vector_for_elem(t_local, t_global);
            const Real tz = t_global[2];
            const Real r_coil = this->constant_value(constant_id, CONST_R);

            Je(0) = (nc / (sc * r_coil)) * tz * (Vn - En);
            break;
        }
        default:
            apec::logging::error_and_abort("Physical region not implemented");
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

void ElementBuilder::build_av2dnodal_tangent_stiffness_residual()
{
    FEBase *fe = &elem_data.element_type->scalar_fe(0);
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
    const std::vector<Point> &xyz = fe->get_xyz();

    fe->reinit(elem_data.elem);

    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs}};
    this->build_matrices(&matrix_sizes, MatricesSet::TANGENT_STIFFNESS_MATRICES);
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};
    this->build_vectors(&vector_sizes, VectorsSet::RESIDUAL_VECTORS);

    const unsigned int n_dofs_var = fe->n_shape_functions();
    const unsigned int nq = fe->n_quadrature_points();

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    auto &Ke = this->matrices[MatricesSet::TANGENT_STIFFNESS_MATRICES][0];
    auto &Re = this->vectors[VectorsSet::RESIDUAL_VECTORS][0];

    if (!this->nonlinear_element())
        return;
    const auto &material_bh = bh_curves_->find(material_id)->second;

    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    const auto &Az = sol_vectors[0];

    const unsigned int axisymm = elem_data.element_type->get_OPTIONs()[3];
    const Real r_sign = axisymm ? -1.0 : 1.0;
    const Real two_pi = 2.0 * libMesh::pi;
    const Real thk = this->constant_value(static_cast<dof_id_type>(this->elem_data.real_constant_num), CONST_THK);

    std::vector<Real> r_qp(nq, 0.0);
    std::vector<Real> r_mult_qp(nq, 1.0);
    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        const Real r = std::sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        r_qp[qp] = r;
        r_mult_qp[qp] = axisymm ? (two_pi * r) : 1.0;
    }

    libMesh::DenseMatrix<libMesh::Number> BAe_qp;
    BAe_qp.resize(2 * nq, n_dofs_var);
    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        const Real r = r_qp[qp];
        for (unsigned int dof_j = 0; dof_j < n_dofs_var; ++dof_j)
        {
            BAe_qp(qp, dof_j) += r_sign * dphi[dof_j][qp](1);
            BAe_qp(qp + nq, dof_j) -= r_sign * (dphi[dof_j][qp](0) + axisymm * phi[dof_j][qp] / r);
        }
    }

    DenseVector<Number> Bq;
    BAe_qp.vector_mult(Bq, Az);
    for (unsigned int qp = 0; qp < nq; ++qp)
    {
        Real B[3] = {static_cast<Real>(Bq(qp)), static_cast<Real>(Bq(qp + nq)), 0.0};
        Real H[3];
        Real jac[3][3];
        this->evaluate_isotropic_bh(material_bh, B, 2, H, jac);

        const Real r_mult = r_mult_qp[qp];

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; ++dof_i)
        {
            const Real bx_i = static_cast<Real>(BAe_qp(qp, dof_i));
            const Real by_i = static_cast<Real>(BAe_qp(qp + nq, dof_i));
            Re(dof_i) += JxW[qp] * (bx_i * H[0] + by_i * H[1]) * r_mult * thk;
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; ++dof_j)
            {
                const Real bx_j = static_cast<Real>(BAe_qp(qp, dof_j));
                const Real by_j = static_cast<Real>(BAe_qp(qp + nq, dof_j));
                const Real hx_j = jac[0][0] * bx_j + jac[0][1] * by_j;
                const Real hy_j = jac[1][0] * bx_j + jac[1][1] * by_j;
                Ke(dof_i, dof_j) += JxW[qp] * (bx_i * hx_j + by_i * hy_j) * r_mult * thk;
            }
        }
    }
}
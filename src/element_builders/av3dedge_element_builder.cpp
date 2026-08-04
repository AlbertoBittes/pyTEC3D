/**
 * @file av3dedge_element_builder.cpp
 * @brief ElementBuilder::build_av3dedge() — element matrices/vectors for the AV3D edge-element
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

#include <cmath>
#include <stdexcept>

#include "element_types.h"
#include "apec_m_system.h"
#include "coordinate_system.h"
#include "element_types_factory.h"
#include "variables_map.h"
#include "apec_m_logging.h"

using libMesh::DenseVector;
using libMesh::dof_id_type;
using libMesh::FEBase;
using libMesh::FEVectorBase;
using libMesh::Node;
using libMesh::Number;
using libMesh::Point;
using libMesh::Real;
using libMesh::RealGradient;
using libMesh::RealVectorValue;

void ElementBuilder::build_av3dedge()
{
    FEBase *fe = &elem_data.element_type->scalar_fe(3);
    FEVectorBase *fe_edge = &elem_data.element_type->vector_fe(0);
    FEVectorBase *fe_vec_lagrange = &elem_data.element_type->vector_fe(2);
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
    const std::vector<std::vector<RealVectorValue>> &phi_edge = fe_edge->get_phi();
    const std::vector<std::vector<RealGradient>> &curl_phi_edge = fe_edge->get_curl_phi();
    const std::vector<std::vector<RealGradient>> &phi_vec_lagrange = fe_vec_lagrange->get_phi();
    // const auto &qp_edge = fe_edge->get_xyz();
    // const auto &qp_vec_lagrange = fe_vec_lagrange->get_xyz();

    fe->reinit(elem_data.elem);
    fe_edge->reinit(elem_data.elem);
    fe_vec_lagrange->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    unsigned int n_dofs_A = fe_edge->n_shape_functions();
    unsigned int n_dofs_V = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();
    unsigned int nq_edge = fe_edge->n_quadrature_points();
    unsigned int nq_lagrange = fe_vec_lagrange->n_quadrature_points();

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

    auto build_ke = [&]()
    {
        for (unsigned int qp = 0; qp < nq_edge; qp++)
        {
            const Real w = JxW[qp];
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                {
                    Real contraction = 0.0;
                    for (unsigned int a = 0; a < 3; ++a)
                        for (unsigned int b = 0; b < 3; ++b)
                            contraction += curl_phi_edge[dof_i][qp](a) * nu_tensor[a][b] * curl_phi_edge[dof_j][qp](b);
                    Ke(dof_i, dof_j) += w * contraction;
                }
        }
    };

    auto build_aae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                        AAe(dof_i + k * n_dofs_V, dof_j) += w * (phi_vec_lagrange[dof_i * 3 + k][qp] * curl_phi_edge[dof_j][qp]);
        }
    };

    auto build_me = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            const Real w = JxW[qp];
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Me(dof_i + k * n_dofs_V, dof_j + l * n_dofs_V) += w * (phi_vec_lagrange[dof_i * 3 + k][qp] * phi_vec_lagrange[dof_j * 3 + l][qp]);
        }
    };

    auto build_bae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    BAe(qp + k * nq, dof_j) += curl_phi_edge[dof_j][qp](k);
    };

    auto build_jae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    JAe(qp + k * nq, dof_j) += phi_edge[dof_j][qp](k);
    };

    auto build_jve = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    JVe(qp + k * nq, dof_j) += dphi[dof_j][qp](k);
    };

    if (!this->nonlinear_element())
        build_ke();
    build_aae();
    build_me();
    build_bae();
    build_jae();
    build_jve();

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]))
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
        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    {
                        Real contraction = 0.0;
                        for (unsigned int a = 0; a < 3; ++a)
                            for (unsigned int b = 0; b < 3; ++b)
                                contraction += phi_edge[dof_i][qp](a) * sigma_tensor[a][b] * phi_edge[dof_j][qp](b);
                        Ce(dof_i, dof_j) += w * contraction;
                    }
            }
        };

        auto build_cav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        Real contraction = 0.0;
                        for (unsigned int a = 0; a < 3; ++a)
                            for (unsigned int b = 0; b < 3; ++b)
                                contraction += phi_edge[dof_i][qp](a) * sigma_tensor[a][b] * dphi[dof_j][qp](b);
                        Ce(dof_i, dof_j + n_dofs_A) += w * contraction;
                    }
            }
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    {
                        Real contraction = 0.0;
                        for (unsigned int a = 0; a < 3; ++a)
                            for (unsigned int b = 0; b < 3; ++b)
                                contraction += dphi[dof_i][qp](a) * sigma_tensor[a][b] * phi_edge[dof_j][qp](b);
                        Ce(dof_i + n_dofs_A, dof_j) += w * contraction;
                    }
            }
        };

        auto build_cvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        Real contraction = 0.0;
                        for (unsigned int a = 0; a < 3; ++a)
                            for (unsigned int b = 0; b < 3; ++b)
                                contraction += dphi[dof_i][qp](a) * sigma_tensor[a][b] * dphi[dof_j][qp](b);
                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += w * contraction;
                    }
            }
        };

        auto build_kav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        Real contraction = 0.0;
                        for (unsigned int a = 0; a < 3; ++a)
                            for (unsigned int b = 0; b < 3; ++b)
                                contraction += phi_edge[dof_i][qp](a) * sigma_tensor[a][b] * dphi[dof_j][qp](b);
                        Ke(dof_i, dof_j + n_dofs_A) += w * contraction;
                    }
            }
        };

        auto build_kvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        Real contraction = 0.0;
                        for (unsigned int a = 0; a < 3; ++a)
                            for (unsigned int b = 0; b < 3; ++b)
                                contraction += dphi[dof_i][qp](a) * sigma_tensor[a][b] * dphi[dof_j][qp](b);
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += w * contraction;
                    }
            }
        };

        build_caa();
        if (elem_data.element_type->get_OPTIONs()[4] == 0 &&
            static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]) == PhysicalRegionEnum::SOLID_CONDUCTOR)
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
        const Real t_local[3] = {
            this->constant_value(constant_id, CONST_T_X),
            this->constant_value(constant_id, CONST_T_Y),
            this->constant_value(constant_id, CONST_T_Z)};
        Real t_global[3] = {0.0, 0.0, 0.0};
        this->rotate_vector_for_elem(t_local, t_global);
        const Real tx = t_global[0];
        const Real ty = t_global[1];
        const Real tz = t_global[2];
        const Real r_coil = this->constant_value(constant_id, CONST_R);
        const Real sym = this->constant_value(constant_id, CONST_S);
        const Real kav_coeff = -(nc / (sc * r_coil));
        const Real kvv_coeff = (1.0 / (r_coil * vc));

        auto build_cav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real t_dot_w =
                            tx * phi_edge[dof_i][qp](0) +
                            ty * phi_edge[dof_i][qp](1) +
                            tz * phi_edge[dof_i][qp](2);
                        const Real kav = w * (kav_coeff * t_dot_w * phi[dof_j][qp]);
                        Ce(dof_i, dof_j + n_dofs_A) += kav;
                        Ce(dof_i, dof_j + n_dofs_A + n_dofs_V) -= kav;
                    }
            }
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    {
                        const Real t_dot_w =
                            tx * phi_edge[dof_j][qp](0) +
                            ty * phi_edge[dof_j][qp](1) +
                            tz * phi_edge[dof_j][qp](2);
                        const Real kav_t = w * (kav_coeff * t_dot_w * phi[dof_i][qp]);
                        Ce(dof_i + n_dofs_A, dof_j) += sym * kav_t;
                        Ce(dof_i + n_dofs_A + n_dofs_V, dof_j) -= kav_t;
                    }
            }
        };

        auto build_cvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real kvv = w * (kvv_coeff * phi[dof_i][qp] * phi[dof_j][qp]);
                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ce(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A + n_dofs_V) -= kvv;
                    }
            }
        };

        auto build_kav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real t_dot_w =
                            tx * phi_edge[dof_i][qp](0) +
                            ty * phi_edge[dof_i][qp](1) +
                            tz * phi_edge[dof_i][qp](2);
                        const Real kav = w * (kav_coeff * t_dot_w * phi[dof_j][qp]);
                        Ke(dof_i, dof_j + n_dofs_A) += kav;
                        Ke(dof_i, dof_j + n_dofs_A + n_dofs_V) -= kav;
                    }
            }
        };

        auto build_kvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                const Real w = JxW[qp];
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real kvv = w * (kvv_coeff * phi[dof_i][qp] * phi[dof_j][qp]);
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ke(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A + n_dofs_V) -= kvv;
                    }
            }
        };

        if (elem_data.element_type->get_OPTIONs()[4] == 0)
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
    default:
        apec::logging::error_and_abort("Physical region not implemented");
    }
    // The penalty values, for solution boundary trace and flux.
}

void ElementBuilder::build_av3dedge_rhs()
{
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};
    this->build_vectors(&vector_sizes, VectorsSet::RHS_VECTORS);
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    FEVectorBase *fe_edge = &elem_data.element_type->vector_fe(1);
    const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_edge->get_phi();
    const std::vector<std::vector<RealVectorValue>> &curl_phi_vec = fe_edge->get_curl_phi();
    const std::vector<Real> &JxW = fe_edge->get_JxW();

    unsigned int n_dofs_A = elem_data.elem->n_edges();
    unsigned int nq;

    Real mg[3] = {0.0, 0.0, 0.0};

    auto physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        fe_edge->reinit(elem_data.elem);
        nq = fe_edge->n_quadrature_points();
    }
    default:
    {
        const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);

        const Real mg_local[3] = {
            this->material_value(material_id, MAT_M_X),
            this->material_value(material_id, MAT_M_Y),
            this->material_value(material_id, MAT_M_Z)};
        this->rotate_vector_for_elem(mg_local, mg);
        if (mg[0] != 0.0 || mg[1] != 0.0 || mg[2] != 0.0)
        {
            fe_edge->reinit(elem_data.elem);
            nq = fe_edge->n_quadrature_points();
        }
    }
    }

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
        const Real js_local[3] = {
            js_array->GetComponent(elem_id, 0),
            js_array->GetComponent(elem_id, 1),
            js_array->GetComponent(elem_id, 2)};
        Real js[3] = {0.0, 0.0, 0.0};
        this->rotate_vector_for_elem(js_local, js);
        auto build_fa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int k = 0; k < 3; k++)
                        Fe(dof_i) += JxW[qp] * phi_vec[dof_i][qp](k) * js[k];
        };
        build_fa();

        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        if (dof_solution_system_->has_variable("VOLT") && dof_solution_system_->has_variable("Az"))
        {
            auto build_fv = [&]()
            {
                const unsigned int volt_var_num = dof_solution_system_->variable_number("VOLT");
                vtkDataArray *amps_array = vtk_grid_->GetPointData()->GetArray("AMPS");

                std::vector<dof_id_type> dof_ids_v;

                dof_solution_dof_map_->dof_indices(elem_data.elem, dof_ids_v, volt_var_num);
                const unsigned int n_dofs_V = static_cast<unsigned int>(dof_ids_v.size());

                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                {
                    const Node &node = elem_data.elem->node_ref(dof_i);
                    const vtkIdType node_id = static_cast<vtkIdType>(node.id());
                    const Real amps_value = amps_array->GetComponent(node_id, 0);
                    const dof_id_type dof_id = node.dof_number(dof_solution_system_->number(), volt_var_num, 0);
                    if (!nodal_loads_dof_cache.insert(dof_id).second)
                        continue;
                    Fe(dof_i + n_dofs_A) = amps_value;
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
        if (mg[0] != 0.0 || mg[1] != 0.0 || mg[2] != 0.0)
        {
            auto build_mg = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                    for (unsigned int k = 0; k < 3; k++)
                        for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                            Fe(dof_i) += JxW[qp] * curl_phi_vec[dof_i][qp](k) * mg[k];
            };
            build_mg();
        }
    }
    }
    return;
}

void ElementBuilder::build_av3dedge_nodal_solution()
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
    const auto &Az = sol_vectors[0];

    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    auto &nodal_vectors = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS];
    nodal_vectors.clear();
    nodal_vectors.resize(api_vars.size());
    DenseVector<Number> Bn;
    DenseVector<Number> Hn;

    unsigned int n_dofs_V = AAe.m() / 3;

    bool have_b = false;
    bool have_h = false;

    auto build_b = [&]()
    {
        DenseVector<Number> Bn_x;
        AAe.vector_mult(Bn_x, Az);
        Me.cholesky_solve(Bn_x, Bn);
        have_b = true;
    };

    auto build_h = [&]()
    {
        if (!have_b)
            build_b();
        Hn.resize(Bn.size());
        for (unsigned int dof_i = 0; dof_i < n_dofs_V; ++dof_i)
            for (unsigned int k = 0; k < 3; ++k)
                for (unsigned int l = 0; l < 3; ++l)
                    Hn(dof_i + k * n_dofs_V) += nu_tensor[k][l] * Bn(dof_i + l * n_dofs_V);
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

void ElementBuilder::build_av3dedge_element_solution()
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
    const auto &Az = sol_vectors[0];

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
        BAe.vector_mult(Bq, Az);
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
        for (unsigned int k = 0; k < 3; ++k)
            for (unsigned int l = 0; l < 3; ++l)
                He(k) += nu_tensor[k][l] * Be(l);
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
            const Real js_local[3] = {
                js_array->GetComponent(elem_id, 0),
                js_array->GetComponent(elem_id, 1),
                js_array->GetComponent(elem_id, 2)};
            Real js[3] = {0.0, 0.0, 0.0};
            this->rotate_vector_for_elem(js_local, js);
            for (unsigned int k = 0; k < 3; k++)
                Je(k) = js[k];
            break;
        }
        case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        case PhysicalRegionEnum::SOLID_CONDUCTOR:
        {
            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
            auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];

            const Real sigma_local[3] = {
                1.0 / this->material_value(material_id, MAT_RHO_X),
                1.0 / this->material_value(material_id, MAT_RHO_Y),
                1.0 / this->material_value(material_id, MAT_RHO_Z)};
            Real sigma_tensor[3][3] = {{0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0},
                                       {0.0, 0.0, 0.0}};
            this->rotate_diag_tensor_for_elem(sigma_local, sigma_tensor);

            DenseVector<Number> Vn;

            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR &&
                elem_data.element_type->get_OPTIONs()[4] == 0)
                Vn = sol_vectors[1];
            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"), 0); // solution_derivative

            const auto &A_dt = sol_vectors_dt[0];

            if (physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR ||
                elem_data.element_type->get_OPTIONs()[4] != 0)
                Vn = sol_vectors_dt[1];

            DenseVector<Number> Jq_A;
            JAe.vector_mult(Jq_A, A_dt);

            DenseVector<Number> Jq_V;
            JVe.vector_mult(Jq_V, Vn);

            Jq_A += Jq_V;

            Real j_avg[3] = {0.0, 0.0, 0.0};
            for (unsigned int l = 0; l < 3; ++l)
            {
                for (unsigned int qp = 0; qp < n_qp; ++qp)
                    j_avg[l] -= Jq_A(qp + l * n_qp);
                j_avg[l] /= n_qp;
            }
            for (unsigned int k = 0; k < 3; ++k)
                for (unsigned int l = 0; l < 3; ++l)
                    Je(k) += sigma_tensor[k][l] * j_avg[l];
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
            const Real r_coil = this->constant_value(constant_id, CONST_R);

            const Real coeff = nc / (sc * r_coil) * (Vn - En);

            for (unsigned int k = 0; k < 3; k++)
                Je(k) = coeff * t_global[k];
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

void ElementBuilder::build_av3dedge_tangent_stiffness_residual()
{
    FEBase *fe = &elem_data.element_type->scalar_fe(3);
    FEVectorBase *fe_edge = &elem_data.element_type->vector_fe(0);
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
    const std::vector<Real> &JxW = fe_edge->get_JxW();
    const std::vector<std::vector<RealVectorValue>> &phi_edge = fe_edge->get_phi();
    const std::vector<std::vector<RealGradient>> &curl_phi_edge = fe_edge->get_curl_phi();

    fe->reinit(elem_data.elem);
    fe_edge->reinit(elem_data.elem);

    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs}};
    this->build_matrices(&matrix_sizes, MatricesSet::TANGENT_STIFFNESS_MATRICES);
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};
    this->build_vectors(&vector_sizes, VectorsSet::RESIDUAL_VECTORS);

    const unsigned int n_dofs_A = fe_edge->n_shape_functions();
    const unsigned int n_dofs_V = fe->n_shape_functions();
    const unsigned int nq_edge = fe_edge->n_quadrature_points();

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const auto sol_vectors =
        this->get_variables_solution(dof_solution_system_->current_local_solution.get(),
                                     ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    const auto &Az = sol_vectors[0];

    auto &Ke = this->matrices[MatricesSet::TANGENT_STIFFNESS_MATRICES][0];
    auto &Re = this->vectors[VectorsSet::RESIDUAL_VECTORS][0];

    if (!this->nonlinear_element())
        return;

    const auto material_bh = &bh_curves_->find(material_id)->second;
    const auto bh_x = material_bh->find(0);
    const auto bh_y = material_bh->find(1);
    const auto bh_z = material_bh->find(2);
    if (bh_x == material_bh->end() || bh_y == material_bh->end() || bh_z == material_bh->end())
        apec::logging::error_and_abort("Missing BH curve component for AV3D edge tangent stiffness material");

    libMesh::Real position[3] = {0.0, 0.0, 0.0};
    bool has_cs = false;
    if (elem_data.coordinate_system)
    {
        has_cs = true;
        const Point center = elem_data.elem->vertex_average();
        position[0] = static_cast<libMesh::Real>(center(0));
        position[1] = static_cast<libMesh::Real>(center(1));
        position[2] = static_cast<libMesh::Real>(center(2));
    }

    // Nonlinear KAA block: build BAe(qp,k;dof) and compute Bq = BAe * Az (3*nq_edge entries).
    libMesh::DenseMatrix<libMesh::Number> BAe_qp;
    BAe_qp.resize(3 * nq_edge, n_dofs_A);
    for (unsigned int qp = 0; qp < nq_edge; ++qp)
        for (unsigned int k = 0; k < 3; ++k)
            for (unsigned int dof_j = 0; dof_j < n_dofs_A; ++dof_j)
                BAe_qp(qp + k * nq_edge, dof_j) += curl_phi_edge[dof_j][qp](k);

    DenseVector<Number> Bq;
    BAe_qp.vector_mult(Bq, Az);

    for (unsigned int qp = 0; qp < nq_edge; ++qp)
    {
        const Real bx_global = static_cast<Real>(Bq(qp + 0 * nq_edge));
        const Real by_global = static_cast<Real>(Bq(qp + 1 * nq_edge));
        const Real bz_global = static_cast<Real>(Bq(qp + 2 * nq_edge));
        Real b_local[3] = {
            bx_global,
            by_global,
            bz_global};
        if (has_cs)
        {
            const libMesh::Real b_global_raw[3] = {
                static_cast<libMesh::Real>(bx_global),
                static_cast<libMesh::Real>(by_global),
                static_cast<libMesh::Real>(bz_global)};
            libMesh::Real b_local_raw[3] = {0.0, 0.0, 0.0};
            elem_data.coordinate_system->vector_components_to_local_raw(
                b_global_raw,
                position,
                1,
                b_local_raw);
            b_local[0] = static_cast<Real>(b_local_raw[0]);
            b_local[1] = static_cast<Real>(b_local_raw[1]);
            b_local[2] = static_cast<Real>(b_local_raw[2]);
        }

        // Nonlinear residual term g^A(A): C_i · H(B(A)).
        const Real h_local[3] = {
            bh_x->second.evaluate_reluctivity_from_b(std::abs(b_local[0])) * b_local[0],
            bh_y->second.evaluate_reluctivity_from_b(std::abs(b_local[1])) * b_local[1],
            bh_z->second.evaluate_reluctivity_from_b(std::abs(b_local[2])) * b_local[2]};
        Real h_global[3] = {h_local[0], h_local[1], h_local[2]};
        if (has_cs)
        {
            const libMesh::Real h_local_raw[3] = {
                static_cast<libMesh::Real>(h_local[0]),
                static_cast<libMesh::Real>(h_local[1]),
                static_cast<libMesh::Real>(h_local[2])};
            libMesh::Real h_global_raw[3] = {0.0, 0.0, 0.0};
            elem_data.coordinate_system->vector_components_to_global_raw(
                h_local_raw,
                position,
                1,
                h_global_raw);
            h_global[0] = static_cast<Real>(h_global_raw[0]);
            h_global[1] = static_cast<Real>(h_global_raw[1]);
            h_global[2] = static_cast<Real>(h_global_raw[2]);
        }

        const Real d0 = bh_x->second.slope_h_vs_b(std::abs(b_local[0]));
        const Real d1 = bh_y->second.slope_h_vs_b(std::abs(b_local[1]));
        const Real d2 = bh_z->second.slope_h_vs_b(std::abs(b_local[2]));

        Real nu_diff[3][3] = {{d0, 0.0, 0.0},
                              {0.0, d1, 0.0},
                              {0.0, 0.0, d2}};
        if (has_cs)
        {
            const libMesh::Real nu_local_raw[9] = {
                d0, 0.0, 0.0,
                0.0, d1, 0.0,
                0.0, 0.0, d2};
            libMesh::Real nu_global_raw[9] = {0.0};
            elem_data.coordinate_system->rotate_tensor_to_global_raw(
                nu_local_raw,
                position,
                nu_global_raw);
            for (unsigned int a = 0; a < 3; ++a)
                for (unsigned int b = 0; b < 3; ++b)
                    nu_diff[a][b] = nu_global_raw[a * 3 + b];
        }

        for (unsigned int dof_i = 0; dof_i < n_dofs_A; ++dof_i)
        {
            Real residual_term = 0.0;
            for (unsigned int a = 0; a < 3; ++a)
                residual_term += curl_phi_edge[dof_i][qp](a) * h_global[a];
            Re(dof_i) += JxW[qp] * residual_term;

            for (unsigned int dof_j = 0; dof_j < n_dofs_A; ++dof_j)
            {
                Real tangential_term = 0.0;
                for (unsigned int a = 0; a < 3; ++a)
                    for (unsigned int b = 0; b < 3; ++b)
                        tangential_term += curl_phi_edge[dof_i][qp](a) * nu_diff[a][b] * curl_phi_edge[dof_j][qp](b);
                Ke(dof_i, dof_j) += JxW[qp] * tangential_term;
            }
        }
    }
}

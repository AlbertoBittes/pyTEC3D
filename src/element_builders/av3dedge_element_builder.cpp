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

#include <stdexcept>

#include "element_types.h"
#include "apec_m_system.h"
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
using libMesh::VectorValue;

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
    const Real nu[3] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z))};

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
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ke(dof_i, dof_j) += JxW[qp] * nu[k] * (curl_phi_edge[dof_i][qp](k) * curl_phi_edge[dof_j][qp](k));
    };

    auto build_aae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                        AAe(dof_i + k * n_dofs_V, dof_j) += JxW[qp] * (phi_vec_lagrange[dof_i * 3 + k][qp] * curl_phi_edge[dof_j][qp]);
    };

    auto build_me = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            Me(dof_i + k * n_dofs_V, dof_j + l * n_dofs_V) += JxW[qp] * (phi_vec_lagrange[dof_i * 3 + k][qp] * phi_vec_lagrange[dof_j * 3 + l][qp]);
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
        const Real sigma_matrix[3] = {
            1.0 / this->material_value(material_id, MAT_RHO_X),
            1.0 / this->material_value(material_id, MAT_RHO_Y),
            1.0 / this->material_value(material_id, MAT_RHO_Z)};
        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ce(dof_i, dof_j) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * phi_edge[dof_j][qp](k));
        };

        auto build_cav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ce(dof_i, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * dphi[dof_j][qp](k));
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ce(dof_i + n_dofs_A, dof_j) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * phi_edge[dof_j][qp](k));
        };

        auto build_cvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
        };

        auto build_kav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ke(dof_i, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * dphi[dof_j][qp](k));
        };

        auto build_kvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                        for (unsigned int k = 0; k < 3; k++)
                            Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
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
        const Real tx = this->constant_value(constant_id, CONST_T_X);
        const Real ty = this->constant_value(constant_id, CONST_T_Y);
        const Real tz = this->constant_value(constant_id, CONST_T_Z);
        const Real r_coil = this->constant_value(constant_id, CONST_R);
        const Real sym = this->constant_value(constant_id, CONST_S);

        auto build_cav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real t_dot_w =
                            tx * phi_edge[dof_i][qp](0) +
                            ty * phi_edge[dof_i][qp](1) +
                            tz * phi_edge[dof_i][qp](2);
                        const Real kav = JxW[qp] * (-(nc / (sc * r_coil)) * t_dot_w * phi[dof_j][qp]);
                        Ce(dof_i, dof_j + n_dofs_A) += kav;
                        Ce(dof_i, dof_j + n_dofs_A + n_dofs_V) -= kav;
                    }
        };

        auto build_cva = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    {
                        const Real t_dot_w =
                            tx * phi_edge[dof_j][qp](0) +
                            ty * phi_edge[dof_j][qp](1) +
                            tz * phi_edge[dof_j][qp](2);
                        const Real kav_t = JxW[qp] * (-(nc / (sc * r_coil)) * t_dot_w * phi[dof_i][qp]);
                        Ce(dof_i + n_dofs_A, dof_j) += sym * kav_t;
                        Ce(dof_i + n_dofs_A + n_dofs_V, dof_j) -= kav_t;
                    }
        };

        auto build_cvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real kvv = JxW[qp] * ((1.0 / (r_coil * vc)) * phi[dof_i][qp] * phi[dof_j][qp]);
                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ce(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A + n_dofs_V) -= kvv;
                    }
        };

        auto build_kav = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real t_dot_w =
                            tx * phi_edge[dof_i][qp](0) +
                            ty * phi_edge[dof_i][qp](1) +
                            tz * phi_edge[dof_i][qp](2);
                        const Real kav = JxW[qp] * (-(nc / (sc * r_coil)) * t_dot_w * phi[dof_j][qp]);
                        Ke(dof_i, dof_j + n_dofs_A) += kav;
                        Ke(dof_i, dof_j + n_dofs_A + n_dofs_V) -= kav;
                    }
        };

        auto build_kvv = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real kvv = JxW[qp] * ((1.0 / (r_coil * vc)) * phi[dof_i][qp] * phi[dof_j][qp]);
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ke(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A + n_dofs_V) -= kvv;
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
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        FEVectorBase *fe_edge = &elem_data.element_type->vector_fe(1);

        const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_edge->get_phi();

        const std::vector<Real> &JxW = fe_edge->get_JxW();
        fe_edge->reinit(elem_data.elem);

        this->prepare_element_matrix_building();

        unsigned int n_dofs_A = fe_edge->n_shape_functions();
        unsigned int nq = fe_edge->n_quadrature_points();

        vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
        const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
        Real js[3] = {
            js_array->GetComponent(elem_id, 0),
            js_array->GetComponent(elem_id, 1),
            js_array->GetComponent(elem_id, 2)};
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
                const unsigned int az_var_num = dof_solution_system_->variable_number("Az");
                vtkDataArray *amps_array = vtk_grid_->GetPointData()->GetArray("AMPS");

                unsigned int n_dofs_A = elem_data.elem->n_dofs(ApecMSystemBase::DOF_SOLUTION_SYSTEM, az_var_num);
                unsigned int n_dofs_V = elem_data.elem->n_dofs(ApecMSystemBase::DOF_SOLUTION_SYSTEM, volt_var_num);

                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                {
                    const Node &node = elem_data.elem->node_ref(dof_i);
                    const vtkIdType node_id = static_cast<vtkIdType>(node.id());
                    const Real amps_value = amps_array->GetComponent(node_id, 0);
                    const dof_id_type dof_id = node.dof_number(ApecMSystemBase::DOF_SOLUTION_SYSTEM, volt_var_num, 0);
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
    return;
}

void ElementBuilder::build_av3dedge_nodal_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::NODAL_SOLUTION);
    if (api_vars.empty())
        return;

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const Real nu[3] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z))};

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
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                Hn(dof_i + k * n_dofs_V) += Bn(dof_i + k * n_dofs_V) * nu[k];
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
    const Real nu[3] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Z))};

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
        for (unsigned int k = 0; k < 3; k++)
            He(k) = Be(k) * nu[k];
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
        {
            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
            auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];

            const Real sigma_matrix[3] = {
                1.0 / this->material_value(material_id, MAT_RHO_X),
                1.0 / this->material_value(material_id, MAT_RHO_Y),
                1.0 / this->material_value(material_id, MAT_RHO_Z)};

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

            for (unsigned int k = 0; k < 3; k++)
            {
                for (unsigned int qp = 0; qp < n_qp; qp++)
                    Je(k) -= Jq_A(qp + k * n_qp);
                Je(k) = Je(k) * sigma_matrix[k] / n_qp;
            }
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
            const VectorValue t_vec(
                this->constant_value(constant_id, CONST_T_X),
                this->constant_value(constant_id, CONST_T_Y),
                this->constant_value(constant_id, CONST_T_Z));
            const Real r_coil = this->constant_value(constant_id, CONST_R);

            const Real coeff = nc / (sc * r_coil) * (Vn - En);

            for (unsigned int k = 0; k < 3; k++)
                Je(k) = coeff * t_vec(k);
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

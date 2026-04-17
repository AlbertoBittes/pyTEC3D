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

// Temporary placeholder for build_av2dnodal
void ElementBuilder::build_av2dnodal()
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
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[0].get());
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    const std::vector<Point> &xyz = fe->get_xyz();

    fe->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    //   const unsigned int n_sf =
    // FEInterface::n_dofs(cfe->get_fe_type(), elem);
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const dof_id_type constant_id = static_cast<dof_id_type>(this->elem_data.real_constant_num);
    const Real nu[2] = {
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_X)),
        1 / (kMu0 * this->material_value(material_id, MAT_MUR_Y))};

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ke];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ce];
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalMe];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][kNodalAAe];
    auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];
    auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJVe];
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementBAe];
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    unsigned int axisymm = elem_data.element_type->get_OPTIONs()[3];

    Real r_sign = pow(-1, axisymm);

    Real thk = this->constant_value(constant_id, CONST_THK);
    if (thk == 0.0)
        thk = 1.0;

    auto build_ke = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(2 * libMesh::pi * r, axisymm);
            Real r_sq = pow(r, 2);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ke(dof_i, dof_j) +=
                        JxW[qp] * (dphi[dof_i][qp](0) * dphi[dof_j][qp](0) * nu[0] + dphi[dof_i][qp](1) * dphi[dof_j][qp](1) * nu[1] + axisymm * phi[dof_i][qp] * phi[dof_j][qp] * nu[1] / r_sq) * r_mult * thk;
        }
    };

    auto build_bae = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            // Here maybe is more suited the centroid radius
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));

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
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(2 * libMesh::pi * r, axisymm);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    AAe(dof_i, dof_j) += JxW[qp] * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](1)) * r_mult * thk;
        }

        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(2 * libMesh::pi * r, axisymm);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    AAe(dof_i + n_dofs_var, dof_j) -=
                        JxW[qp] * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](0) - axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r) * r_mult * thk;
        }
    };

    auto build_me = [&]()
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(2 * libMesh::pi * r, axisymm);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
        }

        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(2 * libMesh::pi * r, axisymm);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Me(dof_i + n_dofs_var, dof_j + n_dofs_var) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
        }
    };

    build_ke();
    build_bae();
    build_jae();
    build_aae();
    build_me();

    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto build_fa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    Fe(dof_i) += JxW[qp] * phi[dof_i][qp] * r_mult * thk; // * js_loads[elem->id()][2];
            }
        };
        build_fa();
        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        const Real sigma = 1.0 / this->material_value(material_id, MAT_RHO_Z);

        auto build_caa = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
            }
            Ce.scale(sigma);
        };

        build_caa();

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR)
        {
            auto build_coupling_timeintegrated = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                {
                    Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                    Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                        for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        {
                            Ce(dof_i, dof_j + n_dofs_var) -=
                                JxW[qp] * sigma * (phi[dof_i][qp] * phi[dof_j][qp]);
                            Ce(dof_i + n_dofs_var, dof_j) -=
                                JxW[qp] * sigma * (phi[dof_i][qp] * phi[dof_j][qp]);
                            Ce(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                                JxW[qp] * sigma * (phi[dof_i][qp] * phi[dof_j][qp]) / r_mult / thk;
                        }
                }
            };
            auto build_coupling_notimeintegrated = [&]()
            {
                for (unsigned int qp = 0; qp < nq; qp++)
                {
                    Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                    Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                        for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        {
                            Ke(dof_i, dof_j + n_dofs_var) -=
                                JxW[qp] * sigma * (phi[dof_i][qp] * phi[dof_j][qp]);
                            Ce(dof_i + n_dofs_var, dof_j) -=
                                JxW[qp] * sigma * (phi[dof_i][qp] * phi[dof_j][qp]);
                            Ke(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                                JxW[qp] * sigma * (phi[dof_i][qp] * phi[dof_j][qp]) / r_mult / thk;
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
        const Real tz = this->constant_value(constant_id, CONST_T_Z);
        const Real r_coil = this->constant_value(constant_id, CONST_R);
        const Real sym = this->constant_value(constant_id, CONST_S);

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL elements.");

        const Real kav_scale = -(nc / (sc * r_coil)) * tz;
        const Real kvv_scale = 1.0 / (r_coil * vc);

        auto build_coil_coupling_timeintegrated = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
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
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult * thk;
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
        libmesh_error_msg("Physical region not implemented");
    }
}

// Temporary placeholder for build_av2dnodal
void ElementBuilder::build_av2dnodal_rhs()
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
    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_OPTIONs()[1]);

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        if (physical_region != PhysicalRegionEnum::STRANDED_CONDUCTOR)
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

        vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
        const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
        const Real js_z = js_array->GetComponent(elem_id, 2);

        auto build_fa = [&]()
        {
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i) *= js_z;
        };
        build_fa();
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    return;
}

void ElementBuilder::build_av2dnodal_nodal_solution()
{
    if (!this->compute_nodal_solution_vectors)
        return;

    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::NODAL_SOLUTION);
    if (api_vars.empty())
        return;

    const dof_id_type material_id = static_cast<dof_id_type>(this->elem_data.material_num);
    const Real nu = 1 / (kMu0 * this->material_value(material_id, MAT_MUR_X));

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
        for (unsigned int dof_i = 0; dof_i < Bn.size(); dof_i++)
            Hn(dof_i) += Bn(dof_i) * nu;
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
    const Real nu = 1 / (kMu0 * this->material_value(material_id, MAT_MUR_X));

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
        for (unsigned int k = 0; k < 2; k++)
            He(k) = Be(k) * nu;
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
            Je(0) = js_array->GetComponent(elem_id, 2);
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

            const Real sigma = 1.0 / this->material_value(material_id, MAT_RHO_Z);

            const auto &An_dt = sol_vectors_dt[0];

            DenseVector<Number> Jq_A(n_qp);
            JAe.vector_mult(Jq_A, An_dt);

            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(0) -= Jq_A(qp);
            Je(0) = Je(0) * sigma / n_qp;

            if (physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR)
                break;

            Real thk = this->constant_value(constant_id, CONST_THK);
            if (thk == 0.0)
                libmesh_error_msg("THK is zero for SOLID_CONDUCTOR elements.");
            Real V_factor = sigma * Vn / thk;
            if (axisymm)
            {
                FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[1].get());
                const std::vector<Point> &xyz = fe->get_xyz();
                fe->reinit(elem_data.elem);
                Real Jv = 0.0;
                for (unsigned int qp = 0; qp < JAe.m(); qp++)
                    for (unsigned int dof_i = 0; dof_i < JAe.n(); dof_i++)
                        Jv += JAe(qp, dof_i) / (2 * libMesh::pi * xyz[qp](0));
                V_factor *= Jv / n_qp;
            }
            Je(0) += V_factor;

            vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
            const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
            Je(0) += js_array->GetComponent(elem_id, 2);
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
            const Real tz = this->constant_value(constant_id, CONST_T_Z);
            const Real r_coil = this->constant_value(constant_id, CONST_R);

            if (sc == 0.0 || r_coil == 0.0)
                libmesh_error_msg("Invalid real constants (SC, R) for STRANDED_COIL elements.");
            Je(0) = (nc / (sc * r_coil)) * tz * (Vn - En);

            vtkDataArray *js_array = vtk_grid_->GetCellData()->GetArray("Js");
            const vtkIdType elem_id = static_cast<vtkIdType>(elem_data.elem->id());
            Je(0) += js_array->GetComponent(elem_id, 2);
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

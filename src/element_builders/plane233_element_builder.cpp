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
// Temporary placeholder for build_plane233
void ElementBuilder::build_plane233()
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

    const auto nu_matrix = this->get_mu_matrix(this->elem_data.material_num);
    const auto sigma_matrix = this->get_sigma_matrix(this->elem_data.material_num);

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

    Real thk = this->get_constants(this->elem_data.real_constant_num)[0];
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
                        JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp] + axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r_sq) * r_mult * thk;
        }
        Ke.scale(nu_matrix[0]);
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
        Ce.scale(sigma_matrix[0]);
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
    {
        build_caa();
        break;
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i) += JxW[qp] * phi[dof_i][qp]; // * js_loads[elem->id()][2];
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        build_caa();

        auto build_coupling = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        Ce(dof_i, dof_j + n_dofs_var) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        Ce(dof_i + n_dofs_var, dof_j) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        Ce(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) / r_mult / thk;
                    }
            }
        };
        build_coupling();
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        build_caa();

        auto build_coupling = [&]()
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(2 * libMesh::pi * r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        Ke(dof_i, dof_j + n_dofs_var) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        Ce(dof_i + n_dofs_var, dof_j) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]);
                        Ke(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) / r_mult / thk;
                    }
            }
        };
        build_coupling();
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        const auto rc = this->get_constants(this->elem_data.real_constant_num);

        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real rad = rc[4];
        Real rad_factor = 1.0;
        if (axisymm)
            rad_factor = 2 * libMesh::pi * rad;
        const Real vc = rad_factor * sc * thk;
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL_SOLENOIDAL elements.");

        const Real kav_scale = -(nc / (sc * r_coil)) * tz;
        const Real kvv_scale = 1.0 / (r_coil * vc);

        auto build_coil_coupling = [&]()
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

                        // Non-time-integrated A-VOLT-EMF formulation (OPTION(2)=0/1).
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
        build_coil_coupling();
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        const auto rc = this->get_constants(this->elem_data.real_constant_num);

        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real rad = rc[4];
        Real rad_factor = 1.0;
        if (axisymm)
            rad_factor = 2 * libMesh::pi * rad;
        const Real vc = rad_factor * sc * thk;
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL_SOLENOIDAL elements.");

        const Real kav_scale = -(nc / (sc * r_coil)) * tz;
        const Real kvv_scale = 1.0 / (r_coil * vc);

        auto build_coil_coupling = [&]()
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

                        // Non-time-integrated A-VOLT-EMF formulation (OPTION(2)=0/1).
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
        build_coil_coupling();
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}

// Temporary placeholder for build_plane233
void ElementBuilder::build_plane233_rhs()
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
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
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
        if (dof_solution_system_->has_variable(volt_var_name) && dof_solution_system_->has_variable("A_z"))
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
                    Fe(dof_i + n_dofs_var) = amps_value;
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

void ElementBuilder::build_plane233_solution_vectors()
{
    this->build_plane233_nodal_solution();
    this->build_plane233_element_solution();
}

void ElementBuilder::build_plane233_nodal_solution()
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
            Hn(dof_i) += Bn(dof_i) * nu_matrix[0];
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

void ElementBuilder::build_plane233_element_solution()
{
    auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::ELEMENT_SOLUTION);
    if (api_vars.empty())
        return;

    std::vector<dof_id_type> dof_indices;

    const auto nu_matrix = this->get_mu_matrix(this->elem_data.material_num);

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
            He(k) = Be(k) * nu_matrix[0];
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
        case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
        {
            Real Vn = 0.0;
            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
                Vn = sol_vectors[1][0];

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR)
                Vn = sol_vectors_dt[1][0];

            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][kElementJAe];

            const auto sigma_matrix = this->get_sigma_matrix(this->elem_data.material_num);

            const auto &An_dt = sol_vectors_dt[0];

            DenseVector<Number> Jq_A(n_qp);
            JAe.vector_mult(Jq_A, An_dt);

            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(0) -= Jq_A(qp);
            Je(0) = Je(0) * sigma_matrix[0] / n_qp;

            if (physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR)
                break;

            Real thk = this->get_constants(this->elem_data.real_constant_num)[0];
            if (thk == 0.0)
                libmesh_error_msg("THK is zero for SOLID_CONDUCTOR_SOLENOIDAL elements.");
            Real V_factor = sigma_matrix[0] * Vn / thk;
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
            break;
        }
        case PhysicalRegionEnum::STRANDED_COIL:
        case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
        {
            Real Vn = 0.0;
            Real En = 0.0;

            if (physical_region == PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL)
            {
                Vn = sol_vectors[1][0];
                En = sol_vectors[2][0];
            }

            const auto sol_vectors_dt =
                get_variables_solution(dof_solution_system_->request_vector("solution_time_derivative"),
                                       ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

            if (physical_region == PhysicalRegionEnum::STRANDED_COIL)
            {
                Vn = sol_vectors_dt[1][0];
                En = sol_vectors_dt[2][0];
            }

            const auto rc = this->get_constants(this->elem_data.real_constant_num);
            const Real sc = rc[1];
            const Real nc = rc[2];
            const Real tz = rc[7];
            const Real r_coil = rc[8];

            if (sc == 0.0 || r_coil == 0.0)
                libmesh_error_msg("Invalid real constants (SC, R) for STRANDED_COIL_SOLENOIDAL elements.");
            Je(0) = (nc / (sc * r_coil)) * tz * (Vn - En);
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

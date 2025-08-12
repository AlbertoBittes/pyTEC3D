#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "loads_object.h"
#include "apec_m_system_objects.h"
#include "element_types_factory.h"

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

    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_var = fe->n_shape_functions();
    unsigned int n_dofs_A = n_dofs_var * 3;
    unsigned int nq = fe->n_quadrature_points();

    elem_builder_misc_data.materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = elem_builder_misc_data.materials_object->get_mu_matrix();

    auto &Ke = this->matrices[0][0];
    auto &Ce = this->matrices[0][1];
    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &JAe = this->matrices[2][0];
    auto &JVe = this->matrices[2][1];
    auto &Fe = this->vectors[0][0];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                for (unsigned int l = 0; l < 3; l++)
                    for (unsigned int k = 0; k < 3; k++)
                    {
                        Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / mu_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
                        AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
                        Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);
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
        elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();
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
        //         for (unsigned int k = 0; k < 3; k++)
        //             Fe(dof_i + k * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
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

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &Fe = this->vectors[0][0];

        unsigned int n_dofs_var = elem_data.elem->n_nodes();

        Real js[3];

        dof_id_type dof_jx = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsx_var_num, 0);
        dof_id_type dof_jy = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsy_var_num, 0);
        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        js[0] = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jx);
        js[1] = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jy);
        js[2] = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);

        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i + k * n_dofs_var) *= js[k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        // Here you will need to implement loads with F
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}

void ElementBuilder::build_solid97_nodal_solution_vector()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

    // Get Ae
    DenseVector<Number> Ae;
    for (unsigned int k = 0; k < 3; k++)
    {
        Ae.append(solution_vector_support.sol_vectors[0][k]);
    }

    auto &AAe = this->matrices[1][1];
    DenseVector<Number> Be_q(AAe.m());
    AAe.vector_mult(Be_q, Ae);
    auto &Me = this->matrices[1][0];
    auto &Be = this->vectors[2][0];
    Be.resize(Me.m());
    Me.cholesky_solve(Be_q, Be);
    // fill_solution_vector_support(main_sys->old_local_solution);

    // // Get Ae
    // DenseVector<Number> Ae_old;
    // for (unsigned int k = 0; k < 3; k++)
    // {
    //     Ae_old.append(solution_vector_support.sol_vectors[k]);
    // }
    // Ae_dt.add(-1.0, Ae_old);
    // Ae_dt.scale(1.0 / time_step->get_active_substep()->get_dt());
}

void ElementBuilder::build_solid97_element_solution_vector()
{
    // Materials could be retrieved in another point

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &J_elem = this->vectors[3][0];
        dof_id_type dof_jx = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsx_var_num, 0);
        dof_id_type dof_jy = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsy_var_num, 0);
        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        J_elem(0) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jx);
        J_elem(1) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jy);
        J_elem(2) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);
        if (!this->elem_builder_misc_data.ramped)
        {
            break;
        }
        else
        {
            if (this->elem_builder_misc_data.old_loads_object != nullptr)
            {
                auto &J_elem_old = this->vectors[4][0];
                // Apply ramping to the current load
                J_elem_old(0) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jx);
                J_elem_old(1) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jy);
                J_elem_old(2) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jz);
            }
            this->ramp_vectors(3, 4, 0, 0);
        }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (elem_builder_misc_data.transient_effect)
        {
            auto &J_elem = this->vectors[3][0];

            elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);

            auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();

            unsigned int solv_var_vector_size = solution_vector_support.sol_vectors[0].size();

            // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

            fill_solution_vector_support(main_sys->older_local_solution, 0); // solution_derivative

            DenseVector<Number> Ae_dt;
            for (unsigned int k = 0; k < 3; k++)
            {
                Ae_dt.append(solution_vector_support.sol_vectors[0][k]);
            }

            auto &JAe = this->matrices[2][0];
            unsigned int j_qp_size = JAe.m();
            DenseVector<Number> J_elem_A(j_qp_size);
            JAe.vector_mult(J_elem_A, Ae_dt);

            DenseVector<Number> J_elem_V(j_qp_size);

            auto &JVe = this->matrices[2][1];
            unsigned int volt_var_num = elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::VOLT].system_var_num;
            auto &Ve_dt = solution_vector_support.sol_vectors[0][volt_var_num];
            JVe.vector_mult(J_elem_V, Ve_dt);
            J_elem_A += J_elem_V;

            unsigned int n_qp = j_qp_size / 3;

            for (unsigned int k = 0; k < 3; k++)
            {
                for (unsigned int qp = 0; qp < n_qp; qp++)
                    J_elem(k) -= J_elem_A(qp + k * n_qp);
                J_elem(k) = J_elem(k) / rs_matrix[k] / n_qp;
            }
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}

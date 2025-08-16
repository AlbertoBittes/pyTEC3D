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
#include "element_loads_object.h"
#include "apec_m_system_objects.h"
#include "element_types_factory.h"
// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13()
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

    elem_builder_misc_data.materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = elem_builder_misc_data.materials_object->get_mu_matrix();

    auto &Ke = this->matrices[0][0];
    auto &Ce = this->matrices[0][1];
    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &JAe = this->matrices[2][0];
    auto &Fe = this->vectors[0][0];

    unsigned int axisymm = elem_data.element_type->keyopts[3];

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
            {
                // Here I should modify the matrices someway for axisymm elements
                AAe(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * dphi[dof_j][qp](1)) * pow(r, axisymm);
                AAe(dof_i + n_dofs_var, dof_j) -= JxW[qp] * (phi[dof_i][qp] * dphi[dof_j][qp](0) + axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r) * pow(r, axisymm);
                Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * pow(r, axisymm);
                Me(dof_i + n_dofs_var, dof_j + n_dofs_var) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * pow(r, axisymm);
            }

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                Ke(dof_i, dof_j) += JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp] + axisymm * phi[dof_i][qp] * phi[dof_j][qp] / pow(r, 2)) * pow(r, axisymm);
    }
    Ke.scale(1.0 / mu_matrix[0]);
    PhysicalRegionEnum physical_region = elem_data.element_type->get_physical_region();

    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i) += JxW[qp] * phi[dof_i][qp] * pow(r, axisymm); // * js_loads[elem->id()][2];
        }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                // matrices[4](0, i) -= scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[i][qp] / scalar_shape_functions->phi[0].size();
                JAe(qp, dof_i) += phi[dof_i][qp];

                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp] * pow(r, axisymm));
            }
        }
        Ce.scale(1.0 / rs_matrix[0]);
        break;
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i) += JxW[qp] * phi[dof_i][qp]; // * js_loads[elem->id()][2];
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}

// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13_rhs()
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
    // Now start logging the element matrix computation

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &Fe = this->vectors[0][0];

        unsigned int n_dofs_var = elem_data.elem->n_nodes();

        Real js;

        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        js = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) *= js;
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

void ElementBuilder::build_plane13_nodal_solution_vector()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

    auto &AAe = this->matrices[1][1];
    DenseVector<Number> Be_q(AAe.m());
    AAe.vector_mult(Be_q, solution_vector_support.sol_vectors[0][0]);
    auto &Me = this->matrices[1][0];
    auto &Be = this->vectors[2][0];
    Be.resize(Me.m());
    Me.cholesky_solve(Be_q, Be);
}

void ElementBuilder::build_plane13_element_solution_vector()
{
    // Materials could be retrieved in another point

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {

        auto &J_elem = this->vectors[3][0];
        dof_id_type dof_jz = elem_data.elem->dof_number(elem_builder_misc_data.loads_sys_number, elem_builder_misc_data.jsz_var_num, 0);

        J_elem(0) = elem_builder_misc_data.loads_map->current_local_solution->el(dof_jz);
        if (!this->elem_builder_misc_data.ramped)
        {
            break;
        }
        else
        {
            if (this->elem_builder_misc_data.old_element_loads_object != nullptr)
            {
                auto &J_elem_old = this->vectors[4][0];
                // Apply ramping to the current load
                J_elem_old(0) = this->elem_builder_misc_data.loads_map->old_local_solution->el(dof_jz);
            }
            // else not needed because resizing already zeros the vector
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
            // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());

            fill_solution_vector_support(main_sys->request_vector(0), 0); // solution_derivative

            auto &JAe = this->matrices[2][0];

            unsigned int j_qp_size = JAe.m();

            DenseVector<Number> J_elem_A(j_qp_size);
            JAe.vector_mult(J_elem_A, solution_vector_support.sol_vectors[0][0]);

            for (unsigned int qp = 0; qp < j_qp_size; qp++)
                J_elem(0) -= J_elem_A(qp);

            J_elem(0) = J_elem(0) / rs_matrix[0] / j_qp_size;
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}
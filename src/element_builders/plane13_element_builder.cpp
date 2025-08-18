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

    // SOmeway you could optimize this by skipping this, is very useles and very rarely used
    auto &Fe = this->vectors[0][0];

    unsigned int n_dofs_var = elem_data.elem->n_nodes();

    std::vector<dof_id_type> dof_indices_nodal_loads;
    std::vector<Real> nodal_loads;
    elem_builder_misc_data.dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_nodal_loads);
    elem_builder_misc_data.nodal_loads_map->current_local_solution->get(dof_indices_nodal_loads, nodal_loads);

    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        Fe(dof_i) = nodal_loads[dof_i];

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

        js = elem_builder_misc_data.element_loads_map->current_local_solution->el(dof_jz);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) *= js;
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    for (dof_id_type &dof_value : dof_indices_nodal_loads)
        elem_builder_misc_data.nodal_loads_map->current_local_solution->set(dof_value, 0);

    return;
}

void ElementBuilder::build_plane13_solution_vectors()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    auto &Bn = this->vectors[2][0];
    auto &Hn = this->vectors[2][1];
    auto &Be = this->vectors[2][3];
    auto &He = this->vectors[2][4];

    auto &Me = this->matrices[1][0];
    auto &AAe = this->matrices[1][1];
    auto &BAe = this->matrices[1][4];

    elem_builder_misc_data.materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = elem_builder_misc_data.materials_object->get_mu_matrix();

    unsigned int n_dofs_var = elem_data.elem->n_nodes();
    unsigned int n_qp = BAe.m();

    auto &An = solution_vector_support.sol_vectors[0][0];
    // Compute B in nodes
    DenseVector<Number> Bn_x(AAe.m());
    AAe.vector_mult(Bn_x, An);
    Bn.resize(Me.m());
    Me.cholesky_solve(Bn_x, Bn);

    // Compute H in nodes
    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        Hn(dof_i) += Bn(dof_i) * mu_matrix[0];

    // Compute B in elements
    DenseVector<Number> Bq(2 * n_qp);
    BAe.vector_mult(Bq, An);
    for (unsigned int k = 0; k < 2; k++)
    {
        for (unsigned int qp = 0; qp < n_qp; qp++)
            Be(k) += Bq(qp + k * n_qp);
        Be(k) = Be(k) / n_qp;
    }

    // Compute H in elements
    for (unsigned int k = 0; k < 2; k++)
        He(k) = Be(k) / mu_matrix[0];

    // Compute J in elements
    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        auto &Je = this->vectors[2][2];
        auto &JHe = this->matrices[1][3];

        DenseVector<Number> Jq(n_qp);
        JHe.vector_mult(Jq, Hn);

        for (unsigned int qp = 0; qp < n_qp; qp++)
            Je(0) += Jq(qp);
        Je(0) = Je(0) / mu_matrix[0] / n_qp;

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (elem_builder_misc_data.transient_effect)
        {
            auto &JAe = this->matrices[1][2];
            auto &Je = this->vectors[2][2];

            elem_builder_misc_data.materials_object->fill_rs_matrix(this->elem_data.material_num);
            auto &rs_matrix = elem_builder_misc_data.materials_object->get_rs_matrix();

            fill_solution_vector_support(main_sys->older_local_solution, 0); // solution_derivative

            // Get derivative of A DOF Solutions (dA/dt, dV/dt)
            auto &An_dt = solution_vector_support.sol_vectors[0][0];

            // Compute dA/dt in quadrature points
            DenseVector<Number> Jq_A(n_qp);
            JAe.vector_mult(Jq_A, An_dt);

            // average on quadrature points and scale by conductivity
            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(0) -= Jq_A(qp);
            Je(0) = Je(0) / rs_matrix[0] / n_qp;
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
}
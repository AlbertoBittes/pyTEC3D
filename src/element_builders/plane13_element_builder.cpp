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

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = materials_object->get_mu_matrix();

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][PLANE13DOFMatrices::Ke_PLANE13];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][PLANE13DOFMatrices::Ce_PLANE13];
    auto &Me = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::Me_PLANE13];
    auto &AAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::AAe_PLANE13];
    auto &JAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::JAe_PLANE13];
    // auto &JHe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::JHe_PLANE13];
    auto &BAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::BAe_PLANE13];
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    unsigned int axisymm = elem_data.element_type->keyopts[3];

    Real r_sign = pow(-1, axisymm);

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        Real r_mult = pow(r, axisymm);
        Real r_sq = pow(r, 2);
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        {
            // Here I should modify the matrices someway for axisymm elements
            // Also, the pow here is not always needed, is not there when you dont integrate
            BAe(qp, dof_i) += r_sign * dphi[dof_i][qp](1);
            BAe(qp + nq, dof_i) -= (r_sign * dphi[dof_i][qp](0) - axisymm * phi[dof_i][qp] / r);
            // JHe(qp, dof_i) += dphi[dof_i][qp](0);
            // JHe(qp, dof_i + n_dofs_var) += (pow(-1, axisymm) * dphi[dof_i][qp](1) + axisymm * phi[dof_i][qp] / r);
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
            {
                // Here I should modify the matrices someway for axisymm elements
                AAe(dof_i, dof_j) += JxW[qp] * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](1)) * r_mult;
                AAe(dof_i + n_dofs_var, dof_j) -= JxW[qp] * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](0) - axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r) * r_mult;
                Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
                Me(dof_i + n_dofs_var, dof_j + n_dofs_var) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
            }
        }
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                Ke(dof_i, dof_j) += JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp] + axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r_sq) * r_mult;
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
        materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = materials_object->get_rs_matrix();
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

    std::vector<dof_id_type> dof_indices_loads;
    std::vector<Real> loads_vec;

    System *nodal_loads_system = systems[Globals::NODAL_LOADS];
    System *element_loads_system = systems[Globals::ELEMENT_LOADS];

    auto dof_map_nodal_loads = dof_maps[Globals::NODAL_LOADS];
    auto dof_map_element_loads = dof_maps[Globals::ELEMENT_LOADS];

    NumericVector<Number> *nodal_loads_vector = nodal_loads_system->request_vector(2);
    NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        dof_map_element_loads->dof_indices(elem_data.elem, dof_indices_loads);
        element_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) *= loads_vec[0];

        // Apply nodal CSG loads
        dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads);
        nodal_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            if (loads_vec[dof_i] != 0)
                Fe(dof_i) = loads_vec[dof_i];

        for (dof_id_type &dof_value : dof_indices_loads)
            nodal_loads_vector->set(dof_value, 0);

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        // Apply nodal CSG loads
        dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads);
        nodal_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) = loads_vec[dof_i];

        for (dof_id_type &dof_value : dof_indices_loads)
            nodal_loads_vector->set(dof_value, 0);

        break;
        // unsigned int volt_var_num = this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::VOLT].system_var_num;
        // unsigned int amp_var_num = this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::AMPS].system_var_num;

        // // Apply nodal AMPS loads
        // dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads);
        // nodal_loads_map->current_local_solution->get(dof_indices_loads, loads_vec);

        // for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //     Fe(dof_i + volt_var_num * n_dofs_var) = loads_vec[dof_i + amp_var_num * n_dofs_var];

        // for (dof_id_type &dof_value : dof_indices_loads)
        //     nodal_loads_map->current_local_solution->set(dof_value, 0);

        // break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    return;
}

void ElementBuilder::build_plane13_solution_vectors()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    DenseVector<Number> Bn;
    DenseVector<Number> Hn;
    DenseVector<Number> Be(2);
    DenseVector<Number> He(2);
    DenseVector<Number> Je(1);

    std::vector<dof_id_type> dof_indices;

    auto &Me = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::Me_PLANE13];
    auto &AAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::AAe_PLANE13];
    auto &BAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::BAe_PLANE13];

    ApecMLinearSystem *main_sys = dynamic_cast<ApecMLinearSystem *>(systems[Globals::DOF_SOLUTION_SYSTEM]);
    System *nodal_solution_system = systems[Globals::NODAL_SOLUTION];
    System *element_solution_system = systems[Globals::ELEMENT_SOLUTION];

    this->fill_solution_vector_support(main_sys->current_local_solution.get(), Globals::DOF_SOLUTION_SYSTEM);

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = materials_object->get_mu_matrix();

    unsigned int n_dofs_var = elem_data.elem->n_nodes();
    unsigned int n_qp = BAe.m() / 2;

    auto &An = solution_vector_support.sol_vectors[0];
    // Compute B in nodes
    DenseVector<Number> Bn_x;
    AAe.vector_mult(Bn_x, An);
    // Bn.resize(Me.m());
    Me.cholesky_solve(Bn_x, Bn);

    // Compute H in nodes
    Hn.resize(Bn.size());
    for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        Hn(dof_i) += Bn(dof_i) / mu_matrix[0];

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
        // auto &JHe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::JHe_PLANE13];

        // DenseVector<Number> Jq(n_qp);
        // JHe.vector_mult(Jq, Hn);

        // for (unsigned int qp = 0; qp < n_qp; qp++)
        //     Je(0) += Jq(qp);
        // Je(0) = Je(0) / n_qp;
        dof_maps[Globals::ELEMENT_LOADS]->dof_indices(elem_data.elem, dof_indices);
        NumericVector<Number> *element_loads_vector = systems[Globals::ELEMENT_LOADS]->request_vector(2);
        element_loads_vector->get(dof_indices, Je.get_values());

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (this->context<TimeStep>()->get_transient_effect())
        {
            auto &JAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::JAe_PLANE13];

            materials_object->fill_rs_matrix(this->elem_data.material_num);
            auto &rs_matrix = materials_object->get_rs_matrix();

            fill_solution_vector_support(main_sys->older_local_solution, Globals::DOF_SOLUTION_SYSTEM); // solution_derivative

            // Get derivative of A DOF Solutions (dA/dt, dV/dt)
            auto &An_dt = solution_vector_support.sol_vectors[0];

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
    // Set B and H in nodes
    this->dof_maps[Globals::NODAL_SOLUTION]->dof_indices(elem_data.elem, dof_indices);
    Bn.append(Hn);
    for (unsigned int i = 0; i < dof_indices.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((nodal_solution_system->solution->first_local_index() <= dof_index) &&
            (dof_index < nodal_solution_system->solution->last_local_index()))
            nodal_solution_system->solution->set(dof_index, Bn(i));
    }

    this->dof_maps[Globals::ELEMENT_SOLUTION]->dof_indices(elem_data.elem, dof_indices);
    Je.append(Be);
    Je.append(He);
    for (unsigned int i = 0; i < dof_indices.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((element_solution_system->solution->first_local_index() <= dof_index) &&
            (dof_index < element_solution_system->solution->last_local_index()))
            element_solution_system->solution->set(dof_index, Je(i));
    }
}
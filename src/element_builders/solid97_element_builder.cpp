#include "element_builder.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/boundary_info.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "apec_m_system_objects.h"
#include <libmesh/boundary_info.h>
#include "element_types_factory.h"
#include "globals.h"

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

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = materials_object->get_mu_matrix();

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][SOLID97DOFMatrices::Ke_SOLID97];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][SOLID97DOFMatrices::Ce_SOLID97];
    auto &Me = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::Me_SOLID97];
    auto &AAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::AAe_SOLID97];
    auto &JAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JAe_SOLID97];
    auto &JVe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JVe_SOLID97];
    auto &BAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::BAe_SOLID97];
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int l = 0; l < 3; l++)
                {
                    BAe(qp + k * nq, dof_i + l * n_dofs_var) += curl_phi_vec[dof_i * 3 + l][qp](k);
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] / mu_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));
                        AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * curl_phi_vec[dof_j * 3 + l][qp]);
                        Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += JxW[qp] * (phi_vec[dof_i * 3 + k][qp] * phi_vec[dof_j * 3 + l][qp]);
                    }
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
        materials_object->fill_rs_matrix(this->elem_data.material_num);
        auto &rs_matrix = materials_object->get_rs_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                for (unsigned int k = 0; k < 3; k++)
                {
                    JVe(qp + k * nq, dof_i) += dphi[dof_i][qp](k);
                    JAe(qp + k * nq, dof_i + k * n_dofs_var) += phi[dof_i][qp];
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
        //         Fe(dof_i + 3 * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
    // The penalty values, for solution boundary trace and flux.

    // Here use the formula for quadrature
    std::vector<boundary_id_type> side_boundary_ids;
    FEBase *fe_face = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[2].get());

    // The following loops over the sides of the element.
    // If the element has no neighbor on a side then that
    // side MUST live on a boundary of the domain.
    Real penalty_base;
    char cont = 0;
    auto boundary_info = &this->get_parent_object()->get_mesh().get_boundary_info();

    for (auto s : elem_data.elem->side_index_range())
        if (elem_data.elem->neighbor_ptr(s) == nullptr)
        {
            boundary_info->boundary_ids(elem_data.elem, s, side_boundary_ids);
            if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), Globals::FLUX_PARALLEL_ID_BC) != side_boundary_ids.end())
            {
                if (cont == 0)
                {
                    penalty_base = 1e7 * Ke.max() + 1e7 * Ce.max();
                    cont++;
                }
                // The value of the shape functions at the quadrature
                // points.
                const std::vector<std::vector<Real>> &phi_face =
                    fe_face->get_phi();

                // The value of the shape function derivatives at the
                // quadrature points.
                const std::vector<std::vector<RealGradient>> &dphi_face =
                    fe_face->get_dphi();

                // The Jacobian * Quadrature Weight at the quadrature
                // points on the face.
                const std::vector<Real> &JxW_face = fe_face->get_JxW();

                // The XYZ locations (in physical space) of the
                // quadrature points on the face.  This is where
                // we will interpolate the boundary value function.
                const std::vector<Point> &qface_point = fe_face->get_xyz();
                const std::vector<std::vector<Point>> &face_tangents = fe_face->get_tangents();

                // Compute the shape function values on the element
                // face.
                fe_face->reinit(elem_data.elem, s);

                // Loop over the face quadrature points for integration.
                for (unsigned int qp = 0; qp < fe_face->n_quadrature_points(); qp++)
                {
                    // Matrix contribution of the L2 projection.
                    // Note that the basis function values are
                    // integrated against test function values while
                    // basis fluxes are integrated against test function
                    // fluxes.
                    for (unsigned int i = 0; i < n_dofs_var; i++)
                        for (unsigned int j = 0; j < n_dofs_var; j++)
                            for (unsigned int l = 0; l < 3; l++)
                                for (unsigned int k = 0; k < 3; k++)
                                    for (unsigned int h = 0; h < face_tangents[0].size(); h++)
                                        Ke(i + k * n_dofs_var, j + l * n_dofs_var) += JxW_face[qp] *
                                                                                      (penalty_base * phi_face[i][qp] *
                                                                                       phi_face[j][qp] * face_tangents[qp][h](k) *
                                                                                       face_tangents[qp][h](l));
                }
            }
            if (std::find(side_boundary_ids.begin(), side_boundary_ids.end(), Globals::FLUX_NORMAL_ID_BC) != side_boundary_ids.end())
            {
                // The value of the shape functions at the quadrature
                // points.
                const std::vector<std::vector<Real>> &phi_face =
                    fe_face->get_phi();

                // The value of the shape function derivatives at the
                // quadrature points.
                const std::vector<std::vector<RealGradient>> &dphi_face =
                    fe_face->get_dphi();

                // The Jacobian * Quadrature Weight at the quadrature
                // points on the face.
                const std::vector<Real> &JxW_face = fe_face->get_JxW();

                // The XYZ locations (in physical space) of the
                // quadrature points on the face.  This is where
                // we will interpolate the boundary value function.
                const std::vector<Point> &qface_point = fe_face->get_xyz();
                const std::vector<Point> &face_normals = fe_face->get_normals();

                // Compute the shape function values on the element
                // face.
                fe_face->reinit(elem_data.elem, s);

                // Loop over the face quadrature points for integration.
                for (unsigned int qp = 0; qp < fe_face->n_quadrature_points(); qp++)
                {
                    // Matrix contribution of the L2 projection.
                    // Note that the basis function values are
                    // integrated against test function values while
                    // basis fluxes are integrated against test function
                    // fluxes.
                    // Compute the shape function values on the element
                    // face.
                    fe_face->reinit(elem_data.elem, s);

                    // Loop over the face quadrature points for integration.
                    for (unsigned int qp = 0; qp < fe_face->n_quadrature_points(); qp++)
                    {
                        // Matrix contribution of the L2 projection.
                        // Note that the basis function values are
                        // integrated against test function values while
                        // basis fluxes are integrated against test function
                        // fluxes.
                        for (unsigned int i = 0; i < n_dofs_var; i++)
                            for (unsigned int j = 0; j < n_dofs_var; j++)
                                for (unsigned int l = 0; l < 3; l++)
                                    for (unsigned int k = 0; k < 3; k++)
                                        Ke(i + k * n_dofs_var, j + l * n_dofs_var) += JxW_face[qp] *
                                                                                      (penalty_base * phi_face[i][qp] *
                                                                                       phi_face[j][qp] * face_normals[qp](k) *
                                                                                       face_normals[qp](l));
                    }
                }
            }
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
        // Pay attention here if you add other variables
        dof_map_element_loads->dof_indices(elem_data.elem, dof_indices_loads);
        element_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i + k * n_dofs_var) *= loads_vec[k];

        // Apply nodal CSG loads
        unsigned int amps_var_num = nodal_loads_system->variable_number(this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::AMPS].global_label);
        dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads, amps_var_num);
        nodal_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int dof_i = 0; dof_i < 3 * n_dofs_var; dof_i++)
            if (loads_vec[dof_i] != 0)
                Fe(dof_i) = loads_vec[dof_i];

        for (dof_id_type &dof_value : dof_indices_loads)
            nodal_loads_vector->set(dof_value, 0);

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        unsigned int amps_var_num = nodal_loads_system->variable_number(this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::AMPS].global_label);

        // Apply nodal AMPS loads
        dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads, amps_var_num);
        nodal_loads_vector->get(dof_indices_loads, loads_vec);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i + 3 * n_dofs_var) = loads_vec[dof_i + amps_var_num * n_dofs_var];

        for (dof_id_type &dof_value : dof_indices_loads)
            nodal_loads_vector->set(dof_value, 0);

        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    return;
}

void ElementBuilder::build_solid97_solution_vectors()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    DenseVector<Number> Bn;
    DenseVector<Number> Hn;
    DenseVector<Number> Be(3);
    DenseVector<Number> He(3);
    DenseVector<Number> Je(3);

    auto &Me = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::Me_SOLID97];
    auto &AAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::AAe_SOLID97];
    auto &BAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::BAe_SOLID97];

    std::vector<dof_id_type> dof_indices;

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = materials_object->get_mu_matrix();

    ApecMLinearSystem *main_sys = dynamic_cast<ApecMLinearSystem *>(systems[Globals::DOF_SOLUTION_SYSTEM]);
    System *nodal_solution_system = systems[Globals::NODAL_SOLUTION];
    System *element_solution_system = systems[Globals::ELEMENT_SOLUTION];

    this->fill_solution_vector_support(main_sys->current_local_solution.get(), Globals::DOF_SOLUTION_SYSTEM);

    unsigned int n_dofs_var = elem_data.elem->n_nodes();
    unsigned int n_qp = BAe.m() / 3;

    // Get DOF Solutions in the nodes
    DenseVector<Number> An;
    for (unsigned int k = 0; k < 3; k++)
    {
        An.append(solution_vector_support.sol_vectors[k]);
    }

    // Compute B in nodes
    DenseVector<Number> Bn_x;
    AAe.vector_mult(Bn_x, An);
    Me.cholesky_solve(Bn_x, Bn);

    // Compute H in nodes
    Hn.resize(Bn.size());
    for (unsigned int k = 0; k < 3; k++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Hn(dof_i + k * n_dofs_var) += Bn(dof_i + k * n_dofs_var) / mu_matrix[k];

    // Compute B in elements
    DenseVector<Number> Bq;
    BAe.vector_mult(Bq, An);
    for (unsigned int k = 0; k < 3; k++)
    {
        for (unsigned int qp = 0; qp < n_qp; qp++)
            Be(k) += Bq(qp + k * n_qp);
        Be(k) = Be(k) / n_qp;
    }

    // Compute H in elements
    for (unsigned int k = 0; k < 3; k++)
        He(k) = Be(k) / mu_matrix[k];

    // Compute J in elements
    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        dof_maps[Globals::ELEMENT_LOADS]->dof_indices(elem_data.elem, dof_indices);
        NumericVector<Number> *element_loads_vector = systems[Globals::ELEMENT_LOADS]->request_vector(2);
        element_loads_vector->get(dof_indices, Je.get_values());
        // DenseVector<Number> Jq;
        // BAe.vector_mult(Jq, Hn);
        // for (unsigned int k = 0; k < 3; k++)
        // {
        //     for (unsigned int qp = 0; qp < n_qp; qp++)
        //         Je(k) += Jq(qp + k * n_qp);
        //     Je(k) = Je(k) / n_qp;
        // }

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (this->context<TimeStep>()->get_transient_effect())
        {
            auto &JAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JAe_SOLID97];
            auto &JVe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JVe_SOLID97];

            materials_object->fill_rs_matrix(this->elem_data.material_num);
            auto &rs_matrix = materials_object->get_rs_matrix();
            unsigned int volt_var_num = main_sys->variable_number(this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::VOLT].global_label);

            fill_solution_vector_support(main_sys->older_local_solution, 0); // solution_derivative

            // Get derivative of A DOF Solutions (dA/dt, dV/dt)
            DenseVector<Number> An_dt;
            for (unsigned int k = 0; k < 3; k++)
            {
                An_dt.append(solution_vector_support.sol_vectors[k]);
            }
            auto &Vn_dt = solution_vector_support.sol_vectors[volt_var_num];

            // Compute dA/dt in quadrature points
            DenseVector<Number> Jq_A;
            JAe.vector_mult(Jq_A, An_dt);

            // Compute dV/dt in quadrature points
            DenseVector<Number> Jq_V;
            JVe.vector_mult(Jq_V, Vn_dt);

            // sum the two contributions
            Jq_A += Jq_V;

            // average on quadrature points and scale by conductivity
            for (unsigned int k = 0; k < 3; k++)
            {
                for (unsigned int qp = 0; qp < n_qp; qp++)
                    Je(k) -= Jq_A(qp + k * n_qp);
                Je(k) = Je(k) / rs_matrix[k] / n_qp;
            }
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

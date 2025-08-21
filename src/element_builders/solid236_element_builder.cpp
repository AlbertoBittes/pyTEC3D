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
#include "element_types_factory.h"

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236()
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

    const std::vector<std::vector<RealTensorValue>> &dphi_vec = fe_vec->get_dphi();

    fe->reinit(elem_data.elem);
    fe_vec->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_A = fe_vec->n_shape_functions();
    unsigned int n_dofs_V = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_mu_matrix(this->elem_data.material_num);
    auto &mu_matrix = materials_object->get_mu_matrix();

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][SOLID97DOFMatrices::Ke_SOLID97];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][SOLID97DOFMatrices::Ce_SOLID97];
    // auto &Me = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::Me_SOLID97];
    // auto &AAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::AAe_SOLID97];
    // auto &JAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JAe_SOLID97];
    // auto &JVe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::JVe_SOLID97];
    // auto &BAe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][SOLID97SolutionMatrices::BAe_SOLID97];
    // auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                    Ke(dof_i, dof_j) += JxW[qp] / mu_matrix[k] * (curl_phi_vec[dof_i][qp](k) * curl_phi_vec[dof_j][qp](k));

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        materials_object->fill_rs_matrix(this->elem_data.material_num);

        auto &rs_matrix = materials_object->get_rs_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i, dof_j) += JxW[qp] / rs_matrix[k] * (phi_vec[dof_i][qp](k) * phi_vec[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i, dof_j + n_dofs_A) += JxW[qp] / rs_matrix[k] * (phi_vec[dof_i][qp](k) * dphi[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i + n_dofs_A, dof_j) += JxW[qp] / rs_matrix[k] * (phi_vec[dof_j][qp](k) * dphi[dof_i][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += JxW[qp] / rs_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));

        // CAV_value = JxW[qp] / rs_matrix[l] * phi[dof_i][qp] * dphi[dof_j][qp](l);
        // Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) += CAV_value;
        // Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) += CAV_value;
        // Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) += JxW[qp] * dphi[dof_i][qp](l) / rs_matrix[l] * dphi[dof_j][qp](l);
        // for (unsigned int k = 0; k < 3; k++)
        // {
        //     JVe(qp + k * n_dofs_var, dof_i) += dphi[dof_i][qp](k);
        //     JAe(qp + k * n_dofs_var, dof_i + k * n_dofs_var) += phi[dof_i][qp];
        // }

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
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236_rhs()
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

    std::vector<dof_id_type> dof_indices_loads;
    std::vector<Real> loads_vec;

    System *main_sys = systems[Globals::DOF_SOLUTION_SYSTEM];
    System *nodal_solution_system = systems[Globals::NODAL_LOADS];
    System *element_solution_system = systems[Globals::ELEMENT_LOADS];

    auto dof_map_nodal_loads = dof_maps[Globals::NODAL_LOADS];
    auto dof_map_element_loads = dof_maps[Globals::ELEMENT_LOADS];

    NumericVector<Number> *nodal_loads_vector = nodal_solution_system->request_vector(2);
    NumericVector<Number> *element_loads_vector = element_solution_system->request_vector(2);

    switch (elem_data.element_type->get_physical_region())
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        FEVectorBase *fe_vec = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[2].get());

        // The velocity shape function gradients at interior
        // quadrature points.
        const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_vec->get_phi();

        const std::vector<Real> &JxW = fe_vec->get_JxW();
        fe_vec->reinit(elem_data.elem);

        this->prepare_element_matrix_building();

        // Get the scalar shape functions and vector shape functions
        unsigned int n_dofs_A = fe_vec->n_shape_functions();
        unsigned int nq = fe_vec->n_quadrature_points();

        dof_map_element_loads->dof_indices(elem_data.elem, dof_indices_loads);
        element_loads_vector->get(dof_indices_loads, loads_vec);
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int k = 0; k < 3; k++)
                    Fe(dof_i) += JxW[qp] * phi_vec[dof_i][qp](k) * loads_vec[k];

        // // Apply nodal CSG loads
        // dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads);
        // nodal_loads_vector->get(dof_indices_loads, loads_vec);

        // for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //     if (loads_vec[dof_i] != 0)
        //         Fe(dof_i) = loads_vec[dof_i];

        // for (dof_id_type &dof_value : dof_indices_loads)
        //     nodal_loads_vector->set(dof_value, 0);

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        unsigned int volt_var_num = main_sys->variable_number(this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::VOLT].global_label);
        unsigned int amp_var_num = main_sys->variable_number(this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::AMPS].global_label);

        // Apply nodal AMPS loads
        dof_map_nodal_loads->dof_indices(elem_data.elem, dof_indices_loads);
        nodal_loads_vector->get(dof_indices_loads, loads_vec);

        unsigned int az_var_num = main_sys->variable_number(this->elem_data.element_type->get_variable_information()[Globals::GlobalVariableName::Az].global_label);

        unsigned int n_dofs_A = elem_data.elem->n_dofs(az_var_num, Globals::DOF_SOLUTION_SYSTEM);
        unsigned int n_dofs_V = elem_data.elem->n_dofs(volt_var_num, Globals::DOF_SOLUTION_SYSTEM);

        for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
            Fe(dof_i + n_dofs_A) = loads_vec[dof_i + n_dofs_A];

        for (dof_id_type &dof_value : dof_indices_loads)
            nodal_loads_vector->set(dof_value, 0);

        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    return;
}

void ElementBuilder::build_solid236_solution_vectors()
{
    return;
}

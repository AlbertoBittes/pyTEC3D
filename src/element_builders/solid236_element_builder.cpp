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
#include "apec_m_system.h"
#include "element_types_factory.h"
#include "variables_map.h"

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
    FEBase *fe = dynamic_cast<FEBase *>(elem_data.element_type->fe_array[3].get());
    FEVectorBase *fe_edge = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[0].get());
    FEVectorBase *fe_vec_lagrange = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[2].get());
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();
    // const std::vector<Real> &JxW_edge = fe_edge->get_JxW();
    // const std::vector<Real> &JxW_vec_lagrange = fe_vec_lagrange->get_JxW();
    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients for the velocity
    // variables evaluated at the quadrature points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // The velocity shape function gradients at interior
    // quadrature points.
    const std::vector<std::vector<RealVectorValue>> &phi_edge = fe_edge->get_phi();

    const std::vector<std::vector<RealGradient>> &curl_phi_edge = fe_edge->get_curl_phi();

    const std::vector<std::vector<RealGradient>> &phi_vec_lagrange = fe_vec_lagrange->get_phi();
    // const std::vector<std::vector<RealTensorValue>> &dphi_vec = fe_vec->get_dphi();

    // const auto &qp_lagr = fe->get_xyz();
    // const auto &qp_edge = fe_edge->get_xyz();
    // const auto &qp_vec_lagrange = fe_vec_lagrange->get_xyz();

    fe->reinit(elem_data.elem);
    fe_edge->reinit(elem_data.elem);
    fe_vec_lagrange->reinit(elem_data.elem);

    this->prepare_element_matrix_building();

    // Get the scalar shape functions and vector shape functions
    unsigned int n_dofs_A = fe_edge->n_shape_functions();
    unsigned int n_dofs_V = fe->n_shape_functions();
    unsigned int nq = fe->n_quadrature_points();
    unsigned int nq_edge = fe_edge->n_quadrature_points();
    unsigned int nq_lagrange = fe_vec_lagrange->n_quadrature_points();

    auto materials_object = this->context<MaterialsObject>();

    materials_object->fill_nu_matrix(this->elem_data.material_num);
    auto &nu_matrix = materials_object->get_nu_matrix();

    auto &Ke = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ke];
    auto &Ce = this->matrices[MatricesSet::DOF_SOLUTION_MATRICES][DOFMatrices::Ce];
    auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::Me];
    auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::AAe];
    auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];
    auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JVe];
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::BAe];
    // auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    for (unsigned int qp = 0; qp < nq_edge; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                for (unsigned int k = 0; k < 3; k++)
                    Ke(dof_i, dof_j) += JxW[qp] * nu_matrix[k] * (curl_phi_edge[dof_i][qp](k) * curl_phi_edge[dof_j][qp](k));

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    AAe(dof_i + k * n_dofs_V, dof_j) += JxW[qp] * (phi_vec_lagrange[dof_i * 3 + k][qp] * curl_phi_edge[dof_j][qp]);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int l = 0; l < 3; l++)
                        Me(dof_i + k * n_dofs_V, dof_j + l * n_dofs_V) += JxW[qp] * (phi_vec_lagrange[dof_i * 3 + k][qp] * phi_vec_lagrange[dof_j * 3 + l][qp]);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                BAe(qp + k * nq, dof_j) += curl_phi_edge[dof_j][qp](k);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                JAe(qp + k * nq, dof_j) += phi_edge[dof_j][qp](k);

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                JVe(qp + k * nq, dof_j) += dphi[dof_j][qp](k);

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);

        auto &sigma_matrix = materials_object->get_sigma_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i, dof_j) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * phi_edge[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * dphi[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i + n_dofs_A, dof_j) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * phi_edge[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));

        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i + 3 * n_dofs_var) += JxW[qp] * phi[dof_i][qp]; //* js_loads[elem->id()][k];
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i, dof_j) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * phi_edge[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_A; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ke(dof_i, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (phi_edge[dof_i][qp](k) * dphi[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_A; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ce(dof_i + n_dofs_A, dof_j) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * phi_edge[dof_j][qp](k));

        for (unsigned int qp = 0; qp < nq; qp++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    for (unsigned int k = 0; k < 3; k++)
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += JxW[qp] * sigma_matrix[k] * (dphi[dof_i][qp](k) * dphi[dof_j][qp](k));

        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        auto real_constants_object = this->context<RealConstantsObject>();
        if (!real_constants_object || !real_constants_object->has_constant(this->elem_data.real_constant_num))
            libmesh_error_msg("Missing real constants for STRANDED_COIL elements.");

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real vc = rc[3];
        const Real tx = rc[5];
        const Real ty = rc[6];
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL elements.");

        if (n_dofs >= n_dofs_A + 2 * n_dofs_V)
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

            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real kvv = JxW[qp] * ((1.0 / (r_coil * vc)) * phi[dof_i][qp] * phi[dof_j][qp]);
                        Ce(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ce(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A + n_dofs_V) -= kvv;
                    }
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        auto real_constants_object = this->context<RealConstantsObject>();
        if (!real_constants_object || !real_constants_object->has_constant(this->elem_data.real_constant_num))
            libmesh_error_msg("Missing real constants for STRANDED_COIL_SOLENOIDAL elements.");

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real vc = rc[3];
        const Real tx = rc[5];
        const Real ty = rc[6];
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL_SOLENOIDAL elements.");

        if (n_dofs >= n_dofs_A + 2 * n_dofs_V)
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

            for (unsigned int qp = 0; qp < nq; qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_V; dof_j++)
                    {
                        const Real kvv = JxW[qp] * ((1.0 / (r_coil * vc)) * phi[dof_i][qp] * phi[dof_j][qp]);
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A) += sym * kvv;
                        Ke(dof_i + n_dofs_A, dof_j + n_dofs_A + n_dofs_V) -= kvv;
                        Ke(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A) -= kvv;
                        Ke(dof_i + n_dofs_A + n_dofs_V, dof_j + n_dofs_A + n_dofs_V) += kvv;
                    }
        }
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

    System *main_sys = &this->get_parent_object()->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    System *nodal_loads_sys = &this->get_parent_object()->get_system(ApecMSystemBase::NODAL_LOADS);
    System *element_loads_sys = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_LOADS);

    auto dof_map_main_system = &main_sys->get_dof_map();
    auto dof_map_nodal_loads = &nodal_loads_sys->get_dof_map();
    auto dof_map_element_loads = &element_loads_sys->get_dof_map();

    NumericVector<Number> *nodal_loads_vector = nodal_loads_sys->request_vector(2);
    NumericVector<Number> *element_loads_vector = element_loads_sys->request_vector(2);

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        FEVectorBase *fe_edge = dynamic_cast<FEVectorBase *>(elem_data.element_type->fe_array[1].get());

        // The velocity shape function gradients at interior
        // quadrature points.
        const std::vector<std::vector<RealVectorValue>> &phi_vec = fe_edge->get_phi();

        const std::vector<Real> &JxW = fe_edge->get_JxW();
        fe_edge->reinit(elem_data.elem);

        this->prepare_element_matrix_building();

        // Get the scalar shape functions and vector shape functions
        unsigned int n_dofs_A = fe_edge->n_shape_functions();
        unsigned int nq = fe_edge->n_quadrature_points();

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
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    case PhysicalRegionEnum::STRANDED_COIL:
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        std::string volt_var_name;
        if (PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL == static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]) ||
            PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL == static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
            volt_var_name = "V";
        else
            volt_var_name = "VOLT";
        if (nodal_loads_sys->has_variable("AMPS") && main_sys->has_variable(volt_var_name) && main_sys->has_variable("Az"))
        {
            unsigned int amp_var_num = nodal_loads_sys->variable_number("AMPS");
            unsigned int volt_var_num = main_sys->variable_number(volt_var_name);
            unsigned int az_var_num = main_sys->variable_number("Az");

            unsigned int n_dofs_A = elem_data.elem->n_dofs(ApecMSystemBase::DOF_SOLUTION_SYSTEM, az_var_num);
            unsigned int n_dofs_V = elem_data.elem->n_dofs(ApecMSystemBase::DOF_SOLUTION_SYSTEM, volt_var_num);

            unsigned int proc_id = main_sys->comm().rank();
            dof_indices_loads.clear();
            for (const auto &node : elem_data.elem->node_ref_range())
                if (node.processor_id() == proc_id)
                    if (node.n_dofs(ApecMSystemBase::NODAL_LOADS, amp_var_num))
                        dof_indices_loads.push_back(node.dof_number(ApecMSystemBase::NODAL_LOADS, amp_var_num, 0));

            nodal_loads_vector->get(dof_indices_loads, loads_vec);

            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
            {
                const dof_id_type dof_id = dof_indices_loads[dof_i];
                if (!nodal_loads_dof_cache.insert(dof_id).second)
                    continue;
                Fe(dof_i + n_dofs_A) = loads_vec[dof_i];
            }
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }
    return;
}

void ElementBuilder::build_solid236_solution_vectors()
{
    std::vector<dof_id_type> dof_indices;

    auto materials_object = this->context<MaterialsObject>();
    materials_object->fill_nu_matrix(this->elem_data.material_num);
    auto &nu_matrix = materials_object->get_nu_matrix();

    auto real_constants_object = this->context<RealConstantsObject>();

    ApecMLinearSystem *main_sys = dynamic_cast<ApecMLinearSystem *>(&this->get_parent_object()->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    System *nodal_solution_system = &this->get_parent_object()->get_system(ApecMSystemBase::NODAL_SOLUTION);
    System *element_solution_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_SOLUTION);

    // Retrieve Az DOF solution
    this->fill_solution_vector_support(main_sys->current_local_solution.get(), ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    auto &Az = solution_vector_support.sol_vectors[0];

    // Nodal solution
    if (this->compute_nodal_solution_vectors)
    {
        auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::Me];
        auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::AAe];
        this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS].resize(NodalSolutionVectors::NUM_NODAL_SOLUTION_VECTORS);
        auto &Bn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::B];
        auto &Hn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::H];

        unsigned int n_dofs_A = AAe.n();
        unsigned int n_dofs_V = AAe.m() / 3;

        // Compute B in nodes
        DenseVector<Number> Bn_x;
        AAe.vector_mult(Bn_x, Az);
        Me.cholesky_solve(Bn_x, Bn);

        // Compute H in nodes
        Hn.resize(Bn.size());
        for (unsigned int k = 0; k < 3; k++)
            for (unsigned int dof_i = 0; dof_i < n_dofs_V; dof_i++)
                Hn(dof_i + k * n_dofs_V) += Bn(dof_i + k * n_dofs_V) * nu_matrix[k];
    }

    // Elemental solution
    this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS].resize(ElementSolutionVectors::NUM_ELEMENT_SOLUTION_VECTORS);
    auto &Be = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::Be];
    auto &He = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::He];
    auto &Je = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::Je];
    Be.resize(3);
    He.resize(3);
    Je.resize(3);

    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::BAe];
    unsigned int n_qp = BAe.m() / 3;
    // Compute B in elements
    DenseVector<Number> Bq;
    BAe.vector_mult(Bq, solution_vector_support.sol_vectors[0]);
    for (unsigned int k = 0; k < 3; k++)
    {
        for (unsigned int qp = 0; qp < n_qp; qp++)
            Be(k) += Bq(qp + k * n_qp);
        Be(k) = Be(k) / n_qp;
    }

    // Compute H in elements
    for (unsigned int k = 0; k < 3; k++)
        He(k) = Be(k) * nu_matrix[k];

    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]);
    // Compute J in elements
    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        System *element_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_LOADS);
        DofMap *dof_map_element_loads = &element_loads_system->get_dof_map();
        dof_map_element_loads->dof_indices(elem_data.elem, dof_indices);
        NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);
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
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];
        auto &JVe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JVe];

        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        DenseVector<Number> Vn;

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
            Vn = solution_vector_support.sol_vectors[1];
        // Get derivative of A DOF Solutions (dA/dt, dV/dt)
        fill_solution_vector_support(main_sys->request_vector("solution_time_derivative"), 0); // solution_derivative

        auto &A_dt = solution_vector_support.sol_vectors[0];

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR || physical_region == PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR)
            Vn = solution_vector_support.sol_vectors[1];

        // Compute dA/dt in quadrature points
        DenseVector<Number> Jq_A;
        JAe.vector_mult(Jq_A, A_dt);

        // Compute dV/dt in quadrature points
        DenseVector<Number> Jq_V;
        JVe.vector_mult(Jq_V, Vn);

        // sum the two contributions
        Jq_A += Jq_V;

        // average on quadrature points and scale by conductivity
        for (unsigned int k = 0; k < 3; k++)
        {
            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(k) -= Jq_A(qp + k * n_qp);
            Je(k) = Je(k) * sigma_matrix[k] / n_qp;
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        Real Vn = 0.0;
        Real En = 0.0;

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
        {
            Vn = solution_vector_support.sol_vectors[1][0];
            En = solution_vector_support.sol_vectors[2][0];
        }

        fill_solution_vector_support(main_sys->request_vector("solution_time_derivative"), ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
        {
            Vn = solution_vector_support.sol_vectors[1][0];
            En = solution_vector_support.sol_vectors[2][0];
        }

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];

        const VectorValue t_vec(rc[5], rc[6], rc[7]);
        const Real r_coil = rc[8];

        const Real coeff = nc / (sc * r_coil) * (Vn - En);
        if (sc == 0.0 || r_coil == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R) for STRANDED_COIL_SOLENOIDAL elements.");

        for (unsigned int k = 0; k < 3; k++)
            Je(k) = coeff * t_vec(k);
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    if (this->compute_nodal_solution_vectors)
    {
        // Set B and H in nodes
        DofMap *dof_map_nodal_solution = &nodal_solution_system->get_dof_map();
        dof_map_nodal_solution->dof_indices(elem_data.elem, dof_indices);
        auto &Bn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::B];
        auto &Hn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::H];
        Bn.append(Hn);
        for (unsigned int i = 0; i < dof_indices.size(); i++)
        {
            dof_id_type dof_index = dof_indices[i];
            if ((nodal_solution_system->solution->first_local_index() <= dof_index) &&
                (dof_index < nodal_solution_system->solution->last_local_index()))
                nodal_solution_system->solution->set(dof_index, Bn(i));
        }
    }

    DofMap *dof_map_element_solution = &element_solution_system->get_dof_map();
    dof_map_element_solution->dof_indices(elem_data.elem, dof_indices);
    Je.append(Be);
    Je.append(He);
    for (unsigned int i = 0; i < dof_indices.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((element_solution_system->solution->first_local_index() <= dof_index) &&
            (dof_index < element_solution_system->solution->last_local_index()))
            element_solution_system->solution->set(dof_index, Je(i));
    }
    return;
}

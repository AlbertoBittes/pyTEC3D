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
    auto &Fe = this->vectors[VectorsSet::RHS_VECTORS][0];

    unsigned int axisymm = elem_data.element_type->get_keyopts()[3];

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
            BAe(qp + nq, dof_i) -= r_sign * (dphi[dof_i][qp](0) + axisymm * phi[dof_i][qp] / r);
            JAe(qp, dof_i) += phi[dof_i][qp];
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
    Ke.scale(nu_matrix[0]);
    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]);

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
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            {
                // matrices[4](0, i) -= scalar_shape_functions->JxW[qp] * scalar_shape_functions->phi[i][qp] / scalar_shape_functions->phi[0].size();
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp] * pow(r, axisymm));
            }
        }
        Ce.scale(sigma_matrix[0]);
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

    std::vector<dof_id_type> dof_indices_loads;
    std::vector<Real> loads_vec;

    System *nodal_loads_system = &this->get_parent_object()->get_system(Globals::NODAL_LOADS);
    System *element_loads_system = &this->get_parent_object()->get_system(Globals::ELEMENT_LOADS);

    DofMap *dof_map_nodal_loads = &nodal_loads_system->get_dof_map();
    DofMap *dof_map_element_loads = &element_loads_system->get_dof_map();

    NumericVector<Number> *nodal_loads_vector = nodal_loads_system->request_vector(2);
    NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        if (element_loads_system->has_variable("Js_z"))
        {
            dof_map_element_loads->dof_indices(elem_data.elem, dof_indices_loads);
            element_loads_vector->get(dof_indices_loads, loads_vec);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                Fe(dof_i) *= loads_vec[0];
        }
        // if (nodal_loads_system->has_variable("CSG_z"))
        // {
        //     // Apply nodal CSG loads
        //     nodal_loads_vector->get(dof_indices_loads, loads_vec);

        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         if (loads_vec[dof_i] != 0)
        //             Fe(dof_i) = loads_vec[dof_i];

        //     for (dof_id_type &dof_value : dof_indices_loads)
        //         nodal_loads_vector->set(dof_value, 0);
        // }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        // Apply nodal CSG loads+
        // if (nodal_loads_system->has_variable("CSG_z"))
        // {
        //     unsigned int csg_var_num = nodal_loads_system->variable_number("CSG_z");

        //     unsigned int proc_id = nodal_loads_system->comm().rank();
        //     dof_indices_loads.clear();
        //     for (const auto &node : elem_data.elem->node_ref_range())
        //         if (node.processor_id() == proc_id)
        //             dof_indices_loads.push_back(node.dof_number(Globals::NODAL_LOADS, csg_var_num, 0));

        //     nodal_loads_vector->get(dof_indices_loads, loads_vec);

        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         if (loads_vec[dof_i] != 0)
        //             Fe(dof_i) = loads_vec[dof_i];

        //     for (dof_id_type &dof_value : dof_indices_loads)
        //         nodal_loads_vector->set(dof_value, 0);
        // }

        break;
        // unsigned int volt_var_num = this->elem_data.element_type->get_variable_information()[GlobalVariableName::VOLT].system_var_num;
        // unsigned int amp_var_num = this->elem_data.element_type->get_variable_information()[GlobalVariableName::AMPS].system_var_num;

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

void ElementBuilder::build_plane233_solution_vectors()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    std::vector<dof_id_type> dof_indices;

    ApecMLinearSystem *main_sys = dynamic_cast<ApecMLinearSystem *>(&this->get_parent_object()->get_system(Globals::DOF_SOLUTION_SYSTEM));
    System *nodal_solution_system = &this->get_parent_object()->get_system(Globals::NODAL_SOLUTION);
    System *element_solution_system = &this->get_parent_object()->get_system(Globals::ELEMENT_SOLUTION);

    auto materials_object = this->context<MaterialsObject>();
    materials_object->fill_nu_matrix(this->elem_data.material_num);
    auto &nu_matrix = materials_object->get_nu_matrix();

    // Retrieve A DOF solution
    this->fill_solution_vector_support(main_sys->current_local_solution.get(), Globals::DOF_SOLUTION_SYSTEM);
    auto &An = solution_vector_support.sol_vectors[0];

    bool axisymm = this->elem_data.element_type->get_keyopts()[3];

    Point centroid = this->elem_data.elem->centroid();

    // Nodal solution
    if (this->compute_nodal_solution_vectors)
    {
        this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS].resize(NodalSolutionVectors::NUM_NODAL_SOLUTION_VECTORS);
        auto &Bn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::B];
        auto &Hn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::H];
        auto &psi = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::PSI];
        auto &Me = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::Me];
        auto &AAe = this->matrices[MatricesSet::NODAL_SOLUTION_MATRICES][NodalSolutionMatrices::AAe];
        unsigned int n_dofs_var = elem_data.elem->n_nodes();

        // Compute B in nodes
        DenseVector<Number> Bn_x;
        AAe.vector_mult(Bn_x, An);
        // Bn.resize(Me.m());
        Me.cholesky_solve(Bn_x, Bn);

        // Compute H in nodes
        Hn.resize(Bn.size());
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Hn(dof_i) += Bn(dof_i) * nu_matrix[0];

        if (axisymm)
        {
            psi.resize(An.size());
            // unsigned int a_var_num = main_sys->variable_number("A_z");
            // unsigned int psi_var_num = nodal_solution_system->variable_number("PSI");
            for (unsigned int i = 0; i < elem_data.elem->n_nodes(); i++)
            {
                auto &node = elem_data.elem->node_ref(i);
                // unsigned int dof_A = node.dof_number(Globals::DOF_SOLUTION_SYSTEM, a_var_num, 0);
                // unsigned int dof_PSI = node.dof_number(Globals::NODAL_SOLUTION, psi_var_num, 0);
                // if ((main_sys->current_local_solution->first_local_index() <= dof_A) &&
                //     (dof_A < main_sys->current_local_solution->last_local_index()) &&
                //     (nodal_solution_system->solution->first_local_index() <= dof_PSI) &&
                //     (dof_PSI < nodal_solution_system->solution->last_local_index()))
                // {
                psi(i) = An(i) * node(0); // r component
                // }
            }
        }
    }

    // Element solution
    this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS].resize(ElementSolutionVectors::NUM_ELEMENT_SOLUTION_VECTORS);
    auto &Be = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::Be];
    auto &He = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::He];
    auto &Je = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::Je];
    auto &PSIe = this->vectors[VectorsSet::ELEMENT_SOLUTION_VECTORS][ElementSolutionVectors::PSIe];
    Be.resize(2);
    He.resize(2);
    Je.resize(1);
    PSIe.resize(1);

    // Compute B in elements
    auto &BAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::BAe];
    unsigned int n_qp = BAe.m() / 2;
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
        He(k) = Be(k) * nu_matrix[0];

    if (axisymm)
    {
        auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];

        DenseVector<Number> Aq(n_qp);
        JAe.vector_mult(Aq, An);

        for (unsigned int qp = 0; qp < n_qp; qp++)
            PSIe(0) += Aq(qp);
        PSIe(0) = PSIe(0) / n_qp * centroid(0); // r component
    }
    // Compute J in elements
    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
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
        auto element_loads_system = &this->get_parent_object()->get_system(Globals::ELEMENT_LOADS);
        auto element_loads_dof_map = &element_loads_system->get_dof_map();
        element_loads_dof_map->dof_indices(elem_data.elem, dof_indices);
        NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);
        element_loads_vector->get(dof_indices, Je.get_values());

        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        if (this->context<TimeStep>()->get_transient_effect())
        {
            auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];

            materials_object->fill_sigma_matrix(this->elem_data.material_num);
            auto &sigma_matrix = materials_object->get_sigma_matrix();

            fill_solution_vector_support(main_sys->request_vector("solution_time_derivative"), Globals::DOF_SOLUTION_SYSTEM); // solution_derivative

            // Get derivative of A DOF Solutions (dA/dt, dV/dt)
            auto &An_dt = solution_vector_support.sol_vectors[0];

            // Compute dA/dt in quadrature points
            DenseVector<Number> Jq_A(n_qp);
            JAe.vector_mult(Jq_A, An_dt);

            // average on quadrature points and scale by conductivity
            for (unsigned int qp = 0; qp < n_qp; qp++)
                Je(0) -= Jq_A(qp);
            Je(0) = Je(0) * sigma_matrix[0] / n_qp;
        }
        break;
    }
    default:
        libmesh_error_msg("Physical region not implemented");
    }

    // Set B and H in nodes
    if (this->compute_nodal_solution_vectors)
    {
        auto &Bn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::B];
        auto &Hn = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::H];
        auto &psi = this->vectors[VectorsSet::NODAL_SOLUTION_VECTORS][NodalSolutionVectors::PSI];
        DofMap *dof_map_nodal_solution = &nodal_solution_system->get_dof_map();
        dof_map_nodal_solution->dof_indices(elem_data.elem, dof_indices);
        Bn.append(Hn);
        if (axisymm)
            Bn.append(psi);
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
    if (axisymm)
        Je.append(PSIe);
    for (unsigned int i = 0; i < dof_indices.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((element_solution_system->solution->first_local_index() <= dof_index) &&
            (dof_index < element_solution_system->solution->last_local_index()))
            element_solution_system->solution->set(dof_index, Je(i));
    }

    // // Compute PSI
    // if (axisymm)
    // {
    //     unsigned int a_var_num = main_sys->variable_number("A_z");
    //     unsigned int psi_var_num = nodal_solution_system->variable_number("PSI");
    //     for (auto &node : elem_data.elem->node_ref_range())
    //     {
    //         unsigned int dof_A = node.dof_number(Globals::DOF_SOLUTION_SYSTEM, a_var_num, 0);
    //         unsigned int dof_PSI = node.dof_number(Globals::NODAL_SOLUTION, psi_var_num, 0);
    //         if ((main_sys->current_local_solution->first_local_index() <= dof_A) &&
    //             (dof_A < main_sys->current_local_solution->last_local_index()) &&
    //             (nodal_solution_system->solution->first_local_index() <= dof_PSI) &&
    //             (dof_PSI < nodal_solution_system->solution->last_local_index()))
    //         {
    //             Real psi_val = main_sys->current_local_solution->el(dof_A) * node(0); // r component
    //             nodal_solution_system->solution->set(dof_PSI, psi_val);
    //         }
    //     }
    // }
}
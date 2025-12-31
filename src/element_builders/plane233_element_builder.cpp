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
#include "apec_m_system.h"
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
    auto real_constants_object = this->context<RealConstantsObject>();

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
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                Ke(dof_i, dof_j) +=
                    JxW[qp] * (dphi[dof_i][qp] * dphi[dof_j][qp] + axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r_sq) * r_mult;
    }
    Ke.scale(nu_matrix[0]);

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        {
            BAe(qp, dof_i) += r_sign * dphi[dof_i][qp](1);
            BAe(qp + nq, dof_i) -= r_sign * (dphi[dof_i][qp](0) + axisymm * phi[dof_i][qp] / r);
        }
    }

    for (unsigned int qp = 0; qp < nq; qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            JAe(qp, dof_i) += phi[dof_i][qp];

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        Real r_mult = pow(r, axisymm);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                AAe(dof_i, dof_j) += JxW[qp] * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](1)) * r_mult;
    }

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        Real r_mult = pow(r, axisymm);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                AAe(dof_i + n_dofs_var, dof_j) -=
                    JxW[qp] * (r_sign * phi[dof_i][qp] * dphi[dof_j][qp](0) - axisymm * phi[dof_i][qp] * phi[dof_j][qp] / r) * r_mult;
    }

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        Real r_mult = pow(r, axisymm);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                Me(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
    }

    for (unsigned int qp = 0; qp < nq; qp++)
    {
        Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
        Real r_mult = pow(r, axisymm);

        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                Me(dof_i + n_dofs_var, dof_j + n_dofs_var) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
    }

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
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();
        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp] * pow(r, axisymm));
        }
        Ce.scale(sigma_matrix[0]);
        break;
        // Here when you will implement loads with F
        // for (unsigned int qp = 0; qp < nq; qp++)
        //     for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
        //         Fe(dof_i) += JxW[qp] * phi[dof_i][qp]; // * js_loads[elem->id()][2];
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        Real thk = 1.0;
        if (real_constants_object && real_constants_object->has_constant(this->elem_data.real_constant_num))
            thk = real_constants_object->get_constant(this->elem_data.real_constant_num)[0];
        if (thk == 0.0)
            libmesh_error_msg("THK is zero for SOLID_CONDUCTOR elements.");

        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(r, axisymm);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
        }
        Ce.scale(sigma_matrix[0]);

        if (n_dofs >= 2 * n_dofs_var)
        {
            const Real inv_thk = 1.0 / thk;
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        Ce(dof_i, dof_j + n_dofs_var) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) * inv_thk;
                        Ce(dof_i + n_dofs_var, dof_j) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
                        Ce(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) * inv_thk;
                    }
            }
        }
        break;
    }
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        Real thk = 1.0;
        if (real_constants_object && real_constants_object->has_constant(this->elem_data.real_constant_num))
            thk = real_constants_object->get_constant(this->elem_data.real_constant_num)[0];
        if (thk == 0.0)
            libmesh_error_msg("THK is zero for SOLID_CONDUCTOR_SOLENOIDAL elements.");

        for (unsigned int qp = 0; qp < nq; qp++)
        {
            Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
            Real r_mult = pow(r, axisymm);

            for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    Ce(dof_i, dof_j) += JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
        }
        Ce.scale(sigma_matrix[0]);

        if (n_dofs >= 2 * n_dofs_var)
        {
            const Real inv_thk = 1.0 / thk;
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        Ke(dof_i, dof_j + n_dofs_var) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) * inv_thk;
                        Ce(dof_i + n_dofs_var, dof_j) -=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
                        Ke(dof_i + n_dofs_var, dof_j + n_dofs_var) +=
                            JxW[qp] * sigma_matrix[0] * (phi[dof_i][qp] * phi[dof_j][qp]) * inv_thk;
                    }
            }
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL:
    {
        if (!real_constants_object || !real_constants_object->has_constant(this->elem_data.real_constant_num))
            libmesh_error_msg("Missing real constants for STRANDED_COIL elements.");

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real vc = rc[3];
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL elements.");

        const Real kav_scale = -(nc / (sc * r_coil)) * tz;
        const Real kvv_scale = 1.0 / (r_coil * vc);

        if (n_dofs >= 3 * n_dofs_var)
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
                        const Real kav = kav_scale * nn;
                        const Real kvv = kvv_scale * nn;

                        // Time-integrated A-VOLT-EMF formulation (KEYOPT(2)=2).
                        Ce(dof_i, dof_j + n_dofs_var) += kav;
                        Ce(dof_i, dof_j + 2 * n_dofs_var) -= kav;
                        Ce(dof_i + n_dofs_var, dof_j) += sym * kav;
                        Ce(dof_i + n_dofs_var, dof_j + n_dofs_var) += sym * kvv;
                        Ce(dof_i + 2 * n_dofs_var, dof_j) -= kav;
                        Ce(dof_i + 2 * n_dofs_var, dof_j + 2 * n_dofs_var) -= kvv;
                    }
            }
        }
        break;
    }
    case PhysicalRegionEnum::STRANDED_COIL_SOLENOIDAL:
    {
        if (!real_constants_object || !real_constants_object->has_constant(this->elem_data.real_constant_num))
            libmesh_error_msg("Missing real constants for STRANDED_COIL_SOLENOIDAL elements.");

        const auto &rc = real_constants_object->get_constant(this->elem_data.real_constant_num);
        const Real sc = rc[1];
        const Real nc = rc[2];
        const Real vc = rc[3];
        const Real tz = rc[7];
        const Real r_coil = rc[8];
        const Real sym = rc[9];

        if (sc == 0.0 || r_coil == 0.0 || vc == 0.0)
            libmesh_error_msg("Invalid real constants (SC, R, VC) for STRANDED_COIL_SOLENOIDAL elements.");

        const Real kav_scale = -(nc / (sc * r_coil)) * tz;
        const Real kvv_scale = 1.0 / (r_coil * vc);

        if (n_dofs >= 3 * n_dofs_var)
        {
            for (unsigned int qp = 0; qp < nq; qp++)
            {
                Real r = sqrt(xyz[qp](2) * xyz[qp](2) + xyz[qp](0) * xyz[qp](0));
                Real r_mult = pow(r, axisymm);

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        const Real nn = JxW[qp] * (phi[dof_i][qp] * phi[dof_j][qp]) * r_mult;
                        const Real kav = kav_scale * nn;
                        const Real kvv = kvv_scale * nn;

                        // Non-time-integrated A-VOLT-EMF formulation (KEYOPT(2)=0/1).
                        Ke(dof_i, dof_j + n_dofs_var) += kav;
                        Ke(dof_i, dof_j + 2 * n_dofs_var) -= kav;
                        Ke(dof_i + n_dofs_var, dof_j + n_dofs_var) += sym * kvv;
                        Ke(dof_i + n_dofs_var, dof_j + 2 * n_dofs_var) -= kvv;
                        Ke(dof_i + 2 * n_dofs_var, dof_j + n_dofs_var) -= kvv;
                        Ke(dof_i + 2 * n_dofs_var, dof_j + 2 * n_dofs_var) += kvv;

                        Ce(dof_i + n_dofs_var, dof_j) += sym * kav;
                        Ce(dof_i + 2 * n_dofs_var, dof_j) -= kav;
                    }
            }
        }
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

    std::vector<dof_id_type> dof_indices_loads;
    std::vector<Real> loads_vec;

    System *main_sys = &this->get_parent_object()->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM);
    System *nodal_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::NODAL_LOADS);
    System *element_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_LOADS);

    DofMap *dof_map_nodal_loads = &nodal_loads_system->get_dof_map();
    DofMap *dof_map_element_loads = &element_loads_system->get_dof_map();

    NumericVector<Number> *nodal_loads_vector = nodal_loads_system->request_vector(2);
    NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);

    switch (static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]))
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
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
        if (nodal_loads_system->has_variable("AMPS") && main_sys->has_variable(volt_var_name) && main_sys->has_variable("A_z"))
        {
            unsigned int amp_var_num = nodal_loads_system->variable_number("AMPS");
            unsigned int v_var_num = main_sys->variable_number(volt_var_name);
            unsigned int az_var_num = main_sys->variable_number("A_z");

            unsigned int n_dofs_A = elem_data.elem->n_dofs(ApecMSystemBase::DOF_SOLUTION_SYSTEM, az_var_num);
            unsigned int n_dofs_V = elem_data.elem->n_dofs(ApecMSystemBase::DOF_SOLUTION_SYSTEM, v_var_num);

            if (n_dofs_V != n_dofs_A)
                libmesh_error_msg("Number of DOFs for V and A_z do not match in stranded coil solenoidal elements");

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

void ElementBuilder::build_plane233_solution_vectors()
{
    // ApecMSolutionSystem *solution_system = dynamic_cast<ApecMSolutionSystem *>(this->get_parent_object());
    std::vector<dof_id_type> dof_indices;

    ApecMLinearSystem *main_sys = dynamic_cast<ApecMLinearSystem *>(&this->get_parent_object()->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    System *nodal_solution_system = &this->get_parent_object()->get_system(ApecMSystemBase::NODAL_SOLUTION);
    System *element_solution_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_SOLUTION);

    auto materials_object = this->context<MaterialsObject>();
    materials_object->fill_nu_matrix(this->elem_data.material_num);
    auto &nu_matrix = materials_object->get_nu_matrix();

    auto real_constants_object = this->context<RealConstantsObject>();

    // Retrieve A DOF solution
    this->fill_solution_vector_support(main_sys->current_local_solution.get(), ApecMSystemBase::DOF_SOLUTION_SYSTEM);
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
                // unsigned int dof_A = node.dof_number(ApecMSystemBase::DOF_SOLUTION_SYSTEM, a_var_num, 0);
                // unsigned int dof_PSI = node.dof_number(ApecMSystemBase::NODAL_SOLUTION, psi_var_num, 0);
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
    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(elem_data.element_type->get_keyopts()[1]);
    switch (physical_region)
    {
    case PhysicalRegionEnum::NON_MAGNETIC:
    case PhysicalRegionEnum::SOFT_MAGNETIC:
    case PhysicalRegionEnum::HARD_MAGNETIC:
        break;
    case PhysicalRegionEnum::STRANDED_CONDUCTOR:
    {
        // auto &JHe = this->matrices[MatricesSet::NODAL_ELEMENT_SOLUTION_MATRICES][PLANE13SolutionMatrices::JHe_PLANE13];

        // DenseVector<Number> Jq(n_qp);
        // JHe.vector_mult(Jq, Hn);

        // for (unsigned int qp = 0; qp < n_qp; qp++)
        //     Je(0) += Jq(qp);
        // Je(0) = Je(0) / n_qp;
        auto element_loads_system = &this->get_parent_object()->get_system(ApecMSystemBase::ELEMENT_LOADS);
        auto element_loads_dof_map = &element_loads_system->get_dof_map();
        element_loads_dof_map->dof_indices(elem_data.elem, dof_indices);
        NumericVector<Number> *element_loads_vector = element_loads_system->request_vector(2);
        element_loads_vector->get(dof_indices, Je.get_values());

        break;
    }
    case PhysicalRegionEnum::SHORT_CIRCUIT_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR:
    case PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL:
    {
        Real Vn = 0.0;
        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR_SOLENOIDAL)
            Vn = solution_vector_support.sol_vectors[1][0];

        fill_solution_vector_support(main_sys->request_vector("solution_time_derivative"), ApecMSystemBase::DOF_SOLUTION_SYSTEM); // solution_derivative

        if (physical_region == PhysicalRegionEnum::SOLID_CONDUCTOR)
            Vn = solution_vector_support.sol_vectors[1][0];

        auto &JAe = this->matrices[MatricesSet::ELEMENT_SOLUTION_MATRICES][ElementSolutionMatrices::JAe];

        materials_object->fill_sigma_matrix(this->elem_data.material_num);
        auto &sigma_matrix = materials_object->get_sigma_matrix();

        // Get derivative of A DOF Solutions (dA/dt, dV/dt)
        auto &An_dt = solution_vector_support.sol_vectors[0];

        // Compute dA/dt in quadrature points
        DenseVector<Number> Jq_A(n_qp);
        JAe.vector_mult(Jq_A, An_dt);

        // average on quadrature points and scale by conductivity
        for (unsigned int qp = 0; qp < n_qp; qp++)
            Je(0) -= Jq_A(qp);
        Je(0) = Je(0) * sigma_matrix[0] / n_qp;

        Real thk = 1.0;
        if (real_constants_object && real_constants_object->has_constant(this->elem_data.real_constant_num))
            thk = real_constants_object->get_constant(this->elem_data.real_constant_num)[0];
        if (thk == 0.0)
            libmesh_error_msg("THK is zero for SOLID_CONDUCTOR_SOLENOIDAL elements.");
        Je(0) += Vn / thk;
        unsigned int axisymm = elem_data.element_type->get_keyopts()[3];
        if (axisymm)
            Je(0) /= 2 * libMesh::pi * centroid(0); // r component
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
    //         unsigned int dof_A = node.dof_number(ApecMSystemBase::DOF_SOLUTION_SYSTEM, a_var_num, 0);
    //         unsigned int dof_PSI = node.dof_number(ApecMSystemBase::NODAL_SOLUTION, psi_var_num, 0);
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

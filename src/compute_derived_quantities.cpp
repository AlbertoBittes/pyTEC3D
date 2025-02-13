// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/explicit_system.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/utility.h"
#include "libmesh/dense_matrix.h"
#include "../include/compute_derived_quantities.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/transient_system.h"
#include "../include/system_builder.h"
#include "libmesh/string_to_enum.h"
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <set>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include "libmesh/enum_solver_package.h"
#include "vtkDataArray.h"
#include "../include/json.hpp"
#include "../include/assembler.h"
#include "../include/element_builder.h" // Include the header file for ElementBuilder
#include <map>                          // Include for std::map
#include <vtkCellData.h>
#include <cmath>
#include <stdexcept>
using namespace libMesh;

// Constructor implementation
DerivedQuantitiesCalculator::DerivedQuantitiesCalculator(std::shared_ptr<AssemblerLinear> assembler, EquationSystems &equation_systems)
    : assembler(assembler), equation_systems(equation_systems) {}

void DerivedQuantitiesCalculator::compute_derived_quantities(std::string sys_name)
{
    const MeshBase &mesh = equation_systems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    ExplicitSystem &derived_quantities_system = equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystem");
    ExplicitSystem &derived_quantities_system_elem = equation_systems.get_system<ExplicitSystem>("DerivedQuantitiesSystemElem");
    TransientLinearImplicitSystem &system = equation_systems.get_system<TransientLinearImplicitSystem>("PlasmaArcLinearSystem");

    const unsigned int a_var = system.variable_number("A");
    unsigned int v_var;
    if (dim == 3)
        v_var = system.variable_number("V");

    const DofMap &dof_map = system.get_dof_map();

    unsigned int b_n_var = derived_quantities_system.variable_number("B");
    unsigned int b_elem_var = derived_quantities_system_elem.variable_number("B_elem");
    unsigned int j_elem_var = derived_quantities_system_elem.variable_number("J_elem");
    unsigned int f_elem_var = derived_quantities_system_elem.variable_number("F_elem");
    const DofMap &derived_quantities_system_dof_map = derived_quantities_system.get_dof_map();
    const DofMap &derived_quantities_system_elem_dof_map = derived_quantities_system_elem.get_dof_map();

    std::shared_ptr<TimeStepper> timeStepper = equation_systems.parameters.get<std::shared_ptr<TimeStepper>>("time_stepper");
    Real dt = timeStepper->get_time_substep_dt();
    bool trans_effects = timeStepper->get_time_steps()[timeStepper->get_current_time_step()].get_transient_effect();
    // Storage for the stress dof indices on each element
    std::vector<dof_id_type> B_indices_var;
    std::vector<dof_id_type> A_indices_var;
    std::vector<dof_id_type> V_indices_var;
    std::vector<dof_id_type> J_indices_var;

    MaterialBuilder mat_builder(*(equation_systems.parameters.get<MaterialsMap *>("materials_properties")), equation_systems.parameters.get<vtkCellData *>("cell_data"));

    DenseMatrix<Number> AAe;
    DenseMatrix<Number> Me;
    DenseMatrix<Number> JAe;
    DenseMatrix<Number> JVe;
    DenseVector<Number> Ae;
    DenseVector<Number> Ae_dt;
    DenseVector<Number> Ae_temp;
    DenseVector<Number> Ve;
    DenseVector<Number> Ve_dt;
    DenseVector<Number> Be;
    DenseVector<Number> Fe;
    DenseVector<Number> Be_q;
    DenseVector<Number> J_elem;
    DenseVector<Number> J_elem_A;
    DenseVector<Number> J_elem_V;
    DenseVector<Number> B_elem;
    DenseVector<Number> F_elem;
    char n_J_dofs;
    if (dim == 2)
    {
        n_J_dofs = 1;
        J_elem.resize(n_J_dofs);
        // B_elem.resize(2);
        // F_elem.resize(2);
    }
    else
    {
        n_J_dofs = 3;
        J_elem.resize(n_J_dofs);
        B_elem.resize(n_J_dofs);
        F_elem.resize(n_J_dofs);
    }
    dof_id_type elem_id;
    char n_A_dofs;
    char n_V_dofs;
    char n_B_dofs;
    bool conductive;
    std::set<subdomain_id_type> *conductive_subdomain = equation_systems.parameters.get<std::set<subdomain_id_type> *>("conductive_subdomain");
    unsigned int local_elem_idx = 0;
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        elem_id = elem->id();

        if (conductive_subdomain->find(elem->subdomain_id()) != conductive_subdomain->end())
            conductive = true;
        else
            conductive = false;

        if (!trans_effects)
            conductive = false;

        mat_builder.compute_material_properties(local_elem_idx, conductive);

        dof_map.dof_indices(elem, A_indices_var, a_var);
        n_A_dofs = std::size(A_indices_var);
        Ae.resize(n_A_dofs);
        Ae_dt.resize(n_A_dofs);

        derived_quantities_system_dof_map.dof_indices(elem, B_indices_var, b_n_var);
        n_B_dofs = std::size(B_indices_var);
        Be_q.resize(n_B_dofs);
        Be.resize(n_B_dofs);

        derived_quantities_system_elem_dof_map.dof_indices(elem, J_indices_var, j_elem_var);
        Ae_temp.resize(n_A_dofs);
        B_elem.zero();
        F_elem.zero();
        if (conductive)
        {
            J_elem.zero();
            if (dim == 3)
            {
                dof_map.dof_indices(elem, V_indices_var, v_var);
                n_V_dofs = std::size(V_indices_var);
                Ve.resize(n_V_dofs);
                Ve_dt.resize(n_V_dofs);
                J_elem_A.zero();
                J_elem_V.zero();
            }
        }

        for (std::size_t j = 0; j < n_A_dofs; j++)
        {
            Ae(j) = system.current_solution(A_indices_var[j]);
            if (conductive)
            {
                Ae_temp(j) = system.old_solution(A_indices_var[j]);
                Ae_dt(j) = (Ae(j) - system.old_solution(A_indices_var[j])) / dt;

                if (dim == 3)
                {
                    if (j < n_V_dofs)
                    {
                        Ve_dt(j) = (system.current_solution(V_indices_var[j]) - system.old_solution(V_indices_var[j])) / dt;
                    }
                }
            }
        }

        AAe = assembler->get_aae_matrix_derived(elem_id);

        AAe.vector_mult(Be_q, Ae);

        // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
        // Note that this gives a discontinuous stress plot on element boundaries, which is
        // appropriate. We then also get the von Mises stress from the projected stress tensor.
        Me = assembler->get_me_matrix_derived(elem_id);

        Me.cholesky_solve(Be_q, Be);

        unsigned int j_qp_size;
        if (conductive)
        {
            JAe = assembler->get_jae_matrix_derived(elem_id);
            if (dim == 2)
            {
                j_qp_size = JAe.m();
                J_elem_A.resize(j_qp_size);
                JAe.vector_mult(J_elem_A, Ae);
                for (char k = 0; k < n_J_dofs; k++)
                    for (char qp = 0; qp < j_qp_size; qp++)
                        J_elem(k) -= J_elem_A(qp + k) / mat_builder.rs_matrix(k) / (j_qp_size / n_J_dofs);
            }
            else
            {
                JVe = assembler->get_jve_matrix_derived(elem_id);
                j_qp_size = JAe.m();
                J_elem_A.resize(j_qp_size);
                J_elem_V.resize(j_qp_size);
                JAe.vector_mult(J_elem_A, Ae_dt);
                JVe.vector_mult(J_elem_V, Ve_dt);
                for (char k = 0; k < n_J_dofs; k++)
                    for (char qp = 0; qp < j_qp_size / n_J_dofs; qp++)
                        J_elem(k) += (-J_elem_A(qp * 3 + k) - J_elem_V(qp * 3 + k)) / mat_builder.rs_matrix(k) / (j_qp_size / n_J_dofs);
            }
        }

        dof_id_type dof_index;
        // Loop over B_indices_var_values
        for (unsigned int i = 0; i < n_B_dofs; i++)
        {
            dof_index = B_indices_var[i];
            if ((derived_quantities_system.solution->first_local_index() <= dof_index) &&
                (dof_index < derived_quantities_system.solution->last_local_index()))
                derived_quantities_system.solution->set(dof_index, Be(i));
        }
        if (conductive)
            for (unsigned int i = 0; i < n_J_dofs; i++)
            {
                dof_index = J_indices_var[i];
                if ((derived_quantities_system.solution->first_local_index() <= dof_index) &&
                    (dof_index < derived_quantities_system.solution->last_local_index()))
                    derived_quantities_system_elem.solution->set(dof_index, J_elem(i));
            }
        // Should call close and update when we set vector entries directly
        derived_quantities_system.solution->close();
        derived_quantities_system.update();
        derived_quantities_system_elem.solution->close();
        derived_quantities_system_elem.update();
        local_elem_idx++;
    }
}
// const MeshBase &mesh = es.get_mesh();
// const unsigned int dim = mesh.mesh_dimension();
// LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>(0);

// const unsigned int a_var = system.variable_number("AZ");

// const DofMap &dof_map = system.get_dof_map();

// FEType fe_type = dof_map.variable_type(a_var);

// std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

// QGauss qrule(dim, fe_type.default_quadrature_order());

// fe->attach_quadrature_rule(&qrule);

// // Here we define some references to cell-specific data that
// // will be used to assemble the linear system.
// // We begin with the element Jacobian * quadrature weight at each
// // integration point.
// const std::vector<Real> &JxW = fe->get_JxW();

// // The element shape functions evaluated at the quadrature points.
// const std::vector<std::vector<Real>> &phi = fe->get_phi();

// // The element shape function gradients evaluated at the quadrature
// // points.
// const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();
// const std::vector<Real> &dxidy = fe->get_dxidy();
// const std::vector<Real> &detady = fe->get_detady();
// const std::vector<Real> &dxidx = fe->get_dxidx();
// const std::vector<Real> &detadx = fe->get_detadx();

// const std::vector<std::vector<Real>> &dphidxi = fe->get_dphidxi();
// const std::vector<std::vector<Real>> &dphideta = fe->get_dphideta();

// // Also, get a reference to the ExplicitSystem
// magnetic_field_system.add_variable("BX", FIRST, L2_LAGRANGE);
// magnetic_field_system.add_variable("BY", FIRST, L2_LAGRANGE);
// magnetic_field_system.init();
// unsigned int bx_n_var = magnetic_field_system.variable_number("BX");
// unsigned int by_n_var = magnetic_field_system.variable_number("BY");
// // magnetic_field_system.add_variable("BY", FIRST, L2_LAGRANGE);
// const DofMap &B_dof_map = magnetic_field_system.get_dof_map();

// // Storage for the stress dof indices on each element
// std::vector<dof_id_type> Bx_indices_var;
// std::vector<dof_id_type> By_indices_var;
// std::vector<dof_id_type> A_indices_var;
// DenseMatrix<Number> Kxe;
// DenseVector<Number> Ae;
// DenseVector<Number> Bxe_q;
// DenseMatrix<Number> Kye;
// DenseVector<Number> Bye_q;
// DenseVector<Number> Bxe;
// DenseVector<Number> Bye;
// for (const auto &elem : mesh.active_local_element_ptr_range())
// {
//     dof_map.dof_indices(elem, A_indices_var);

//     const unsigned int n_var_dofs = A_indices_var.size();
//     const unsigned int n_dofs = cast_int<unsigned int>(n_var_dofs);

//     fe->reinit(elem);

//     Kxe.resize(n_dofs, n_dofs);
//     Ae.resize(n_dofs);
//     Bxe_q.resize(n_dofs);
//     Kye.resize(n_dofs, n_dofs);
//     Bye_q.resize(n_dofs);
//     Bxe.resize(n_dofs);
//     Bye.resize(n_dofs);
//     for (std::size_t j = 0; j < n_dofs; j++)
//         Ae(j) = system.current_solution(A_indices_var[j]);
//     // Build a vector containing dxixy, detady
//     std::vector<std::vector<Real>> jacob_derivatives_bx = {dxidy, detady};
//     std::vector<std::vector<Real>> jacob_derivatives_by = {dxidx, detadx};
//     std::vector<std::vector<std::vector<Real>>> phi_derivatives = {dphidxi, dphideta};
//     for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//     {
//         // Row is variable u1, u2, or u3, column is x, y, or z
//         for (std::size_t i = 0; i < n_dofs; i++)
//             for (std::size_t j = 0; j < n_dofs; j++)
//                 for (std::size_t l = 0; l < dim; l++)
//                 {
//                     Kxe(i, j) += JxW[qp] * phi[i][qp] * (jacob_derivatives_bx[l][qp] * phi_derivatives[l][j][qp]);
//                     Kye(i, j) -= JxW[qp] * phi[i][qp] * (jacob_derivatives_by[l][qp] * phi_derivatives[l][j][qp]);
//                 }
//     }

//     Kxe.vector_mult(Bxe_q, Ae);
//     Kye.vector_mult(Bye_q, Ae);
//     B_dof_map.dof_indices(elem, Bx_indices_var, bx_n_var);
//     B_dof_map.dof_indices(elem, By_indices_var, by_n_var);

//     // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
//     // Note that this gives a discontinuous stress plot on element boundaries, which is
//     // appropriate. We then also get the von Mises stress from the projected stress tensor.
//     DenseMatrix<Real> Me_bx(n_dofs, n_dofs);
//     DenseMatrix<Real> Me_by(n_dofs, n_dofs);
//     for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//     {
//         for (unsigned int i = 0; i < n_dofs; i++)
//             for (unsigned int j = 0; j < n_dofs; j++)
//             {
//                 Me_bx(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
//                 Me_by(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
//             }
//     }

//     Me_bx.cholesky_solve(Bxe_q, Bxe);
//     Me_by.cholesky_solve(Bye_q, Bye);

//     dof_id_type dof_index;
//     // Loop over B_indices_var_values
//     for (unsigned int i = 0; i < n_dofs; i++)
//     {
//         dof_index = Bx_indices_var[i];
//         Real val;
//         if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
//             (dof_index < magnetic_field_system.solution->last_local_index()))
//             magnetic_field_system.solution->set(dof_index, Bxe(i));
//         dof_index = By_indices_var[i];
//         if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
//             (dof_index < magnetic_field_system.solution->last_local_index()))
//             magnetic_field_system.solution->set(dof_index, Bye(i));
//     }
// }
// // Should call close and update when we set vector entries directly
// magnetic_field_system.solution->close();
// magnetic_field_system.update();
// }

// Post-process the solution to compute stresses
// void compute_magnetic_field_static_3D(EquationSystems &es)
// {
//     const MeshBase &mesh = es.get_mesh();
//     const unsigned int dim = mesh.mesh_dimension();
//     ExplicitSystem &magnetic_field_system = es.add_system<ExplicitSystem>("MagneticFieldSystem");
//     LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>(0);

//     const unsigned int a_var = system.variable_number("A");

//     const DofMap &dof_map = system.get_dof_map();

//     FEType fe_type = dof_map.variable_type(a_var);

//     std::unique_ptr<FEVectorBase> fe(FEVectorBase::build(dim, fe_type));

//     QGauss qrule(dim, fe_type.default_quadrature_order());

//     fe->attach_quadrature_rule(&qrule);

//     const std::vector<Real> &JxW = fe->get_JxW();
//     const std::vector<std::vector<RealGradient>> &phi = fe->get_phi();
//     const std::vector<std::vector<RealGradient>> &curl_phi = fe->get_curl_phi();

//     // Also, get a reference to the ExplicitSystem
//     magnetic_field_system.add_variable("B", FIRST, L2_LAGRANGE_VEC);
//     magnetic_field_system.init();
//     unsigned int b_n_var = magnetic_field_system.variable_number("B");
//     const DofMap &B_dof_map = magnetic_field_system.get_dof_map();

//     // Storage for the stress dof indices on each element
//     std::vector<dof_id_type> B_indices_var;
//     std::vector<dof_id_type> A_indices_var;
//     DenseMatrix<Number> Ke;
//     DenseVector<Number> Ae;
//     DenseVector<Number> Be;
//     DenseVector<Number> Be_q;
//     for (const auto &elem : mesh.active_local_element_ptr_range())
//     {
//         dof_map.dof_indices(elem, A_indices_var);

//         const unsigned int n_var_dofs = A_indices_var.size();
//         const unsigned int n_dofs = cast_int<unsigned int>(n_var_dofs);

//         fe->reinit(elem);

//         Ke.resize(n_dofs, n_dofs);
//         Ae.resize(n_dofs);
//         Be_q.resize(n_dofs);

//         for (std::size_t j = 0; j < n_dofs; j++)
//             Ae(j) = system.current_solution(A_indices_var[j]);

//         for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//         {
//             // Row is variable u1, u2, or u3, column is x, y, or z
//             for (std::size_t i = 0; i < n_dofs; i++)
//                 for (std::size_t j = 0; j < n_dofs; j++)
//                     Ke(i, j) += JxW[qp] * (phi[i][qp] * curl_phi[j][qp]);
//         }

//         Ke.vector_mult(Be_q, Ae);
//         B_dof_map.dof_indices(elem, B_indices_var);

//         // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
//         // Note that this gives a discontinuous stress plot on element boundaries, which is
//         // appropriate. We then also get the von Mises stress from the projected stress tensor.
//         DenseMatrix<Real> Me(n_dofs, n_dofs);
//         for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//         {
//             for (unsigned int i = 0; i < n_dofs; i++)
//                 for (unsigned int j = 0; j < n_dofs; j++)
//                 {
//                     Me(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
//                 }
//         }

//         Me.cholesky_solve(Be_q, Be);

//         dof_id_type dof_index;
//         // Loop over B_indices_var_values
//         for (unsigned int i = 0; i < n_dofs; i++)
//         {
//             dof_index = B_indices_var[i];
//             if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
//                 (dof_index < magnetic_field_system.solution->last_local_index()))
//                 magnetic_field_system.solution->set(dof_index, Be(i));
//         }
//     }
//     // Should call close and update when we set vector entries directly
//     magnetic_field_system.solution->close();
//     magnetic_field_system.update();
// }

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
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/utility.h"
#include "../include/loads_bc_object.h"
#include "../include/global_resources.h"

using namespace libMesh;

// Post-process the solution to compute stresses
// void compute_magnetic_field_static(EquationSystems &es)
// {
//     const MeshBase &mesh = es.get_mesh();
//     const unsigned int dim = mesh.mesh_dimension();
//     ExplicitSystem &magnetic_field_system = es.add_system<ExplicitSystem>("MagneticFieldSystem");
//     LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>("Elasticity");

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
//     const DofMap &stress_dof_map = magnetic_field_system.get_dof_map();
//     magnetic_field_system.add_variable("B", FIRST, L2_LAGRANGE_VEC);

//     // Storage for the stress dof indices on each element
//     std::vector<std::vector<dof_id_type>> B_indices_var(system.n_vars());
//     std::vector<dof_id_type> stress_dof_indices_var;
//     std::vector<dof_id_type> vonmises_dof_indices_var;

//     for (const auto &elem : mesh.active_local_element_ptr_range())
//     {
//         for (unsigned int var = 0; var < 3; var++)
//             dof_map.dof_indices(elem, dof_indices_var[var], displacement_vars[var]);

//         const unsigned int n_var_dofs = dof_indices_var[0].size();

//         fe->reinit(elem);

//         std::vector<TensorValue<Number>> stress_tensor_qp(qrule.n_points());
//         for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//         {
//             // Row is variable u1, u2, or u3, column is x, y, or z
//             TensorValue<Number> grad_u;
//             for (unsigned int var_i = 0; var_i < 3; var_i++)
//                 for (unsigned int var_j = 0; var_j < 3; var_j++)
//                     for (unsigned int j = 0; j < n_var_dofs; j++)
//                         grad_u(var_i, var_j) += dphi[j][qp](var_j) * system.current_solution(dof_indices_var[var_i][j]);

//             for (unsigned int var_i = 0; var_i < 3; var_i++)
//                 for (unsigned int var_j = 0; var_j < 3; var_j++)
//                     for (unsigned int k = 0; k < 3; k++)
//                         for (unsigned int l = 0; l < 3; l++)
//                             stress_tensor_qp[qp](var_i, var_j) += elasticity_tensor(var_i, var_j, k, l) * grad_u(k, l);
//         }

//         stress_dof_map.dof_indices(elem, vonmises_dof_indices_var, vonMises_var);
//         std::vector<TensorValue<Number>> elem_sigma_vec(vonmises_dof_indices_var.size());

//         // Below we project each component of the stress tensor onto a L2_LAGRANGE discretization.
//         // Note that this gives a discontinuous stress plot on element boundaries, which is
//         // appropriate. We then also get the von Mises stress from the projected stress tensor.
//         unsigned int stress_var_index = 0;
//         for (unsigned int var_i = 0; var_i < 3; var_i++)
//             for (unsigned int var_j = var_i; var_j < 3; var_j++)
//             {
//                 stress_dof_map.dof_indices(elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

//                 const unsigned int n_proj_dofs = stress_dof_indices_var.size();

//                 DenseMatrix<Real> Me(n_proj_dofs, n_proj_dofs);
//                 for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//                 {
//                     for (unsigned int i = 0; i < n_proj_dofs; i++)
//                         for (unsigned int j = 0; j < n_proj_dofs; j++)
//                         {
//                             Me(i, j) += JxW[qp] * (phi[i][qp] * phi[j][qp]);
//                         }
//                 }

//                 DenseVector<Number> Fe(n_proj_dofs);
//                 for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
//                     for (unsigned int i = 0; i < n_proj_dofs; i++)
//                     {
//                         Fe(i) += JxW[qp] * stress_tensor_qp[qp](var_i, var_j) * phi[i][qp];
//                     }

//                 DenseVector<Number> projected_data;
//                 Me.cholesky_solve(Fe, projected_data);

//                 for (unsigned int index = 0; index < n_proj_dofs; index++)
//                 {
//                     dof_id_type dof_index = stress_dof_indices_var[index];
//                     if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
//                         (dof_index < magnetic_field_system.solution->last_local_index()))
//                         magnetic_field_system.solution->set(dof_index, projected_data(index));

//                     elem_sigma_vec[index](var_i, var_j) = projected_data(index);
//                 }

//                 stress_var_index++;
//             }

//         for (std::size_t index = 0; index < elem_sigma_vec.size(); index++)
//         {
//             elem_sigma_vec[index](1, 0) = elem_sigma_vec[index](0, 1);
//             elem_sigma_vec[index](2, 0) = elem_sigma_vec[index](0, 2);
//             elem_sigma_vec[index](2, 1) = elem_sigma_vec[index](1, 2);

//             // Get the von Mises stress from the projected stress tensor
//             Number vonMises_value = std::sqrt(0.5 * (Utility::pow<2>(elem_sigma_vec[index](0, 0) - elem_sigma_vec[index](1, 1)) +
//                                                      Utility::pow<2>(elem_sigma_vec[index](1, 1) - elem_sigma_vec[index](2, 2)) +
//                                                      Utility::pow<2>(elem_sigma_vec[index](2, 2) - elem_sigma_vec[index](0, 0)) +
//                                                      6. * (Utility::pow<2>(elem_sigma_vec[index](0, 1)) +
//                                                            Utility::pow<2>(elem_sigma_vec[index](1, 2)) +
//                                                            Utility::pow<2>(elem_sigma_vec[index](2, 0)))));

//             dof_id_type dof_index = vonmises_dof_indices_var[index];

//             if ((magnetic_field_system.solution->first_local_index() <= dof_index) &&
//                 (dof_index < magnetic_field_system.solution->last_local_index()))
//                 magnetic_field_system.solution->set(dof_index, vonMises_value);
//         }
//     }

//     // Should call close and update when we set vector entries directly
//     magnetic_field_system.solution->close();
//     magnetic_field_system.update();
// }

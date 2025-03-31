#include "../include/element_builder.h"
#include "libmesh/equation_systems.h"
#include "libmesh/string_to_enum.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "../include/session.h" // Ensure the Session class is included
#include <cmath>
#include "../include/element_types.h" // Include the header defining ElementType

// Constructor
ElementBuilder::ElementBuilder(std::shared_ptr<ElementDataManager> element_data_manager)
    : element_data_manager(element_data_manager)
{
}

// Method

void ElementBuilder::build_solid97(const ElementType &element_type, bool build_matrix)
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
    ScalarShapeFunctions shape_functions = element_type.get_scalar_shape_functions();
    VectorShapeFunctions vector_shape_functions = element_type.get_vector_shape_functions();

    unsigned int n_dofs_var = shape_functions.phi.size();
    unsigned int n_dofs_A = n_dofs_var * 3;

    if (build_matrix)
    {
        if (element_type.is_conductive())
        {
            // Now start logging the element matrix computation
            Real CAV_value;

            for (unsigned int qp = 0; qp < element_type.fe->n_quadrature_points(); qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                {
                    JAe(qp, dof_i) += shape_functions.phi[dof_i][qp];

                    for (unsigned int k = 0; k < 3; k++)
                        JVe(qp + k * n_dofs_var, dof_i) += shape_functions.dphi[dof_i][qp](k);

                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                        {
                            Ce(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] / element_data_manager->rs_matrix[l] * (shape_functions.phi[dof_i][qp] * shape_functions.phi[dof_j][qp]);
                            CAV_value = shape_functions.JxW[qp] / element_data_manager->rs_matrix[l] * shape_functions.phi[dof_i][qp] * shape_functions.dphi[dof_j][qp](l);
                            Ce(dof_i + l * n_dofs_var, n_dofs_A + dof_j) += CAV_value;
                            Ce(n_dofs_A + dof_j, dof_i + l * n_dofs_var) += CAV_value;
                            Ce(n_dofs_A + dof_i, n_dofs_A + dof_j) += shape_functions.JxW[qp] * shape_functions.dphi[dof_i][qp](l) / element_data_manager->rs_matrix[l] * shape_functions.dphi[dof_j][qp](l);
                            for (unsigned int k = 0; k < 3; k++)
                            {
                                Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] / element_data_manager->mu_matrix[k] * (shape_functions.dphi[dof_i][qp](k) * shape_functions.dphi[dof_j][qp](k));
                                AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] * (vector_shape_functions.phi_vec[dof_i * 3 + k][qp] * vector_shape_functions.curl_phi_vec[dof_j * 3 + l][qp]);
                                Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] * (vector_shape_functions.phi_vec[dof_i * 3 + k][qp] * vector_shape_functions.phi_vec[dof_j * 3 + l][qp]);
                            }
                        }
                }
        }
        else
            for (unsigned int qp = 0; qp < element_type.fe->n_quadrature_points(); qp++)
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            for (unsigned int k = 0; k < 3; k++)
                            {
                                Ke(dof_i + l * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] / element_data_manager->mu_matrix[k] * (shape_functions.dphi[dof_i][qp](k) * shape_functions.dphi[dof_j][qp](k));
                                AAe(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] * (vector_shape_functions.phi_vec[dof_i * 3 + k][qp] * vector_shape_functions.curl_phi_vec[dof_j * 3 + l][qp]);
                                Me(dof_i + k * n_dofs_var, dof_j + l * n_dofs_var) += shape_functions.JxW[qp] * (vector_shape_functions.phi_vec[dof_i * 3 + k][qp] * vector_shape_functions.phi_vec[dof_j * 3 + l][qp]);
                            }
    }

    for (unsigned int qp = 0; qp < element_type.fe->n_quadrature_points(); qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var * 3; dof_i++)
            for (unsigned int k = 0; k < 3; k++)
                Fe(dof_i + k * n_dofs_var) += shape_functions.JxW[qp] * (shape_functions.phi[dof_i][qp]);
}

// Temporary placeholder for build_solid236
void ElementBuilder::build_solid236(const ElementType &element_type, bool build_matrix)
{
    libmesh_error_msg("build_solid236 is not implemented yet.");
}

// Temporary placeholder for build_plane13
void ElementBuilder::build_plane13(const ElementType &element_type, bool build_matrix)
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
    ScalarShapeFunctions shape_functions = element_type.get_scalar_shape_functions();
    VectorShapeFunctions vector_shape_functions = element_type.get_vector_shape_functions();

    unsigned int n_dofs_var = shape_functions.phi.size();

    if (build_matrix)
    {
        if (element_type.is_conductive())
        {
            // Now start logging the element matrix computation
            Real CAV_value;

            for (unsigned int qp = 0; qp < element_type.fe->n_quadrature_points(); qp++)
            {
                for (unsigned int dof_i = 0; dof_i < 2 * n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < 2 * n_dofs_var; dof_j++)
                    {
                        if (dof_j < n_dofs_var)
                            AAe(dof_i, dof_j) += shape_functions.JxW[qp] * std::pow(-1, dof_i) * (shape_functions.phi[dof_i / 2][qp] * shape_functions.dphi[dof_j][qp](1 - (dof_i % 2)));

                        if ((dof_i + dof_j) % 2 == 0)
                            Me(dof_i, dof_j) += shape_functions.JxW[qp] * (shape_functions.phi[dof_i / 2][qp] * shape_functions.phi[dof_j / 2][qp]);
                    }
                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                {
                    // JAe(0, i) -= shape_functions.JxW[qp] * shape_functions.phi[i][qp] / shape_functions.phi[0].size();
                    JAe(qp, dof_i) -= shape_functions.phi[dof_i][qp];

                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                    {
                        Ke(dof_i, dof_j) += shape_functions.JxW[qp] * (shape_functions.dphi[dof_i][qp] * shape_functions.dphi[dof_j][qp]) / element_data_manager->mu_matrix[0];
                        Ce(dof_i, dof_j) += shape_functions.JxW[qp] / element_data_manager->rs_matrix[0] * (shape_functions.phi[dof_i][qp] * shape_functions.phi[dof_j][qp]);
                    }
                }
            }
        }
        else
            for (unsigned int qp = 0; qp < element_type.fe->n_quadrature_points(); qp++)
            {
                for (unsigned int dof_i = 0; dof_i < 2 * n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < 2 * n_dofs_var; dof_j++)
                    {
                        if (dof_j < n_dofs_var)
                            AAe(dof_i, dof_j) += shape_functions.JxW[qp] * std::pow(-1, dof_i) * (shape_functions.phi[dof_i / 2][qp] * shape_functions.dphi[dof_j][qp](1 - (dof_i % 2)));

                        if ((dof_i + dof_j) % 2 == 0)
                            Me(dof_i, dof_j) += shape_functions.JxW[qp] * (shape_functions.phi[dof_i / 2][qp] * shape_functions.phi[dof_j / 2][qp]);
                    }

                for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
                    for (unsigned int dof_j = 0; dof_j < n_dofs_var; dof_j++)
                        for (unsigned int l = 0; l < 3; l++)
                            for (unsigned int k = 0; k < 3; k++)
                                Ke(dof_i, dof_j) += shape_functions.JxW[qp] * (shape_functions.dphi[dof_i][qp] * shape_functions.dphi[dof_j][qp]) / element_data_manager->mu_matrix[0];
            }
    }

    for (unsigned int qp = 0; qp < element_type.fe->n_quadrature_points(); qp++)
        for (unsigned int dof_i = 0; dof_i < n_dofs_var; dof_i++)
            Fe(dof_i) += shape_functions.JxW[qp] * shape_functions.phi[dof_i][qp];
}

void ElementBuilder::build_element(const ElementType &element_type, bool build_matrix)
{
    if (element_type.properties.element_type == "SOLID97")
    {
        build_solid97(element_type, build_matrix);
    }
    else if (element_type.properties.element_type == "SOLID236")
    {
        build_solid236(element_type, build_matrix);
    }
    else if (element_type.properties.element_type == "PLANE13")
    {
        build_plane13(element_type, build_matrix);
    }
    else
    {
        libmesh_error_msg("Element type not supported: " + element_type.properties.element_type);
    }
}

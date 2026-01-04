#include "element_types.h"
#include "element_builder.h"
// Direct constructor
SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(
          keyopts,
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_solid97_fe_structs(fe_structs, keyopts);
          })
{
}

// Chain-prolonging constructor for derived types
SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts,
                         std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV3DElementType(
          keyopts,
          [fe_collector, keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_solid97_fe_structs(fe_structs, keyopts);
          })
{
}

// Helper to add two FEDataShapeFunctions structs
void SOLID97Type::add_solid97_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                         const std::vector<unsigned int> &keyopts)
{
    // Example: add a scalar and a vector FE struct
    // fe_structs.push_back({
    //     static_cast<Order>(keyopts[2]), // order
    //     LAGRANGE,                       // family
    //     3,                              // dimension
    //     INVALID_ORDER,                  // quadrature_order
    //     3,                              // quadrature_dim
    //     false,                          // is_infinite
    //     false                           // is_vector
    // });
    fe_structs.push_back({
        static_cast<Order>(keyopts[2]), // order
        LAGRANGE_VEC,                   // family
        3,                              // dimension
        INVALID_ORDER,                  // quadrature_order (use default)
        3,                              // quadrature_dim
        false,                          // is_infinite
        true                            // is_vector (vector FE)
    });
    fe_structs.push_back({
        static_cast<Order>(keyopts[2]), // order
        LAGRANGE,                       // family
        3,                              // dimension
        INVALID_ORDER,                  // quadrature_order (use default)
        2,                              // quadrature_dim
        false,                          // is_infinite
        false                           // is_vector (vector FE)
    });
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID97Type::get_solution_matrix_idx_size_map(bool nodal_solution)
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[2]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[2]->n_quadrature_points();

    unsigned int num_matrices = nodal_solution ? ElementBuilder::kNumNodalSolutionMatrices
                                               : ElementBuilder::kNumElementSolutionMatrices;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(num_matrices, std::vector<unsigned int>(2, 0));

    if (nodal_solution)
    {
        // Me
        (*solution_matrix_idx_size_map)[0] = {n_sf_vec, n_sf_vec};
        // AAe
        (*solution_matrix_idx_size_map)[1] = {n_sf_vec, n_sf_vec};
    }
    else
    {
        // JAe
        (*solution_matrix_idx_size_map)[0] = {3 * n_qp, n_sf_vec};
        // JVe
        (*solution_matrix_idx_size_map)[1] = {3 * n_qp, n_sf};
        // BAe
        (*solution_matrix_idx_size_map)[2] = {3 * n_qp, n_sf_vec};
    }
    return solution_matrix_idx_size_map;
}

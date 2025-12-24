#include "element_types.h"
#include "element_builder.h"
#include "variables_map.h"
// Direct constructor
SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(
          keyopts,
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_solid236_fe_structs(fe_structs, keyopts);
          })
{
    set_variables(VariablesMap::get_variables_for_element(*this));
}

// Chain-prolonging constructor for derived types
SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts,
                           std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV3DElementType(
          keyopts,
          [fe_collector, keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_solid236_fe_structs(fe_structs, keyopts);
          })
{
    set_variables(VariablesMap::get_variables_for_element(*this));
}

// Helper to add two FEDataShapeFunctions structs
void SOLID236Type::add_solid236_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                           const std::vector<unsigned int> &keyopts)
{
    fe_structs.push_back({
        FIRST,                                  // order
        NEDELEC_ONE,                            // family
        3,                                      // dimension
        static_cast<Order>(2 * keyopts[2] + 1), // quadrature_order (use default)
        3,                                      // quadrature_dim
        false,                                  // is_infinite
        true                                    // is_vector (scalar FE)
    });
    fe_structs.push_back({
        FIRST,                                  // order
        NEDELEC_ONE,                            // family
        3,                                      // dimension
        static_cast<Order>(2 * keyopts[2] + 1), // quadrature_order (use default)
        3,                                      // quadrature_dim
        false,                                  // is_infinite
        true                                    // is_vector (scalar FE)
    });
    fe_structs.push_back({
        static_cast<Order>(keyopts[2]), // order
        LAGRANGE_VEC,                   // family
        3,                              // dimension
        INVALID_ORDER,                  // quadrature_order (use default)
        3,                              // quadrature_dim
        false,                          // is_infinite
        true                            // is_vector (vector FE)
    });
    // fe_structs.push_back({
    //     SECOND,        // order
    //     LAGRANGE,      // family
    //     3,             // dimension
    //     INVALID_ORDER, // quadrature_order (use default)
    //     3,             // quadrature_dim
    //     false,         // is_infinite
    //     false          // is_vector (vector FE)
    // });
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::get_solution_matrix_idx_size_map(bool nodal_solution)
{
    // TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!     /////
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[3]->n_shape_functions();
    unsigned int n_sf_vec = 3 * n_sf;
    unsigned int n_sf_edge = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();

    unsigned int num_matrices = nodal_solution ? ElementBuilder::NUM_NODAL_SOLUTION_MATRICES
                                               : ElementBuilder::NUM_ELEMENT_SOLUTION_MATRICES;
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(num_matrices, std::vector<unsigned int>(2, 0));

    if (nodal_solution)
    {
        // Me
        (*solution_matrix_idx_size_map)[ElementBuilder::NodalSolutionMatrices::Me] = {n_sf_vec, n_sf_vec};
        // AAe
        (*solution_matrix_idx_size_map)[ElementBuilder::NodalSolutionMatrices::AAe] = {n_sf_vec, n_sf_edge};
    }
    else
    {
        // JAe
        (*solution_matrix_idx_size_map)[ElementBuilder::ElementSolutionMatrices::JAe] = {3 * n_qp, n_sf_edge};
        // JVe
        (*solution_matrix_idx_size_map)[ElementBuilder::ElementSolutionMatrices::JVe] = {3 * n_qp, n_sf};
        // BAe
        (*solution_matrix_idx_size_map)[ElementBuilder::ElementSolutionMatrices::BAe] = {3 * n_qp, n_sf_edge};
    }
    return solution_matrix_idx_size_map;
}

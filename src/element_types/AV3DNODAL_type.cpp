#include "element_types.h"
#include "element_builder.h"
// Direct constructor
AV3DNODALType::AV3DNODALType(const std::vector<unsigned int> &OPTIONs)
    : AV3DElementType(
          OPTIONs,
          [OPTIONs](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_av3dnodal_fe_structs(fe_structs, OPTIONs);
          })
{
}

// Chain-prolonging constructor for derived types
AV3DNODALType::AV3DNODALType(const std::vector<unsigned int> &OPTIONs,
                         std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV3DElementType(
          OPTIONs,
          [fe_collector, OPTIONs](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_av3dnodal_fe_structs(fe_structs, OPTIONs);
          })
{
}

// Helper to add two FEDataShapeFunctions structs
void AV3DNODALType::add_av3dnodal_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                         const std::vector<unsigned int> &OPTIONs)
{
    fe_structs.push_back({
        static_cast<Order>(OPTIONs[2]), // order
        LAGRANGE_VEC,                   // family
        3,                              // dimension
        INVALID_ORDER,                  // quadrature_order (use default)
        3,                              // quadrature_dim
        false,                          // is_infinite
        true                            // is_vector (vector FE)
    });
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DNODALType::get_nodal_solution_matrix_idx_size_map()
{
    unsigned int n_sf = fe_array[1]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[0]->n_shape_functions();

    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(ElementBuilder::kNumNodalSolutionMatrices, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[0] = {n_sf_vec, n_sf_vec};
    // AAe
    (*solution_matrix_idx_size_map)[1] = {n_sf_vec, n_sf_vec};
    return solution_matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV3DNODALType::get_element_solution_matrix_idx_size_map()
{
    unsigned int n_sf = fe_array[1]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[1]->n_quadrature_points();

    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(ElementBuilder::kNumElementSolutionMatrices, std::vector<unsigned int>(2, 0));

    // JAe
    (*solution_matrix_idx_size_map)[0] = {3 * n_qp, n_sf_vec};
    // JVe
    (*solution_matrix_idx_size_map)[1] = {3 * n_qp, n_sf};
    // BAe
    (*solution_matrix_idx_size_map)[2] = {3 * n_qp, n_sf_vec};
    return solution_matrix_idx_size_map;
}

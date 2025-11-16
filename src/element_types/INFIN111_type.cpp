#include "element_types.h"
#include "globals.h"

using namespace Globals;

INFIN111Type::INFIN111Type(const std::vector<unsigned int> &keyopts)
    : SOLID97Type(
          keyopts,
          [](std::vector<Globals::GlobalVariableName> &)
          {
              // No INFIN111-specific variables
          },
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              // Add one FE struct for INFIN111
              fe_structs.push_back({
                  static_cast<Order>(keyopts[2]), // order
                  LAGRANGE,                       // family
                  3,                              // dimension (plane element)
                  INVALID_ORDER,                  // quadrature_order (use default)
                  3,                              // quadrature_dim
                  true,                           // is_infinite
                  false                           // is_vector (scalar FE)
              });
          })
{
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> INFIN111Type::get_solution_matrix_idx_size_map(bool nodal_solution)
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[2]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[2]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(1, std::vector<unsigned int>(2, 1));

    return solution_matrix_idx_size_map;
}
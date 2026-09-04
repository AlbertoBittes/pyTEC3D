/**
 * @file AV2DNODAL_type.cpp
 * @brief AV2DNODALType: FE descriptor set for the AV2D nodal (scalar Az + V) formulation.
 */
#include "apec_m_solver/element_types.h"
#include "apec_m_solver/element_builder.h"

AV2DNODALType::AV2DNODALType(const std::vector<unsigned int> &OPTIONs)
    : AV2DElementType(
          OPTIONs,
          [OPTIONs](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_av2dnodal_fe_structs(fe_structs, OPTIONs);
          })
{
}

AV2DNODALType::AV2DNODALType(const std::vector<unsigned int> &OPTIONs,
                             std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV2DElementType(
          OPTIONs,
          [fe_collector, OPTIONs](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_av2dnodal_fe_structs(fe_structs, OPTIONs);
          })
{
}

void AV2DNODALType::add_av2dnodal_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                             const std::vector<unsigned int> &OPTIONs)
{
    // Example: add a scalar and a vector FE struct
    fe_structs.push_back({
        static_cast<libMesh::Order>(OPTIONs[2]), // order
        libMesh::LAGRANGE,                       // family
        2,                                       // dimension (plane element)
        libMesh::INVALID_ORDER,                  // quadrature_order (use default)
        2,                                       // quadrature_dim
        false,                                   // is_infinite
        false                                    // is_vector (scalar FE)
    });
    fe_structs.push_back({
        static_cast<libMesh::Order>(OPTIONs[2]), // order
        libMesh::LAGRANGE,                       // family
        2,                                       // dimension (plane element)
        libMesh::INVALID_ORDER,                  // quadrature_order (use default)
        2,                                       // quadrature_dim
        false,                                   // is_infinite
        false                                    // is_vector (scalar FE)
    });
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DNODALType::get_nodal_solution_matrix_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(ElementBuilder::kNumNodalSolutionMatrices, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[0] = {2 * n_sf, 2 * n_sf};
    // AAe
    (*solution_matrix_idx_size_map)[1] = {2 * n_sf, n_sf};
    return solution_matrix_idx_size_map;
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> AV2DNODALType::get_element_solution_matrix_idx_size_map()
{
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();

    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(ElementBuilder::kNumElementSolutionMatrices, std::vector<unsigned int>(2, 0));

    // JAe
    (*solution_matrix_idx_size_map)[0] = {n_qp, n_sf};
    // JHe
    (*solution_matrix_idx_size_map)[1] = {1, 1};
    // BAe
    (*solution_matrix_idx_size_map)[2] = {2 * n_qp, n_sf};
    return solution_matrix_idx_size_map;
}

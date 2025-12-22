#include "element_types.h"
#include "element_builder.h"
#include "globals.h"

using namespace Globals;

// Direct constructor
PLANE233Type::PLANE233Type(const std::vector<unsigned int> &keyopts)
    : AV2DElementType(
          keyopts,
          [keyopts](std::vector<Globals::GlobalVariableName> &vars)
          { add_plane233_variables(keyopts, vars); },
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_plane233_fe_structs(fe_structs, keyopts);
          })
{
}

// Chain-prolonging constructor for derived types
PLANE233Type::PLANE233Type(const std::vector<unsigned int> &keyopts,
                           std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                           std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV2DElementType(
          keyopts,
          [variable_collector, keyopts](std::vector<Globals::GlobalVariableName> &vars)
          {
              variable_collector(vars);
              add_plane233_variables(keyopts, vars);
          },
          [fe_collector, keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_plane233_fe_structs(fe_structs, keyopts);
          })
{
}

// Helper to add PLANE233 variables
void PLANE233Type::add_plane233_variables(const std::vector<unsigned int> &keyopts,
                                          std::vector<Globals::GlobalVariableName> &vars)
{
    vars.insert(vars.end(), {GlobalVariableName::A_z,
                             GlobalVariableName::CSG_z});

    PhysicalRegionEnum physical_region = static_cast<PhysicalRegionEnum>(keyopts[1]);
    switch (physical_region)
    {
    case STRANDED_CONDUCTOR:
        vars.insert(vars.end(), {GlobalVariableName::Js_z});
        break;
    case SOLID_CONDUCTOR:
    case STRANDED_COIL:
        vars.insert(vars.end(), {
                                    GlobalVariableName::VOLT,
                                    GlobalVariableName::AMPS,
                                });
        if (physical_region == STRANDED_COIL)
            vars.insert(vars.end(), {GlobalVariableName::EMF,
                                     GlobalVariableName::CURT});
        break;
    case STRANDED_COIL_CIRCUIT_COUPLED:
    case SOLID_CONDUCTOR_CIRCUIT_COUPLED:
        vars.insert(vars.end(), {GlobalVariableName::V,
                                 GlobalVariableName::AMPS});
        if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
            vars.insert(vars.end(), {GlobalVariableName::EMF,
                                     GlobalVariableName::CURT});
        break;
    default:
        break;
    }
    if (keyopts[3] == 1) // PSI option
    {
        vars.insert(vars.end(), {GlobalVariableName::PSI});
        vars.insert(vars.end(), {GlobalVariableName::PSIe});
    }
}

// Helper to add two FEDataShapeFunctions structs
void PLANE233Type::add_plane233_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                           const std::vector<unsigned int> &keyopts)
{
    // Example: add a scalar and a vector FE struct
    fe_structs.push_back({
        static_cast<Order>(keyopts[2]), // order
        LAGRANGE,                       // family
        2,                              // dimension (plane element)
        INVALID_ORDER,                  // quadrature_order (use default)
        2,                              // quadrature_dim
        false,                          // is_infinite
        false                           // is_vector (scalar FE)
    });
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> PLANE233Type::get_solution_matrix_idx_size_map(bool nodal_solution)
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();

    unsigned int num_matrices = nodal_solution ? ElementBuilder::NUM_NODAL_SOLUTION_MATRICES
                                               : ElementBuilder::NUM_ELEMENT_SOLUTION_MATRICES;

    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(num_matrices, std::vector<unsigned int>(2, 0));

    if (nodal_solution)
    {
        // Me
        (*solution_matrix_idx_size_map)[ElementBuilder::NodalSolutionMatrices::Me] = {2 * n_sf, 2 * n_sf};
        // AAe
        (*solution_matrix_idx_size_map)[ElementBuilder::NodalSolutionMatrices::AAe] = {2 * n_sf, n_sf};
    }
    else
    {
        // JAe
        (*solution_matrix_idx_size_map)[ElementBuilder::ElementSolutionMatrices::JAe] = {n_qp, n_sf};
        // JHe
        (*solution_matrix_idx_size_map)[ElementBuilder::ElementSolutionMatrices::JVe] = {1, 1};
        // BAe
        (*solution_matrix_idx_size_map)[ElementBuilder::ElementSolutionMatrices::BAe] = {2 * n_qp, n_sf};
    }
    return solution_matrix_idx_size_map;
}
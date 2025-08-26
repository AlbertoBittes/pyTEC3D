#include "element_types.h"
#include "element_builder.h"
#include "globals.h"

using namespace Globals;

// Direct constructor
SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(
          keyopts,
          [keyopts](std::vector<Globals::GlobalVariableName> &vars)
          { add_solid236_variables(static_cast<PhysicalRegionEnum>(keyopts[1]), vars); },
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_solid236_fe_structs(fe_structs, keyopts);
          })
{
}

// Chain-prolonging constructor for derived types
SOLID236Type::SOLID236Type(const std::vector<unsigned int> &keyopts,
                           std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                           std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV3DElementType(
          keyopts,
          [variable_collector, keyopts](std::vector<Globals::GlobalVariableName> &vars)
          {
              variable_collector(vars);
              add_solid236_variables(static_cast<PhysicalRegionEnum>(keyopts[1]), vars);
          },
          [fe_collector, keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_solid236_fe_structs(fe_structs, keyopts);
          })
{
}

// Helper to add SOLID236 variables
void SOLID236Type::add_solid236_variables(const PhysicalRegionEnum physical_region,
                                          std::vector<Globals::GlobalVariableName> &vars)
{
    vars.insert(vars.end(), {GlobalVariableName::Az});

    switch (physical_region)
    {
    case STRANDED_CONDUCTOR:
        vars.insert(vars.end(), {GlobalVariableName::Js_x,
                                 GlobalVariableName::Js_y,
                                 GlobalVariableName::Js_z});
        break;
    case SOLID_CONDUCTOR:
    case STRANDED_COIL:
        vars.insert(vars.end(), {GlobalVariableName::VOLT,
                                 GlobalVariableName::AMPS});
        if (physical_region == STRANDED_COIL)
        {
            vars.insert(vars.end(), {GlobalVariableName::EMF,
                                     GlobalVariableName::CURT});
        }
        break;
    default:
        break;
    }
    // Add any further SOLID236-specific logic here if needed
}

// Helper to add two FEDataShapeFunctions structs
void SOLID236Type::add_solid236_fe_structs(std::vector<FEDataShapeFunctions> &fe_structs,
                                           const std::vector<unsigned int> &keyopts)
{
    fe_structs.push_back({
        FIRST,       // order
        NEDELEC_ONE, // family
        3,           // dimension
        THIRD,       // quadrature_order (use default)
        3,           // quadrature_dim
        false,       // is_infinite
        true         // is_vector (scalar FE)
    });
    fe_structs.push_back({
        FIRST,       // order
        NEDELEC_ONE, // family
        3,           // dimension
        THIRD,       // quadrature_order (use default)
        3,           // quadrature_dim
        false,       // is_infinite
        true         // is_vector (scalar FE)
    });
    fe_structs.push_back({
        SECOND,        // order
        LAGRANGE_VEC,  // family
        3,             // dimension
        INVALID_ORDER, // quadrature_order (use default)
        3,             // quadrature_dim
        false,         // is_infinite
        true           // is_vector (vector FE)
    });
}

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID236Type::get_solution_matrix_idx_size_map()
{
    // TEMPORARY !!!!!!!!!!!!!!!!!!!!!!!!     /////
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[0]->n_shape_functions();
    unsigned int n_sf_vec = 3 * n_sf;
    unsigned int n_sf_edge = fe_array[1]->n_shape_functions();
    unsigned int n_qp = fe_array[0]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(5, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::Me_SOLID97] = {n_sf_vec, n_sf_vec};
    // AAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::AAe_SOLID97] = {n_sf_vec, n_sf_edge};
    // JAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JAe_SOLID97] = {3 * n_qp, n_sf_edge};
    // JVe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JVe_SOLID97] = {3 * n_qp, n_sf};
    // BAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::BAe_SOLID97] = {3 * n_qp, n_sf_edge};

    return solution_matrix_idx_size_map;
}
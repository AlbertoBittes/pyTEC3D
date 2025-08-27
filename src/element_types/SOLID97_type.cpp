#include "element_types.h"
#include "element_builder.h"
#include "globals.h"

using namespace Globals;
// Direct constructor
SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts)
    : AV3DElementType(
          keyopts,
          [keyopts](std::vector<Globals::GlobalVariableName> &vars)
          { add_solid97_variables(static_cast<PhysicalRegionEnum>(keyopts[1]), vars); },
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              add_solid97_fe_structs(fe_structs, keyopts);
          })
{
}

// Chain-prolonging constructor for derived types
SOLID97Type::SOLID97Type(const std::vector<unsigned int> &keyopts,
                         std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                         std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AV3DElementType(
          keyopts,
          [variable_collector, keyopts](std::vector<Globals::GlobalVariableName> &vars)
          {
              variable_collector(vars);
              add_solid97_variables(static_cast<PhysicalRegionEnum>(keyopts[1]), vars);
          },
          [fe_collector, keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              add_solid97_fe_structs(fe_structs, keyopts);
          })
{
}

// Helper to add SOLID97 variables
void SOLID97Type::add_solid97_variables(const PhysicalRegionEnum physical_region,
                                        std::vector<GlobalVariableName> &vars)
{
    vars.insert(vars.end(), {GlobalVariableName::A_x,
                             GlobalVariableName::A_y,
                             GlobalVariableName::A_z,
                             GlobalVariableName::CSG_x,
                             GlobalVariableName::CSG_y,
                             GlobalVariableName::CSG_z});

    switch (physical_region)
    {
    case SOLID_CONDUCTOR:
    case SOLID_CONDUCTOR_CIRCUIT_COUPLED:
        vars.insert(vars.end(), {GlobalVariableName::VOLT,
                                 GlobalVariableName::AMPS});
        if (physical_region == SOLID_CONDUCTOR_CIRCUIT_COUPLED)
        {
            vars.insert(vars.end(), {GlobalVariableName::EMF,
                                     GlobalVariableName::CURT});
        }
        break;
    case STRANDED_COIL_CIRCUIT_COUPLED:
    case STRANDED_COIL_VOLTAGE_FED:
        vars.insert(vars.end(), {GlobalVariableName::CURR,
                                 GlobalVariableName::VLTG});
        if (physical_region == STRANDED_COIL_CIRCUIT_COUPLED)
        {
            vars.insert(vars.end(), {GlobalVariableName::EMF,
                                     GlobalVariableName::CURT});
        }
        break;
    case STRANDED_CONDUCTOR:
        vars.insert(vars.end(), {GlobalVariableName::Js_x,
                                 GlobalVariableName::Js_y,
                                 GlobalVariableName::Js_z});
        break;
    default:
        break;
    }
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

std::unique_ptr<std::vector<std::vector<unsigned int>>> SOLID97Type::get_solution_matrix_idx_size_map()
{
    // Check this when you have periodic BCs
    unsigned int n_sf = fe_array[2]->n_shape_functions();
    unsigned int n_sf_vec = fe_array[0]->n_shape_functions();
    unsigned int n_qp = fe_array[2]->n_quadrature_points();
    std::unique_ptr<std::vector<std::vector<unsigned int>>> solution_matrix_idx_size_map =
        std::make_unique<std::vector<std::vector<unsigned int>>>(5, std::vector<unsigned int>(2, 0));

    // Me
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::Me_SOLID97] = {n_sf_vec, n_sf_vec};
    // AAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::AAe_SOLID97] = {n_sf_vec, n_sf_vec};
    // JAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JAe_SOLID97] = {3 * n_qp, n_sf_vec};
    // JVe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::JVe_SOLID97] = {3 * n_qp, n_sf};
    // BAe
    (*solution_matrix_idx_size_map)[ElementBuilder::SOLID97SolutionMatrices::BAe_SOLID97] = {3 * n_qp, n_sf_vec};

    return solution_matrix_idx_size_map;
}
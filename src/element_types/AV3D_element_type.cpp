#include "element_types.h"
#include "globals.h"

using namespace Globals;

AV3DElementType::AV3DElementType(const std::vector<unsigned int> &keyopts,
                                 std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AVElementType(
          keyopts,
          [variable_collector](std::vector<Globals::GlobalVariableName> &vars)
          {
              variable_collector(vars);
              vars.insert(vars.end(), {GlobalVariableName::B_z,
                                       GlobalVariableName::H_z,
                                       GlobalVariableName::Jt_x,
                                       GlobalVariableName::Jt_y,
                                       GlobalVariableName::Be_z,
                                       GlobalVariableName::He_z});
          },
          [fe_collector, keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              fe_structs.push_back({
                  static_cast<Order>(keyopts[2]), // order
                  LAGRANGE,                       // family
                  3,                              // dimension
                  INVALID_ORDER,                  // quadrature_order
                  3,                              // quadrature_dim
                  false,                          // is_infinite
                  false                           // is_vector
              });
          })
{
    // Optionally call set_fe_array() or other initialization here
}
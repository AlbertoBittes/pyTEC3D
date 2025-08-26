#include "element_types.h"
#include "globals.h"

using namespace Globals;

INFIN110Type::INFIN110Type(const std::vector<unsigned int> &keyopts)
    : PLANE233Type(
          keyopts,
          [](std::vector<Globals::GlobalVariableName> &)
          {
              // No INFIN110-specific variables
          },
          [keyopts](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              // Add one FE struct for INFIN110
              fe_structs.push_back({
                  static_cast<Order>(keyopts[2]), // order
                  LAGRANGE,                       // family
                  2,                              // dimension (plane element)
                  INVALID_ORDER,                  // quadrature_order (use default)
                  2,                              // quadrature_dim
                  true,                           // is_infinite
                  false                           // is_vector (scalar FE)
              });
          })
{
}
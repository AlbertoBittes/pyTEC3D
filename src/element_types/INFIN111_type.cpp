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
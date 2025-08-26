#include "element_types.h"
#include "globals.h"

using namespace Globals;

INFIN112Type::INFIN112Type(const std::vector<unsigned int> &keyopts)
    : SOLID236Type(
          keyopts,
          [](std::vector<Globals::GlobalVariableName> &)
          {
              // No INFIN112-specific variables
          },
          [](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              // Add one FE struct for INFIN112
              fe_structs.push_back({
                  FIRST,         // order
                  NEDELEC_ONE,   // family
                  3,             // dimension (plane element)
                  INVALID_ORDER, // quadrature_order (use default)
                  3,             // quadrature_dim
                  true,          // is_infinite
                  false          // is_vector (scalar FE)
              });
          })
{
}
#include "element_types.h"
INFIN112Type::INFIN112Type(const std::vector<unsigned int> &keyopts)
    : SOLID236Type(
          keyopts,
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

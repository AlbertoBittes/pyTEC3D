#include "element_types.h"

AV3DElementType::AV3DElementType(const std::vector<unsigned int> &OPTIONs,
                                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AVElementType(
          OPTIONs,
          [fe_collector, OPTIONs](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              fe_structs.push_back({
                  static_cast<Order>(OPTIONs[2]), // order
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

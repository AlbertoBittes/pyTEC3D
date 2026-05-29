/**
 * @file AV3D_element_type.cpp
 * @brief AV3DElementType: base class for all 3D AV formulation element types.
 */
#include "element_types.h"

AV3DElementType::AV3DElementType(const std::vector<unsigned int> &OPTIONs,
                                 std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : AVElementType(
          OPTIONs,
          [fe_collector, OPTIONs](std::vector<FEDataShapeFunctions> &fe_structs)
          {
              fe_collector(fe_structs);
              fe_structs.push_back({
                  static_cast<libMesh::Order>(OPTIONs[2]), // order
                  libMesh::LAGRANGE,                       // family
                  3,                                       // dimension
                  libMesh::INVALID_ORDER,                  // quadrature_order
                  3,                                       // quadrature_dim
                  false,                                   // is_infinite
                  false                                    // is_vector
              });
          })
{
    // Optionally call set_fe_array() or other initialization here
}

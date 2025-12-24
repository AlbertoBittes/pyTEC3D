#include "element_types.h"
// Constructor
AVElementType::AVElementType(const std::vector<unsigned int> &keyopts,
                             std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : ElementType(keyopts, fe_collector)
{
    // Optionally call set_fe_array() or other initialization here
}

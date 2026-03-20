#include "element_types.h"
// Constructor
AVElementType::AVElementType(const std::vector<unsigned int> &OPTIONs,
                             std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : ElementType(OPTIONs, fe_collector)
{
    // Optionally call set_fe_array() or other initialization here
}

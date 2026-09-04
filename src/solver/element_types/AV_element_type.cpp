/**
 * @file AV_element_type.cpp
 * @brief AVElementType: intermediate base for AV formulation element types.
 */
#include "apec_m_solver/element_types.h"

AVElementType::AVElementType(const std::vector<unsigned int> &OPTIONs,
                             std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : ElementType(OPTIONs, fe_collector)
{
    // Optionally call set_fe_array() or other initialization here
}

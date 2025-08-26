#include "element_types.h"
#include "globals.h"

using namespace Globals;
// Constructor
AVElementType::AVElementType(const std::vector<unsigned int> &keyopts,
                             std::function<void(std::vector<Globals::GlobalVariableName> &)> variable_collector,
                             std::function<void(std::vector<FEDataShapeFunctions> &)> fe_collector)
    : ElementType(
          keyopts,
          [variable_collector](std::vector<Globals::GlobalVariableName> &vars)
          {
              variable_collector(vars);
              vars.insert(vars.end(), {GlobalVariableName::B_x,
                                       GlobalVariableName::B_y,
                                       GlobalVariableName::H_x,
                                       GlobalVariableName::H_y,
                                       GlobalVariableName::Jt_z,
                                       GlobalVariableName::Be_x,
                                       GlobalVariableName::Be_y,
                                       GlobalVariableName::He_x,
                                       GlobalVariableName::He_y});
          },
          fe_collector // Just forward, do not add anything here
      )
{
    // Optionally call set_fe_array() or other initialization here
}
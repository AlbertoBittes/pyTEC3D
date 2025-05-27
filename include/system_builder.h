#ifndef SYSTEM_BUILDER_H
#define SYSTEM_BUILDER_H

#include <string>
#include "libmesh/equation_systems.h"
#include "fe_data_object.h"
#include "fem_service_abstract.h"
#include "globals.h"

using namespace libMesh;

class SystemBuilder : public System::Initialization, public FEDataObjectBase<System>, public FEMServiceAbstract
{
public:
    // Default constructor
    using FEDataObjectBase<System>::FEDataObjectBase;

    // Override the initialize method
    void initialize() override;
};

#endif // SYSTEM_BUILDER_H
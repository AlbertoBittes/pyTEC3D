#ifndef APEC_M_LINEAR_SYSTEM_H
#define APEC_M_LINEAR_SYSTEM_H

#include <memory>
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "fem_service_abstract.h"

// Add other necessary includes here

class ApecMLinearSystem : public FEMServiceAbstract, public TransientSystem<LinearImplicitSystem>
{
public:
    using TransientSystem<LinearImplicitSystem>::TransientSystem;

private:
    // Add private members here
};

#endif // APEC_M_LINEAR_SYSTEM_H
#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/system.h>
#include <libmesh/equation_systems.h>
#include "fem_service_abstract.h"

using namespace libMesh;

class BoundaryConditionsObject : public System::Constraint, public FEMServiceAbstract
{
public:
    // Default constructor
    BoundaryConditionsObject() = default;

    // Destructor
    virtual ~BoundaryConditionsObject() = default;

    // Other methods related to boundary conditions can be added here
    void constrain() override;

    bool get_update_boundary_conditions() const
    {
        // This method should return whether the boundary conditions need to be updated
        return update_boundary_conditions;
    }

private:
    bool update_boundary_conditions = false; // Flag to indicate if boundary conditions need updating
    // Other private members related to boundary conditions can be added here
};
#endif // BOUNDARY_CONDITIONS_OBJECT_H
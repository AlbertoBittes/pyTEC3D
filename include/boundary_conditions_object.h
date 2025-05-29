#ifndef BOUNDARY_CONDITIONS_OBJECT_H
#define BOUNDARY_CONDITIONS_OBJECT_H

#include <libmesh/system.h>

class BoundaryConditionsObject : public System::Constraint
{
public:
    // Default constructor
    BoundaryConditionsObject() = default;

    // Destructor
    virtual ~BoundaryConditionsObject() = default;

    // Other methods related to boundary conditions can be added here
    void constrain() override;
};
#endif // BOUNDARY_CONDITIONS_OBJECT_H
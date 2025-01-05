#ifndef NODE_DOF_CONSTRAINTS_H
#define NODE_DOF_CONSTRAINTS_H

#include "vtkPointData.h"
#include <libmesh/equation_systems.h>
#include <string>
#include "../include/time_stepper.h"

// Forward declaration from libMesh
using namespace libMesh;

class NodalBCConstraint : public System::Constraint
{
private:
    System &_sys;
    vtkPointData *_pointData; // New attribute to store the vtkPointData pointer
    unsigned int _dim;
    std::shared_ptr<TimeStepper> timeStepper;

public:
    // Modified constructor to take a vtkPointData pointer
    NodalBCConstraint(TransientLinearImplicitSystem &system, vtkPointData *pointData, std::shared_ptr<TimeStepper> timeStepper)
        : _sys(system), _pointData(pointData), timeStepper(timeStepper)
    {
        _dim = _sys.get_mesh().mesh_dimension();
    }

    virtual ~NodalBCConstraint() {}

    void constrain();
};

#endif // NODE_DOF_CONSTRAINTS_H
#ifndef NODE_DOF_CONSTRAINTS_H
#define NODE_DOF_CONSTRAINTS_H

#include "vtkPointData.h"
#include <libmesh/equation_systems.h>
#include <string>

// Forward declaration from libMesh
using namespace libMesh;

class NodalBCConstraint : public System::Constraint
{
private:
    System &_sys;
    vtkPointData *_pointData; // New attribute to store the vtkPointData pointer
    unsigned int _dim;

public:
    // Modified constructor to take a vtkPointData pointer
    NodalBCConstraint(System &sys, vtkPointData *pointData)
        : Constraint(), _sys(sys), _pointData(pointData)
    {
        _dim = _sys.get_mesh().mesh_dimension();
    }

    virtual ~NodalBCConstraint() {}

    void constrain();
};

#endif // NODE_DOF_CONSTRAINTS_H
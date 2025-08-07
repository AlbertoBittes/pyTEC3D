#include <libmesh/libmesh.h>

#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/tensor_value.h"

#ifndef AZIMUTHAL_PERIODIC_BOUNDARY_H
#define AZIMUTHAL_PERIODIC_BOUNDARY_H

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifdef LIBMESH_ENABLE_PERIODIC
// Here we define the azimuthal periodic boundary condition class.
// This class assumes that (u,v,w) are variables 0, 1, and 2 in
// the System, but this could be made more general if needed.
class AzimuthalPeriodicBoundary : public PeriodicBoundaryBase
{
public:
    AzimuthalPeriodicBoundary(
        Point center,
        Point axis,
        Real angle)
        : PeriodicBoundaryBase(),
          _center(center),
          _axis(axis),
          _theta(angle)
    {
        set_variable(0);
        // set_variable(1);
        // set_variable(2);
        set_up_rotation_matrix();
    }

    AzimuthalPeriodicBoundary(const AzimuthalPeriodicBoundary &o, TransformationType t = FORWARD)
        : PeriodicBoundaryBase(o),
          _center(o._center),
          _axis(o._axis),
          _theta(o._theta)
    {
        if (t == INVERSE)
        {
            std::swap(myboundary, pairedboundary);
            _theta *= -1.0;
        }

        set_variable(0);
        // set_variable(1);
        // set_variable(2);
        set_up_rotation_matrix();
    }

    virtual ~AzimuthalPeriodicBoundary() = default;

    void set_up_rotation_matrix();

    virtual Point get_corresponding_pos(const Point &pt) const override
    {
        // Note that since _theta defines the angle from "paired boundary" to
        // "my boundary", and we want the inverse of that here, we multiply
        // by the transpose of the transformation matrix below.
        return _R.left_multiply(pt - _center) + _center;
    }

    virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const override
    {
        return std::make_unique<AzimuthalPeriodicBoundary>(*this, t);
    }

private:
    // Define the properties needed for the azimuthal periodic boundary.
    // Note that _theta specifies the angle from "paired boundary" to
    // "my boundary".
    Point _center;
    Point _axis;
    Real _theta;
    RealTensorValue _R;
};

#endif // LIBMESH_ENABLE_PERIODIC

#endif // AZIMUTHAL_PERIODIC_BOUNDARY_H
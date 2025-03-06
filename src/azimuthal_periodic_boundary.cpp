#include "../include/azimuthal_periodic_boundary.h"

void AzimuthalPeriodicBoundary::set_up_rotation_matrix()
{
    // Formula for rotation matrix about an axis is given on wikipedia:
    // en.wikipedia.org/wiki/Rotation_matrix
    // We rotate by angle theta about the axis defined by u, which is a
    // unit vector in the direction of _axis.
    // Note: this example already requires LIBMESH_DIM==3 so we don't
    // explicitly check for that here.
    Point u = _axis.unit();
    Real u_x = u(0);
    Real u_y = u(1);
    Real u_z = u(2);
    Real cost = std::cos(_theta);
    Real sint = std::sin(_theta);
    _R = RealTensorValue(cost + u_x * u_x * (1.0 - cost), u_x * u_y * (1.0 - cost) - u_z * sint, u_x * u_z * (1.0 - cost) + u_y * sint,
                         u_y * u_x * (1.0 - cost) + u_z * sint, cost + u_y * u_y * (1.0 - cost), u_y * u_z * (1.0 - cost) - u_x * sint,
                         u_z * u_x * (1.0 - cost) - u_y * sint, u_z * u_y * (1.0 - cost) + u_x * sint, cost + u_z * u_z * (1.0 - cost));

    // For generality, PeriodicBoundaryBase::set_transformation_matrix() takes an
    // (n_vars * n_vars) DenseMatrix which we construct on the fly.
    // this->set_transformation_matrix(DenseMatrix<Real>(3, 3, {_R(0, 0), _R(0, 1), _R(0, 2), _R(1, 0), _R(1, 1), _R(1, 2), _R(2, 0), _R(2, 1), _R(2, 2)}));
}

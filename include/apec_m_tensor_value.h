#ifndef APEC_M_TENSOR_VALUE_H
#define APEC_M_TENSOR_VALUE_H

#include <array>
#include <cmath>

#include <libmesh/libmesh.h>
#include <libmesh/tensor_value.h>

using namespace libMesh;

class ApecMTensorValue : public TensorValue<Real>
{
public:
    using TensorValue<Real>::TensorValue;

    static Real snap_trig_value(Real value)
    {
        const Real eps = 1e-12;
        if (std::abs(value) < eps)
            return 0.;
        if (std::abs(value - 1.) < eps)
            return 1.;
        if (std::abs(value + 1.) < eps)
            return -1.;
        return value;
    }

    static ApecMTensorValue build_ansys_rotation_matrix(Real angle_xy_deg,
                                                        Real angle_yz_deg,
                                                        Real angle_xz_deg)
    {
        const Real to_rad = pi / 180.0;
        const Real angle_xy = angle_xy_deg * to_rad;
        const Real angle_yz = angle_yz_deg * to_rad;
        const Real angle_xz = angle_xz_deg * to_rad;

        const Real c1 = snap_trig_value(std::cos(angle_xy));
        const Real s1 = snap_trig_value(std::sin(angle_xy));
        const Real c2 = snap_trig_value(std::cos(angle_yz));
        const Real s2 = snap_trig_value(std::sin(angle_yz));
        const Real c3 = snap_trig_value(std::cos(angle_xz));
        const Real s3 = snap_trig_value(std::sin(angle_xz));

        ApecMTensorValue rotation;
        // ANSYS convention: intrinsic rotations about Z, then X', then Y''
        rotation(0, 0) = c1 * c3 - s1 * s2 * s3;
        rotation(0, 1) = -s1 * c2;
        rotation(0, 2) = c1 * s3 + s1 * s2 * c3;
        rotation(1, 0) = s1 * c3 + c1 * s2 * s3;
        rotation(1, 1) = c1 * c2;
        rotation(1, 2) = s1 * s3 - c1 * s2 * c3;
        rotation(2, 0) = -c2 * s3;
        rotation(2, 1) = s2;
        rotation(2, 2) = c2 * c3;
        return rotation;
    }
};

#endif // APEC_M_TENSOR_VALUE_H

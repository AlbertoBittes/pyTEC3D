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

    static Real det2x2(Real a00, Real a01, Real a10, Real a11)
    {
        return a00 * a11 - a01 * a10;
    }

    static bool invert2x2(Real a00,
                          Real a01,
                          Real a10,
                          Real a11,
                          std::array<Real, 4> &inv,
                          Real tol = 1e-12)
    {
        const Real det = det2x2(a00, a01, a10, a11);
        if (std::abs(det) < tol)
            return false;
        inv[0] = a11 / det;
        inv[1] = -a01 / det;
        inv[2] = -a10 / det;
        inv[3] = a00 / det;
        return true;
    }

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

class ApecMTensor2Value
{
public:
    ApecMTensor2Value()
        : data{{0., 0., 0., 0.}}
    {
    }

    ApecMTensor2Value(Real a00, Real a01, Real a10, Real a11)
        : data{{a00, a01, a10, a11}}
    {
    }

    ApecMTensor2Value(const ApecMTensorValue &tensor,
                      unsigned int row0,
                      unsigned int row1,
                      unsigned int col0,
                      unsigned int col1)
        : data{{tensor(row0, col0), tensor(row0, col1), tensor(row1, col0), tensor(row1, col1)}}
    {
    }

    Real &operator()(unsigned int row, unsigned int col)
    {
        return data[row * 2 + col];
    }

    Real operator()(unsigned int row, unsigned int col) const
    {
        return data[row * 2 + col];
    }

    Real det() const
    {
        return data[0] * data[3] - data[1] * data[2];
    }

    bool inverse(ApecMTensor2Value &inv, Real tol = 1e-12) const
    {
        const Real determinant = det();
        if (std::abs(determinant) < tol)
            return false;

        inv(0, 0) = data[3] / determinant;
        inv(0, 1) = -data[1] / determinant;
        inv(1, 0) = -data[2] / determinant;
        inv(1, 1) = data[0] / determinant;
        return true;
    }

private:
    std::array<Real, 4> data;
};

class ApecMVector2Value
{
public:
    ApecMVector2Value()
        : data{{0., 0.}}
    {
    }

    ApecMVector2Value(Real x, Real y)
        : data{{x, y}}
    {
    }

    Real &operator()(unsigned int i)
    {
        return data[i];
    }

    Real operator()(unsigned int i) const
    {
        return data[i];
    }

private:
    std::array<Real, 2> data;
};

#endif // APEC_M_TENSOR_VALUE_H

#ifndef COORDINATE_SYSTEM_H
#define COORDINATE_SYSTEM_H

#include <memory>

#include <libmesh/point.h>
#include <libmesh/vector_value.h>

#include "apec_m_tensor_value.h"

class CoordinateSystem
{
public:
    enum class Type
    {
        Cartesian = 0,
        Cylindrical = 1,
        Spherical = 2
    };

    virtual ~CoordinateSystem() = default;

    Type type() const { return type_; }
    const Point &origin() const { return origin_; }

    Point axis_direction() const;
    virtual Point to_local_coordinates(const Point &global_point) const = 0;

    virtual ApecMTensorValue global_to_local_rotation(const Point &global_point) const = 0;

    ApecMTensorValue local_to_global_rotation(const Point &global_point) const
    {
        return global_to_local_rotation(global_point).transpose();
    }

protected:
    CoordinateSystem(Type type,
                     const Point &origin,
                     const Point &nxax,
                     const Point &nxypl);

    CoordinateSystem(Type type,
                     const Point &origin,
                     Real thxy_deg,
                     Real thyz_deg,
                     Real thzx_deg);

    const ApecMTensorValue &local_to_global_cartesian() const { return local_to_global_cart_; }
    VectorValue<Real> to_local_vector(const Point &global_point) const;

private:
    Type type_;
    Point origin_;
    ApecMTensorValue local_to_global_cart_;
};

class CartesianCoordinateSystem : public CoordinateSystem
{
public:
    CartesianCoordinateSystem(const Point &origin = Point(0.0, 0.0, 0.0),
                              const Point &nxax = Point(1.0, 0.0, 0.0),
                              const Point &nxypl = Point(0.0, 1.0, 0.0));

    CartesianCoordinateSystem(const Point &origin,
                              Real thxy_deg,
                              Real thyz_deg,
                              Real thzx_deg);

    ApecMTensorValue global_to_local_rotation(const Point &global_point) const override;
    Point to_local_coordinates(const Point &global_point) const override;
};

class CylindricalCoordinateSystem : public CoordinateSystem
{
public:
    CylindricalCoordinateSystem(const Point &origin = Point(0.0, 0.0, 0.0),
                                const Point &nxax = Point(1.0, 0.0, 0.0),
                                const Point &nxypl = Point(0.0, 1.0, 0.0));

    CylindricalCoordinateSystem(const Point &origin,
                                Real thxy_deg,
                                Real thyz_deg,
                                Real thzx_deg);

    ApecMTensorValue global_to_local_rotation(const Point &global_point) const override;
    Point to_local_coordinates(const Point &global_point) const override;
};

class SphericalCoordinateSystem : public CoordinateSystem
{
public:
    SphericalCoordinateSystem(const Point &origin = Point(0.0, 0.0, 0.0),
                              const Point &nxax = Point(1.0, 0.0, 0.0),
                              const Point &nxypl = Point(0.0, 1.0, 0.0));

    SphericalCoordinateSystem(const Point &origin,
                              Real thxy_deg,
                              Real thyz_deg,
                              Real thzx_deg);

    ApecMTensorValue global_to_local_rotation(const Point &global_point) const override;
    Point to_local_coordinates(const Point &global_point) const override;
};

std::unique_ptr<CoordinateSystem> make_coordinate_system_local(int kcs,
                                                               const Point &origin,
                                                               Real thxy_deg,
                                                               Real thyz_deg,
                                                               Real thzx_deg);

std::unique_ptr<CoordinateSystem> make_coordinate_system_cs(int kcs,
                                                            const Point &norig,
                                                            const Point &nxax,
                                                            const Point &nxypl);

#endif // COORDINATE_SYSTEM_H

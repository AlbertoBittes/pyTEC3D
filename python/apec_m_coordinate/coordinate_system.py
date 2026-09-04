"""Coordinate-system API backed by the standalone C++ coordinate extension."""

from .apec_m_coordinate_system import (
    CoordinateSystem,
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    XAxisCylindricalCoordinateSystem,
    YAxisCylindricalCoordinateSystem,
    ZAxisCylindricalCoordinateSystem,
    SphericalCoordinateSystem,
)

__all__ = [
    "CoordinateSystem",
    "CartesianCoordinateSystem",
    "CylindricalCoordinateSystem",
    "XAxisCylindricalCoordinateSystem",
    "YAxisCylindricalCoordinateSystem",
    "ZAxisCylindricalCoordinateSystem",
    "SphericalCoordinateSystem",
    "COORDINATE_CLASS_TO_ID",
    "COORDINATE_ID_TO_CLASS",
    "coordinate_class_from_id",
]


COORDINATE_CLASS_TO_ID = {
    "CartesianCoordinateSystem": 1,
    "CylindricalCoordinateSystem": 2,
    "XAxisCylindricalCoordinateSystem": 3,
    "YAxisCylindricalCoordinateSystem": 4,
    "ZAxisCylindricalCoordinateSystem": 5,
    "SphericalCoordinateSystem": 6,
}

COORDINATE_ID_TO_CLASS = {
    1: CartesianCoordinateSystem,
    2: CylindricalCoordinateSystem,
    3: XAxisCylindricalCoordinateSystem,
    4: YAxisCylindricalCoordinateSystem,
    5: ZAxisCylindricalCoordinateSystem,
    6: SphericalCoordinateSystem,
}


def coordinate_class_from_id(class_id: int):
    return COORDINATE_ID_TO_CLASS.get(int(class_id))

"""Coordinate systems: Cartesian, cylindrical and spherical, placed anywhere in space.

A coordinate system is an origin and a rotation: the 3x3 matrix ``R`` whose columns are its local x, y, z axes in global
coordinates (``local_to_global_cart``). On that local Cartesian frame the cylindrical system uses (r, theta, z) and
the spherical one (r, theta, phi); angles are in degrees.

The orientation is given by one of:

* three angles in degrees ``thxy_deg, thyz_deg, thzx_deg`` (the APDL ``LOCAL`` convention);
* two points: ``nxax`` on the local x axis, ``nxypl`` in the local xy plane;
* the rotation matrix ``local_to_global_cart``;
* nothing: the global axes, rotated as the type does by default.

The model data stores a coordinate system as ``{"type", "origin", "matrix"}`` (:func:`coordinate_system_to_dict`);
:func:`coordinate_system_from_dict` reads it back, and also accepts ``angles`` or ``points``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar, Self

import numpy as np

#: Orientation angles (thxy, thyz, thzx) of the cylindrical systems whose axis is the global x, y or z axis.
X_AXIS_CYLINDRICAL_ANGLES = (0.0, 0.0, 90.0)
Y_AXIS_CYLINDRICAL_ANGLES = (0.0, -90.0, 0.0)
Z_AXIS_CYLINDRICAL_ANGLES = (0.0, 0.0, 0.0)

#: Trigonometric values closer than this to 0 or +-1 are snapped, so that right angles give exact matrices.
_SNAP_TOLERANCE = 1e-12


def _vector3(values: Any, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (3,):
        raise ValueError(f"{name} must have 3 components, got shape {array.shape}.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite.")
    return array


def _points3(values: Any, name: str) -> np.ndarray:
    """Points or vectors as a finite (N, 3) array; a single (3,) input becomes (1, 3)."""
    array = np.asarray(values, dtype=float)
    if array.ndim == 1 and array.shape[0] == 3:
        array = array.reshape(1, 3)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"{name} must have shape (N, 3) or (3,), got {array.shape}.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite.")
    return array


def _rotation_matrix(values: Any, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3), got {array.shape}.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite.")
    if not np.allclose(array.T @ array, np.eye(3), atol=1e-8):
        raise ValueError(f"{name} must be orthonormal (R^T R = I).")
    return array


def _unit(vector: np.ndarray, message: str) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if norm <= 0.0:
        raise ValueError(message)
    return vector / norm


def rotation_coefficients_from_angles(angles_deg: Any) -> np.ndarray:
    """Local-to-global rotation matrices of angle triplets, flattened row by row.

    Args:
        angles_deg: Shape (N, 3): ``thxy, thyz, thzx`` in degrees (APDL ``LOCAL`` convention).

    Returns:
        Shape (N, 9): ``R[0, 0], R[0, 1], ..., R[2, 2]`` for each triplet.

    Raises:
        ValueError: If the shape is not (N, 3) or a value is not finite.
    """
    angles = np.asarray(angles_deg, dtype=float)
    if angles.ndim != 2 or angles.shape[1] != 3:
        raise ValueError(f"angles_deg must have shape (N, 3), got {angles.shape}.")
    if not np.all(np.isfinite(angles)):
        raise ValueError("angles_deg must be finite.")

    def snapped(values: np.ndarray) -> np.ndarray:
        out = values.copy()
        out[np.abs(out) < _SNAP_TOLERANCE] = 0.0
        out[np.abs(out - 1.0) < _SNAP_TOLERANCE] = 1.0
        out[np.abs(out + 1.0) < _SNAP_TOLERANCE] = -1.0
        return out

    radians = np.deg2rad(angles)
    c1, s1 = snapped(np.cos(radians[:, 0])), snapped(np.sin(radians[:, 0]))
    c2, s2 = snapped(np.cos(radians[:, 1])), snapped(np.sin(radians[:, 1]))
    c3, s3 = snapped(np.cos(radians[:, 2])), snapped(np.sin(radians[:, 2]))

    coefficients = np.empty((angles.shape[0], 9), dtype=float)
    coefficients[:, 0] = c1 * c3 - s1 * s2 * s3
    coefficients[:, 1] = -s1 * c2
    coefficients[:, 2] = c1 * s3 + s1 * s2 * c3
    coefficients[:, 3] = s1 * c3 + c1 * s2 * s3
    coefficients[:, 4] = c1 * c2
    coefficients[:, 5] = s1 * s3 - c1 * s2 * c3
    coefficients[:, 6] = -c2 * s3
    coefficients[:, 7] = s2
    coefficients[:, 8] = c2 * c3
    return coefficients


class CoordinateSystem(ABC):
    """A coordinate system: an origin, a rotation, and the coordinates of the type on the rotated frame.

    Subclasses define :attr:`type`, the default orientation (:attr:`DEFAULT_NXYPL`) and the conversion between the
    local Cartesian frame and their own coordinates and vector components.
    """

    #: Name of the type in the model data ("cartesian", "cylindrical", "spherical").
    TYPE: ClassVar[str]
    #: Point of the local xy plane used when no orientation is given (with ``nxax = (1, 0, 0)``).
    DEFAULT_NXYPL: ClassVar[tuple[float, float, float]] = (0.0, 0.0, -1.0)

    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        nxax: Sequence[float] | None = None,
        nxypl: Sequence[float] | None = None,
        thxy_deg: float | None = None,
        thyz_deg: float | None = None,
        thzx_deg: float | None = None,
        local_to_global_cart: Any = None,
    ) -> None:
        """Create a coordinate system.

        Args:
            origin: Origin, in global coordinates.
            nxax: A point of the local x axis (with ``nxypl``).
            nxypl: A point of the local xy plane (with ``nxax``).
            thxy_deg: Rotation angles in degrees, all three or none.
            thyz_deg: See ``thxy_deg``.
            thzx_deg: See ``thxy_deg``.
            local_to_global_cart: The rotation matrix, columns = local axes in global coordinates.

        Raises:
            ValueError: If the orientation is given more than one way, only some of the angles are given, or the
                values are invalid (not finite, collinear points, non-orthonormal matrix).
        """
        self._origin = _vector3(origin, "origin")
        self._rotation = self._resolve_rotation(nxax, nxypl, (thxy_deg, thyz_deg, thzx_deg), local_to_global_cart)

    def _resolve_rotation(
        self,
        nxax: Sequence[float] | None,
        nxypl: Sequence[float] | None,
        angles: tuple[float | None, float | None, float | None],
        matrix: Any,
    ) -> np.ndarray:
        given_angles = [angle is not None for angle in angles]
        if any(given_angles) and not all(given_angles):
            raise ValueError("thxy_deg, thyz_deg and thzx_deg must be given together.")
        if matrix is not None:
            if nxax is not None or nxypl is not None or any(given_angles):
                raise ValueError("local_to_global_cart cannot be combined with nxax / nxypl or angles.")
            return _rotation_matrix(matrix, "local_to_global_cart")
        if all(given_angles):
            if nxax is not None or nxypl is not None:
                raise ValueError("The angles cannot be combined with nxax / nxypl.")
            return CoordinateSystem.build_local_to_global_matrix(*angles)  # type: ignore[arg-type]
        x_point = np.asarray((1.0, 0.0, 0.0) if nxax is None else nxax, dtype=float)
        plane_point = np.asarray(self.DEFAULT_NXYPL if nxypl is None else nxypl, dtype=float)
        x_axis = _unit(x_point - self._origin, "The local x axis is degenerate (nxax at the origin).")
        plane = _unit(plane_point - self._origin, "The local xy plane is degenerate (nxypl at the origin).")
        z_axis = _unit(np.cross(x_axis, plane), "nxax and nxypl are collinear with the origin.")
        y_axis = _unit(np.cross(z_axis, x_axis), "The local y axis is degenerate.")
        return np.column_stack((x_axis, y_axis, z_axis))

    # ------------------------------------------------------------------------------------------- constructors

    @classmethod
    def from_angles(cls, origin: Sequence[float], angles_deg: Sequence[float]) -> Self:
        """A coordinate system of this type from its origin and its angles ``(thxy, thyz, thzx)`` in degrees."""
        if len(angles_deg) != 3:
            raise ValueError("angles_deg must have 3 values (thxy, thyz, thzx).")
        thxy, thyz, thzx = angles_deg
        return cls(origin, thxy_deg=thxy, thyz_deg=thyz, thzx_deg=thzx)

    @classmethod
    def from_points(cls, origin: Sequence[float], nxax: Sequence[float], nxypl: Sequence[float]) -> Self:
        """A coordinate system of this type from its origin, a point of its x axis and a point of its xy plane."""
        return cls(origin, nxax=nxax, nxypl=nxypl)

    @staticmethod
    def build_local_to_global_matrix(angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float) -> np.ndarray:
        """The local-to-global rotation matrix of three angles in degrees."""
        return rotation_coefficients_from_angles([[angle_xy_deg, angle_yz_deg, angle_xz_deg]])[0].reshape(3, 3)

    @staticmethod
    def build_global_to_local_matrix(angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float) -> np.ndarray:
        """The global-to-local rotation matrix of three angles in degrees (the transpose of the local-to-global one)."""
        return CoordinateSystem.build_local_to_global_matrix(angle_xy_deg, angle_yz_deg, angle_xz_deg).T

    @staticmethod
    def build_local_to_global_coefficients(angles_deg: Any) -> np.ndarray:
        """Flattened local-to-global matrices of an (N, 3) array of angles in degrees.

        See :func:`rotation_coefficients_from_angles`; called by the C++ solver for the ``node_angles`` point array.
        """
        return rotation_coefficients_from_angles(angles_deg)

    # ------------------------------------------------------------------------------------------- attributes

    @property
    def type(self) -> str:
        """Name of the type in the model data: ``"cartesian"``, ``"cylindrical"`` or ``"spherical"``."""
        return self.TYPE

    @property
    def origin(self) -> np.ndarray:
        """Origin in global coordinates (a copy)."""
        return self._origin.copy()

    @property
    def local_to_global_cart(self) -> np.ndarray:
        """Rotation matrix whose columns are the local x, y, z axes in global coordinates (a copy)."""
        return self._rotation.copy()

    def axis_direction(self) -> np.ndarray:
        """Global direction of the local z axis (the axis of a cylindrical system)."""
        return self._rotation[:, 2].copy()

    def __repr__(self) -> str:
        return f"{type(self).__name__}(origin={self._origin.tolist()}, local_to_global_cart={self._rotation.tolist()})"

    # ------------------------------------------------------------------------------------------- points

    def _local_cartesian(self, global_points: Any) -> np.ndarray:
        return (_points3(global_points, "global_points") - self._origin) @ self._rotation

    def to_local(self, global_points: Any) -> np.ndarray:
        """Coordinates of the type (x y z, r theta z, or r theta phi; angles in degrees) of global points.

        Args:
            global_points: Shape (N, 3) or (3,).

        Returns:
            Shape (N, 3).
        """
        return self._from_cartesian(self._local_cartesian(global_points))

    def to_global(self, local_points: Any) -> np.ndarray:
        """Global points of coordinates of the type; the inverse of :meth:`to_local`.

        Args:
            local_points: Shape (N, 3) or (3,).

        Returns:
            Shape (N, 3).
        """
        return self._origin + self._to_cartesian(_points3(local_points, "local_points")) @ self._rotation.T

    # ------------------------------------------------------------------------------------------- vectors

    def vector_components_to_local(self, global_vectors: Any, global_positions: Any) -> np.ndarray:
        """Components of global Cartesian vectors in the local basis of the type at their positions.

        Args:
            global_vectors: Vectors, shape (N, 3) or (3,).
            global_positions: Where each vector is, shape (N, 3) or (3,).

        Returns:
            Local components (x y z, r theta z, or r theta phi), shape (N, 3).

        Raises:
            ValueError: If the shapes are invalid or the row counts differ.
        """
        vectors, positions = self._check_pairs(global_vectors, global_positions, "global_vectors")
        return self._components_to_local(vectors @ self._rotation, self._local_cartesian(positions))

    def vector_components_to_global(self, local_components: Any, global_positions: Any) -> np.ndarray:
        """Global Cartesian vectors of components in the local basis of the type; inverse of the above.

        Args:
            local_components: Components, shape (N, 3) or (3,).
            global_positions: Where each vector is, shape (N, 3) or (3,).

        Returns:
            Global vectors, shape (N, 3).
        """
        components, positions = self._check_pairs(local_components, global_positions, "local_components")
        return self._components_to_cartesian(components, self._local_cartesian(positions)) @ self._rotation.T

    def rotate_vector_to_global(self, local_components: Any, global_positions: Any) -> np.ndarray:
        """Same as :meth:`vector_components_to_global`."""
        return self.vector_components_to_global(local_components, global_positions)

    def local_to_global_matrices(self, global_positions: Any) -> np.ndarray:
        """Rotation from the local basis of the type to the global one at each position.

        Column ``k`` of each matrix is the k-th local unit vector (e.g. e_r, e_theta, e_z) in global coordinates.
        Called by the C++ solver for the elements of a COORDINATE_SYSTEM.

        Args:
            global_positions: Shape (N, 3) or (3,).

        Returns:
            Shape (N, 3, 3).
        """
        positions = _points3(global_positions, "global_positions")
        matrices = np.zeros((positions.shape[0], 3, 3), dtype=float)
        for axis in range(3):
            unit = np.zeros_like(positions)
            unit[:, axis] = 1.0
            matrices[:, :, axis] = self.vector_components_to_global(unit, positions)
        return matrices

    def local_to_global_matrix(self, global_position: Any) -> np.ndarray:
        """:meth:`local_to_global_matrices` at one position, shape (3, 3)."""
        return self.local_to_global_matrices(global_position)[0]

    def rotate_tensor_to_global(self, local_tensors: Any, global_positions: Any) -> np.ndarray:
        """Rotate tensors from the local basis of the type to the global one: ``Q T Q^T`` at each position.

        Args:
            local_tensors: Shape (3, 3), or (N, 3, 3); a single tensor is used at every position.
            global_positions: Shape (N, 3) or (3,).

        Returns:
            Shape (3, 3) for a single tensor input, (N, 3, 3) otherwise.
        """
        positions = _points3(global_positions, "global_positions")
        matrices = self.local_to_global_matrices(positions)
        tensors = np.asarray(local_tensors, dtype=float)
        single = tensors.shape == (3, 3)
        if single:
            tensors = tensors[np.newaxis]
        if tensors.ndim != 3 or tensors.shape[1:] != (3, 3):
            raise ValueError(f"local_tensors must have shape (3, 3) or (N, 3, 3), got {tensors.shape}.")
        if tensors.shape[0] != positions.shape[0]:
            if tensors.shape[0] != 1:
                raise ValueError("local_tensors and global_positions must have the same number of rows.")
            tensors = np.repeat(tensors, positions.shape[0], axis=0)
        rotated = np.einsum("nij,njk,nlk->nil", matrices, tensors, matrices)
        return rotated[0] if single else rotated

    @staticmethod
    def _check_pairs(vectors: Any, positions: Any, name: str) -> tuple[np.ndarray, np.ndarray]:
        vectors_array = _points3(vectors, name)
        positions_array = _points3(positions, "global_positions")
        if vectors_array.shape[0] != positions_array.shape[0]:
            raise ValueError(f"{name} and global_positions must have the same number of rows.")
        return vectors_array, positions_array

    # ------------------------------------------------------------------------------------------- per type

    @abstractmethod
    def _from_cartesian(self, local: np.ndarray) -> np.ndarray:
        """Coordinates of the type of local Cartesian points."""

    @abstractmethod
    def _to_cartesian(self, coordinates: np.ndarray) -> np.ndarray:
        """Local Cartesian points of coordinates of the type."""

    @abstractmethod
    def _components_to_local(self, local_vectors: np.ndarray, local: np.ndarray) -> np.ndarray:
        """Components in the basis of the type of vectors given in the local Cartesian basis, at points ``local``."""

    @abstractmethod
    def _components_to_cartesian(self, components: np.ndarray, local: np.ndarray) -> np.ndarray:
        """Local Cartesian vectors of components in the basis of the type, at points ``local``."""


class CartesianCoordinateSystem(CoordinateSystem):
    """Cartesian coordinates (x, y, z) on the rotated frame. By default the global axes."""

    TYPE = "cartesian"
    DEFAULT_NXYPL = (0.0, 1.0, 0.0)

    def _from_cartesian(self, local: np.ndarray) -> np.ndarray:
        return local

    def _to_cartesian(self, coordinates: np.ndarray) -> np.ndarray:
        return coordinates

    def _components_to_local(self, local_vectors: np.ndarray, local: np.ndarray) -> np.ndarray:
        return local_vectors

    def _components_to_cartesian(self, components: np.ndarray, local: np.ndarray) -> np.ndarray:
        return components


class CylindricalCoordinateSystem(CoordinateSystem):
    """Cylindrical coordinates (r, theta in degrees, z) around the local z axis.

    By default the axis is the global -y axis (local xy plane through (0, 0, -1)). On the axis (r = 0) the basis is
    the local Cartesian one.
    """

    TYPE = "cylindrical"

    def _from_cartesian(self, local: np.ndarray) -> np.ndarray:
        x, y = local[:, 0], local[:, 1]
        r = np.sqrt(x * x + y * y)
        theta = np.zeros_like(r)
        nonzero = r > 0.0
        theta[nonzero] = np.arctan2(y[nonzero], x[nonzero])
        return np.stack((r, np.rad2deg(theta), local[:, 2]), axis=1)

    def _to_cartesian(self, coordinates: np.ndarray) -> np.ndarray:
        r = coordinates[:, 0]
        theta = np.deg2rad(coordinates[:, 1])
        return np.stack((r * np.cos(theta), r * np.sin(theta), coordinates[:, 2]), axis=1)

    def _components_to_local(self, local_vectors: np.ndarray, local: np.ndarray) -> np.ndarray:
        x, y = local[:, 0], local[:, 1]
        r = np.sqrt(x * x + y * y)
        components = local_vectors.copy()
        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            vx = local_vectors[nonzero, 0]
            vy = local_vectors[nonzero, 1]
            components[nonzero, 0] = (x[nonzero] * vx + y[nonzero] * vy) * inv_r
            components[nonzero, 1] = (-y[nonzero] * vx + x[nonzero] * vy) * inv_r
        return components

    def _components_to_cartesian(self, components: np.ndarray, local: np.ndarray) -> np.ndarray:
        x, y = local[:, 0], local[:, 1]
        r = np.sqrt(x * x + y * y)
        comp_r, comp_theta = components[:, 0], components[:, 1]
        vectors = components.copy()
        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            xr = x[nonzero] * inv_r
            yr = y[nonzero] * inv_r
            vectors[nonzero, 0] = xr * comp_r[nonzero] - yr * comp_theta[nonzero]
            vectors[nonzero, 1] = yr * comp_r[nonzero] + xr * comp_theta[nonzero]
        return vectors


class SphericalCoordinateSystem(CoordinateSystem):
    """Spherical coordinates (r, theta, phi in degrees): theta from the local z axis, phi around it.

    At the origin the basis is (e_z, e_x, e_y); on the z axis e_r is radial and (e_theta, e_phi) = (e_x, e_y).
    """

    TYPE = "spherical"

    @staticmethod
    def _basis(local: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Unit vectors e_r, e_theta, e_phi at each point, in local Cartesian components."""
        x, y, z = local[:, 0], local[:, 1], local[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        rho = np.sqrt(x * x + y * y)

        e_r = np.zeros_like(local)
        e_theta = np.zeros_like(local)
        e_phi = np.zeros_like(local)
        e_r[:, 2] = 1.0
        e_theta[:, 0] = 1.0
        e_phi[:, 1] = 1.0

        nonzero_r = r > 0.0
        if np.any(nonzero_r):
            inv_r = 1.0 / r[nonzero_r]
            e_r[nonzero_r, 0] = x[nonzero_r] * inv_r
            e_r[nonzero_r, 1] = y[nonzero_r] * inv_r
            e_r[nonzero_r, 2] = z[nonzero_r] * inv_r

            nonzero_rho = nonzero_r & (rho > 0.0)
            if np.any(nonzero_rho):
                inv_r_all = np.zeros_like(r)
                inv_r_all[nonzero_r] = inv_r
                inv_rho = 1.0 / rho[nonzero_rho]
                cos_theta = z[nonzero_rho] * inv_r_all[nonzero_rho]
                sin_theta = rho[nonzero_rho] * inv_r_all[nonzero_rho]
                cos_phi = x[nonzero_rho] * inv_rho
                sin_phi = y[nonzero_rho] * inv_rho
                e_theta[nonzero_rho, 0] = cos_theta * cos_phi
                e_theta[nonzero_rho, 1] = cos_theta * sin_phi
                e_theta[nonzero_rho, 2] = -sin_theta
                e_phi[nonzero_rho, 0] = -sin_phi
                e_phi[nonzero_rho, 1] = cos_phi
                e_phi[nonzero_rho, 2] = 0.0
        return e_r, e_theta, e_phi

    def _from_cartesian(self, local: np.ndarray) -> np.ndarray:
        x, y, z = local[:, 0], local[:, 1], local[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        theta = np.zeros_like(r)
        phi = np.zeros_like(r)
        nonzero = r > 0.0
        theta[nonzero] = np.arccos(z[nonzero] / r[nonzero])
        phi[nonzero] = np.arctan2(y[nonzero], x[nonzero])
        return np.stack((r, np.rad2deg(theta), np.rad2deg(phi)), axis=1)

    def _to_cartesian(self, coordinates: np.ndarray) -> np.ndarray:
        r = coordinates[:, 0]
        theta = np.deg2rad(coordinates[:, 1])
        phi = np.deg2rad(coordinates[:, 2])
        sin_theta = np.sin(theta)
        return np.stack((r * sin_theta * np.cos(phi), r * sin_theta * np.sin(phi), r * np.cos(theta)), axis=1)

    def _components_to_local(self, local_vectors: np.ndarray, local: np.ndarray) -> np.ndarray:
        e_r, e_theta, e_phi = self._basis(local)
        return np.stack(
            (
                np.einsum("ij,ij->i", e_r, local_vectors),
                np.einsum("ij,ij->i", e_theta, local_vectors),
                np.einsum("ij,ij->i", e_phi, local_vectors),
            ),
            axis=1,
        )

    def _components_to_cartesian(self, components: np.ndarray, local: np.ndarray) -> np.ndarray:
        e_r, e_theta, e_phi = self._basis(local)
        return (
            e_r * components[:, 0, np.newaxis]
            + e_theta * components[:, 1, np.newaxis]
            + e_phi * components[:, 2, np.newaxis]
        )


class _AxisCylindricalCoordinateSystem(CylindricalCoordinateSystem):
    """A cylindrical system with a fixed orientation: only the origin (or an explicit matrix) can be given."""

    ANGLES: ClassVar[tuple[float, float, float]]

    def __init__(self, origin: Sequence[float] = (0.0, 0.0, 0.0), local_to_global_cart: Any = None) -> None:
        if local_to_global_cart is not None:
            super().__init__(origin=origin, local_to_global_cart=local_to_global_cart)
        else:
            thxy, thyz, thzx = self.ANGLES
            super().__init__(origin=origin, thxy_deg=thxy, thyz_deg=thyz, thzx_deg=thzx)


class XAxisCylindricalCoordinateSystem(_AxisCylindricalCoordinateSystem):
    """Cylindrical coordinates around the global x axis."""

    ANGLES = X_AXIS_CYLINDRICAL_ANGLES


class YAxisCylindricalCoordinateSystem(_AxisCylindricalCoordinateSystem):
    """Cylindrical coordinates around the global y axis (the axis of a 2D axisymmetric model)."""

    ANGLES = Y_AXIS_CYLINDRICAL_ANGLES


class ZAxisCylindricalCoordinateSystem(_AxisCylindricalCoordinateSystem):
    """Cylindrical coordinates around the global z axis."""

    ANGLES = Z_AXIS_CYLINDRICAL_ANGLES


#: The ``type`` names of the model data and their classes.
COORDINATE_SYSTEM_TYPES: dict[str, type[CoordinateSystem]] = {
    "cartesian": CartesianCoordinateSystem,
    "cylindrical": CylindricalCoordinateSystem,
    "x_axis_cylindrical": XAxisCylindricalCoordinateSystem,
    "y_axis_cylindrical": YAxisCylindricalCoordinateSystem,
    "z_axis_cylindrical": ZAxisCylindricalCoordinateSystem,
    "spherical": SphericalCoordinateSystem,
}


def coordinate_system_from_dict(spec: Mapping[str, Any]) -> CoordinateSystem:
    """Build a coordinate system from its attributes in the model data.

    Attributes of ``spec``:

    * ``type``: a key of :data:`COORDINATE_SYSTEM_TYPES` (required);
    * ``origin``: default (0, 0, 0);
    * at most one of ``matrix`` (3x3 local-to-global rotation), ``angles`` (``[thxy, thyz, thzx]`` in degrees) and
      ``points`` (``{"nxax": [...], "nxypl": [...]}``); none gives the default orientation of the type. The x / y /
      z axis cylindrical types only accept ``matrix``.

    Raises:
        ValueError: On an unknown type or attribute, or an invalid combination or value.
    """
    if not isinstance(spec, Mapping):
        raise ValueError(f"A coordinate system must be an object, got {type(spec).__name__}.")
    unknown = set(spec) - {"type", "origin", "matrix", "angles", "points"}
    if unknown:
        raise ValueError(f"Unknown coordinate system attributes: {sorted(unknown)}.")

    type_name = spec.get("type")
    cls = COORDINATE_SYSTEM_TYPES.get(type_name) if isinstance(type_name, str) else None
    if cls is None:
        raise ValueError(f"type must be one of {sorted(COORDINATE_SYSTEM_TYPES)}, got {type_name!r}.")

    orientations = [key for key in ("matrix", "angles", "points") if spec.get(key) is not None]
    if len(orientations) > 1:
        raise ValueError(f"Give one of matrix, angles or points, not {orientations}.")
    if issubclass(cls, _AxisCylindricalCoordinateSystem) and orientations not in ([], ["matrix"]):
        raise ValueError(f"The type {type_name!r} has a fixed orientation: only 'matrix' is accepted.")

    kwargs: dict[str, Any] = {}
    if "origin" in spec:
        kwargs["origin"] = spec["origin"]
    if orientations == ["matrix"]:
        kwargs["local_to_global_cart"] = np.asarray(spec["matrix"], dtype=float)
    elif orientations == ["angles"]:
        angles = spec["angles"]
        if len(angles) != 3:
            raise ValueError("angles must be [thxy, thyz, thzx] in degrees.")
        kwargs.update(thxy_deg=angles[0], thyz_deg=angles[1], thzx_deg=angles[2])
    elif orientations == ["points"]:
        points = spec["points"]
        if not isinstance(points, Mapping) or set(points) != {"nxax", "nxypl"}:
            raise ValueError("points must be an object with 'nxax' and 'nxypl'.")
        kwargs.update(nxax=points["nxax"], nxypl=points["nxypl"])
    return cls(**kwargs)


def coordinate_system_to_dict(coordinate_system: CoordinateSystem) -> dict[str, Any]:
    """The model data attributes of a coordinate system: ``{"type", "origin", "matrix"}``."""
    return {
        "type": coordinate_system.type,
        "origin": coordinate_system.origin.tolist(),
        "matrix": coordinate_system.local_to_global_cart.tolist(),
    }

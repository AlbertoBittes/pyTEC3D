import math
from typing import Optional, Sequence

import numpy as np


def _normalize_or_error(vec: np.ndarray, error_message: str) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm <= 0.0:
        raise ValueError(error_message)
    return vec / norm


def _build_rotation_from_basis(e0: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    rotation = np.zeros((3, 3), dtype=float)
    rotation[0, :] = e0
    rotation[1, :] = e1
    rotation[2, :] = e2
    return rotation


def _build_local_to_global_from_axes(
    x_axis: np.ndarray, y_axis: np.ndarray, z_axis: np.ndarray
) -> np.ndarray:
    rotation = np.zeros((3, 3), dtype=float)
    rotation[:, 0] = x_axis
    rotation[:, 1] = y_axis
    rotation[:, 2] = z_axis
    return rotation


def _build_local_to_global_from_points(
    origin: np.ndarray, nxax: np.ndarray, nxypl: np.ndarray
) -> np.ndarray:
    origin_vec = np.asarray(origin, dtype=float)
    x_axis = _normalize_or_error(nxax - origin_vec, "Coordinate system x-axis is invalid.")
    plane_vec = _normalize_or_error(nxypl - origin_vec, "Coordinate system xy-plane vector is invalid.")
    z_axis = _normalize_or_error(
        np.cross(x_axis, plane_vec), "Coordinate system axes are collinear."
    )
    y_axis = _normalize_or_error(np.cross(z_axis, x_axis), "Coordinate system y-axis is invalid.")
    return _build_local_to_global_from_axes(x_axis, y_axis, z_axis)


def build_ansys_rotation_matrix(
    angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float
) -> np.ndarray:
    def snap_trig_value(value: float) -> float:
        eps = 1e-12
        if abs(value) < eps:
            return 0.0
        if abs(value - 1.0) < eps:
            return 1.0
        if abs(value + 1.0) < eps:
            return -1.0
        return value

    to_rad = math.pi / 180.0
    angle_xy = angle_xy_deg * to_rad
    angle_yz = angle_yz_deg * to_rad
    angle_xz = angle_xz_deg * to_rad

    c1 = snap_trig_value(math.cos(angle_xy))
    s1 = snap_trig_value(math.sin(angle_xy))
    c2 = snap_trig_value(math.cos(angle_yz))
    s2 = snap_trig_value(math.sin(angle_yz))
    c3 = snap_trig_value(math.cos(angle_xz))
    s3 = snap_trig_value(math.sin(angle_xz))

    rotation = np.zeros((3, 3), dtype=float)
    rotation[0, 0] = c1 * c3 - s1 * s2 * s3
    rotation[0, 1] = -s1 * c2
    rotation[0, 2] = c1 * s3 + s1 * s2 * c3
    rotation[1, 0] = s1 * c3 + c1 * s2 * s3
    rotation[1, 1] = c1 * c2
    rotation[1, 2] = s1 * s3 - c1 * s2 * c3
    rotation[2, 0] = -c2 * s3
    rotation[2, 1] = s2
    rotation[2, 2] = c2 * c3
    return rotation


class CoordinateSystem:
    def __init__(self, origin: Sequence[float], local_to_global_cart: np.ndarray, cs_type: str):
        self._origin = np.asarray(origin, dtype=float).reshape(3)
        self._local_to_global_cart = np.asarray(local_to_global_cart, dtype=float).reshape(3, 3)
        self._type = cs_type

    @property
    def type(self) -> str:
        return self._type

    @property
    def origin(self) -> np.ndarray:
        return self._origin

    def axis_direction(self) -> np.ndarray:
        return self._local_to_global_cart[:, 2].copy()

    def _to_local_vector(self, global_point: np.ndarray) -> np.ndarray:
        rel = np.asarray(global_point, dtype=float) - self._origin
        return rel @ self._local_to_global_cart

    def _to_local_vectors(self, global_points: np.ndarray) -> np.ndarray:
        rel = np.asarray(global_points, dtype=float) - self._origin
        return rel @ self._local_to_global_cart

    def to_local_coordinates(self, global_point: Sequence[float]) -> np.ndarray:
        raise NotImplementedError

    def to_local_coordinates_array(self, global_points: np.ndarray) -> np.ndarray:
        return np.vstack([self.to_local_coordinates(point) for point in global_points])

    def global_to_local_rotation(self, global_point: Sequence[float]) -> np.ndarray:
        raise NotImplementedError

    def global_to_local_rotation_array(self, global_points: np.ndarray) -> np.ndarray:
        return np.stack([self.global_to_local_rotation(point) for point in global_points], axis=0)

    def local_to_global_rotation(self, global_point: Sequence[float]) -> np.ndarray:
        return self.global_to_local_rotation(global_point).T

    def local_to_global_rotation_array(self, global_points: np.ndarray) -> np.ndarray:
        rotations = self.global_to_local_rotation_array(global_points)
        return np.transpose(rotations, axes=(0, 2, 1))


class CartesianCoordinateSystem(CoordinateSystem):
    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        nxax: Optional[Sequence[float]] = None,
        nxypl: Optional[Sequence[float]] = None,
        thxy_deg: Optional[float] = None,
        thyz_deg: Optional[float] = None,
        thzx_deg: Optional[float] = None,
    ):
        if thxy_deg is not None and thyz_deg is not None and thzx_deg is not None:
            rotation = build_ansys_rotation_matrix(thxy_deg, thyz_deg, thzx_deg)
        else:
            if nxax is None:
                nxax = (1.0, 0.0, 0.0)
            if nxypl is None:
                nxypl = (0.0, 1.0, 0.0)
            rotation = _build_local_to_global_from_points(
                np.asarray(origin, dtype=float),
                np.asarray(nxax, dtype=float),
                np.asarray(nxypl, dtype=float),
            )
        super().__init__(origin, rotation, "cartesian")

    def global_to_local_rotation(self, global_point: Sequence[float]) -> np.ndarray:
        return self._local_to_global_cart.T.copy()

    def to_local_coordinates(self, global_point: Sequence[float]) -> np.ndarray:
        local = self._to_local_vector(global_point)
        return np.asarray(local, dtype=float)


class CylindricalCoordinateSystem(CoordinateSystem):
    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        nxax: Optional[Sequence[float]] = None,
        nxypl: Optional[Sequence[float]] = None,
        thxy_deg: Optional[float] = None,
        thyz_deg: Optional[float] = None,
        thzx_deg: Optional[float] = None,
    ):
        if thxy_deg is not None and thyz_deg is not None and thzx_deg is not None:
            rotation = build_ansys_rotation_matrix(thxy_deg, thyz_deg, thzx_deg)
        else:
            if nxax is None:
                nxax = (1.0, 0.0, 0.0)
            if nxypl is None:
                nxypl = (0.0, 1.0, 0.0)
            rotation = _build_local_to_global_from_points(
                np.asarray(origin, dtype=float),
                np.asarray(nxax, dtype=float),
                np.asarray(nxypl, dtype=float),
            )
        super().__init__(origin, rotation, "cylindrical")

    def global_to_local_rotation(self, global_point: Sequence[float]) -> np.ndarray:
        local = self._to_local_vector(global_point)
        x = float(local[0])
        z = float(local[2])
        r = math.sqrt(x * x + z * z)

        e_r_local = np.array([1.0, 0.0, 0.0], dtype=float)
        e_theta_local = np.array([0.0, 0.0, -1.0], dtype=float)
        if r > 0.0:
            inv_r = 1.0 / r
            e_r_local = np.array([x * inv_r, 0.0, z * inv_r], dtype=float)
            e_theta_local = np.array([z * inv_r, 0.0, -x * inv_r], dtype=float)
        e_y_local = np.array([0.0, 1.0, 0.0], dtype=float)
        e_z_local = -e_theta_local

        e_r_global = e_r_local @ self._local_to_global_cart.T
        e_y_global = e_y_local @ self._local_to_global_cart.T
        e_z_global = e_z_local @ self._local_to_global_cart.T
        return _build_rotation_from_basis(e_r_global, e_y_global, e_z_global)

    def global_to_local_rotation_array(self, global_points: np.ndarray) -> np.ndarray:
        local = self._to_local_vectors(global_points)
        x = local[:, 0]
        z = local[:, 2]
        r = np.sqrt(x * x + z * z)

        e_r_local = np.zeros_like(local)
        e_theta_local = np.zeros_like(local)
        e_r_local[:, 0] = 1.0
        e_theta_local[:, 2] = -1.0

        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            e_r_local[nonzero, 0] = x[nonzero] * inv_r
            e_r_local[nonzero, 2] = z[nonzero] * inv_r
            e_theta_local[nonzero, 0] = z[nonzero] * inv_r
            e_theta_local[nonzero, 2] = -x[nonzero] * inv_r

        e_y_local = np.zeros_like(local)
        e_y_local[:, 1] = 1.0
        e_z_local = -e_theta_local

        e_r_global = e_r_local @ self._local_to_global_cart.T
        e_y_global = e_y_local @ self._local_to_global_cart.T
        e_z_global = e_z_local @ self._local_to_global_cart.T

        rotations = np.zeros((local.shape[0], 3, 3), dtype=float)
        rotations[:, 0, :] = e_r_global
        rotations[:, 1, :] = e_y_global
        rotations[:, 2, :] = e_z_global
        return rotations

    def to_local_coordinates(self, global_point: Sequence[float]) -> np.ndarray:
        local = self._to_local_vector(global_point)
        x = float(local[0])
        z = float(local[2])
        r = math.sqrt(x * x + z * z)
        theta = math.atan2(-z, x) if r > 0.0 else 0.0
        return np.array([r, theta, float(local[1])], dtype=float)

    def to_local_coordinates_array(self, global_points: np.ndarray) -> np.ndarray:
        local = self._to_local_vectors(global_points)
        x = local[:, 0]
        z = local[:, 2]
        r = np.sqrt(x * x + z * z)
        theta = np.zeros_like(r)
        nonzero = r > 0.0
        theta[nonzero] = np.arctan2(-z[nonzero], x[nonzero])
        return np.stack((r, theta, local[:, 1]), axis=1)


class SphericalCoordinateSystem(CoordinateSystem):
    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        nxax: Optional[Sequence[float]] = None,
        nxypl: Optional[Sequence[float]] = None,
        thxy_deg: Optional[float] = None,
        thyz_deg: Optional[float] = None,
        thzx_deg: Optional[float] = None,
    ):
        if thxy_deg is not None and thyz_deg is not None and thzx_deg is not None:
            rotation = build_ansys_rotation_matrix(thxy_deg, thyz_deg, thzx_deg)
        else:
            if nxax is None:
                nxax = (1.0, 0.0, 0.0)
            if nxypl is None:
                nxypl = (0.0, 1.0, 0.0)
            rotation = _build_local_to_global_from_points(
                np.asarray(origin, dtype=float),
                np.asarray(nxax, dtype=float),
                np.asarray(nxypl, dtype=float),
            )
        super().__init__(origin, rotation, "spherical")

    def global_to_local_rotation(self, global_point: Sequence[float]) -> np.ndarray:
        local = self._to_local_vector(global_point)
        x = float(local[0])
        y = float(local[1])
        z = float(local[2])
        r = math.sqrt(x * x + y * y + z * z)
        rho = math.sqrt(x * x + y * y)

        e_r_local = np.array([0.0, 0.0, 1.0], dtype=float)
        e_theta_local = np.array([1.0, 0.0, 0.0], dtype=float)
        e_phi_local = np.array([0.0, 1.0, 0.0], dtype=float)

        if r > 0.0:
            inv_r = 1.0 / r
            e_r_local = np.array([x * inv_r, y * inv_r, z * inv_r], dtype=float)
            if rho > 0.0:
                inv_rho = 1.0 / rho
                cos_theta = z * inv_r
                sin_theta = rho * inv_r
                cos_phi = x * inv_rho
                sin_phi = y * inv_rho
                e_theta_local = np.array(
                    [cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta], dtype=float
                )
                e_phi_local = np.array([-sin_phi, cos_phi, 0.0], dtype=float)

        e_r_global = e_r_local @ self._local_to_global_cart.T
        e_theta_global = e_theta_local @ self._local_to_global_cart.T
        e_phi_global = e_phi_local @ self._local_to_global_cart.T
        return _build_rotation_from_basis(e_r_global, e_theta_global, e_phi_global)

    def global_to_local_rotation_array(self, global_points: np.ndarray) -> np.ndarray:
        local = self._to_local_vectors(global_points)
        x = local[:, 0]
        y = local[:, 1]
        z = local[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        rho = np.sqrt(x * x + y * y)

        e_r_local = np.zeros_like(local)
        e_theta_local = np.zeros_like(local)
        e_phi_local = np.zeros_like(local)
        e_r_local[:, 2] = 1.0
        e_theta_local[:, 0] = 1.0
        e_phi_local[:, 1] = 1.0

        nonzero_r = r > 0.0
        if np.any(nonzero_r):
            inv_r = 1.0 / r[nonzero_r]
            e_r_local[nonzero_r, 0] = x[nonzero_r] * inv_r
            e_r_local[nonzero_r, 1] = y[nonzero_r] * inv_r
            e_r_local[nonzero_r, 2] = z[nonzero_r] * inv_r

            nonzero_rho = nonzero_r & (rho > 0.0)
            if np.any(nonzero_rho):
                inv_rho = 1.0 / rho[nonzero_rho]
                cos_theta = z[nonzero_rho] * inv_r[nonzero_rho]
                sin_theta = rho[nonzero_rho] * inv_r[nonzero_rho]
                cos_phi = x[nonzero_rho] * inv_rho
                sin_phi = y[nonzero_rho] * inv_rho
                e_theta_local[nonzero_rho, 0] = cos_theta * cos_phi
                e_theta_local[nonzero_rho, 1] = cos_theta * sin_phi
                e_theta_local[nonzero_rho, 2] = -sin_theta
                e_phi_local[nonzero_rho, 0] = -sin_phi
                e_phi_local[nonzero_rho, 1] = cos_phi
                e_phi_local[nonzero_rho, 2] = 0.0

        e_r_global = e_r_local @ self._local_to_global_cart.T
        e_theta_global = e_theta_local @ self._local_to_global_cart.T
        e_phi_global = e_phi_local @ self._local_to_global_cart.T

        rotations = np.zeros((local.shape[0], 3, 3), dtype=float)
        rotations[:, 0, :] = e_r_global
        rotations[:, 1, :] = e_theta_global
        rotations[:, 2, :] = e_phi_global
        return rotations

    def to_local_coordinates(self, global_point: Sequence[float]) -> np.ndarray:
        local = self._to_local_vector(global_point)
        x = float(local[0])
        y = float(local[1])
        z = float(local[2])
        r = math.sqrt(x * x + y * y + z * z)
        if r <= 0.0:
            return np.array([0.0, 0.0, 0.0], dtype=float)
        theta = math.acos(z / r)
        phi = math.atan2(y, x)
        return np.array([r, theta, phi], dtype=float)

    def to_local_coordinates_array(self, global_points: np.ndarray) -> np.ndarray:
        local = self._to_local_vectors(global_points)
        x = local[:, 0]
        y = local[:, 1]
        z = local[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        theta = np.zeros_like(r)
        phi = np.zeros_like(r)
        nonzero = r > 0.0
        theta[nonzero] = np.arccos(z[nonzero] / r[nonzero])
        phi[nonzero] = np.arctan2(y[nonzero], x[nonzero])
        return np.stack((r, theta, phi), axis=1)


def make_coordinate_system_local(
    kcs: int,
    origin: Sequence[float],
    thxy_deg: float,
    thyz_deg: float,
    thzx_deg: float,
) -> CoordinateSystem:
    if kcs == 0:
        return CartesianCoordinateSystem(origin, thxy_deg=thxy_deg, thyz_deg=thyz_deg, thzx_deg=thzx_deg)
    if kcs == 1:
        return CylindricalCoordinateSystem(origin, thxy_deg=thxy_deg, thyz_deg=thyz_deg, thzx_deg=thzx_deg)
    if kcs == 2:
        return SphericalCoordinateSystem(origin, thxy_deg=thxy_deg, thyz_deg=thyz_deg, thzx_deg=thzx_deg)
    raise ValueError("Unsupported coordinate system type.")


def make_coordinate_system_cs(
    kcs: int,
    origin: Sequence[float],
    nxax: Sequence[float],
    nxypl: Sequence[float],
) -> CoordinateSystem:
    if kcs == 0:
        return CartesianCoordinateSystem(origin, nxax=nxax, nxypl=nxypl)
    if kcs == 1:
        return CylindricalCoordinateSystem(origin, nxax=nxax, nxypl=nxypl)
    if kcs == 2:
        return SphericalCoordinateSystem(origin, nxax=nxax, nxypl=nxypl)
    raise ValueError("Unsupported coordinate system type.")

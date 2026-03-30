import math
from abc import ABC, abstractmethod
from typing import Any, Optional, Sequence, Type, TypeVar, cast

import numpy as np

TCoordinateSystem = TypeVar("TCoordinateSystem", bound="CoordinateSystem")

class CoordinateSystem(ABC):
    @abstractmethod
    def __init__(
        self, origin: Sequence[float], local_to_global_cart: np.ndarray, cs_type: str
    ):
        self._origin = np.asarray(origin, dtype=float).reshape(3)
        self._local_to_global_cart = np.asarray(
            local_to_global_cart, dtype=float
        ).reshape(3, 3)
        self._type = cs_type

    @property
    def type(self) -> str:
        return self._type

    @property
    def origin(self) -> np.ndarray:
        return self._origin

    def axis_direction(self) -> np.ndarray:
        return self._local_to_global_cart[:, 2].copy()

    def to_local(self, global_points: np.ndarray) -> np.ndarray:
        return (global_points - self._origin) @ self._local_to_global_cart

    def to_global(self, local_points: np.ndarray) -> np.ndarray:
        return self._origin + local_points @ self._local_to_global_cart.T

    def vector_components_to_local(
        self, global_vectors: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        return global_vectors @ self._local_to_global_cart

    def vector_components_to_global(
        self, local_components: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        return local_components @ self._local_to_global_cart.T

    @staticmethod
    def _deg_to_rad(values: np.ndarray) -> np.ndarray:
        return np.deg2rad(values)

    @staticmethod
    def _rad_to_deg(values: np.ndarray) -> np.ndarray:
        return np.rad2deg(values)

    @staticmethod
    def _normalize_or_error(vec: np.ndarray, error_message: str) -> np.ndarray:
        norm = np.linalg.norm(vec)
        if norm <= 0.0:
            raise ValueError(error_message)
        return vec / norm

    @staticmethod
    def _build_local_to_global_from_axes(
        x_axis: np.ndarray, y_axis: np.ndarray, z_axis: np.ndarray
    ) -> np.ndarray:
        rotation = np.zeros((3, 3), dtype=float)
        rotation[:, 0] = x_axis
        rotation[:, 1] = y_axis
        rotation[:, 2] = z_axis
        return rotation

    @staticmethod
    def _build_local_to_global_from_points(
        origin: np.ndarray, nxax: np.ndarray, nxypl: np.ndarray
    ) -> np.ndarray:
        origin_vec = np.asarray(origin, dtype=float)
        x_axis = CoordinateSystem._normalize_or_error(
            nxax - origin_vec, "Coordinate system x-axis is invalid."
        )
        plane_vec = CoordinateSystem._normalize_or_error(
            nxypl - origin_vec, "Coordinate system xy-plane vector is invalid."
        )
        z_axis = CoordinateSystem._normalize_or_error(
            np.cross(x_axis, plane_vec), "Coordinate system axes are collinear."
        )
        y_axis = CoordinateSystem._normalize_or_error(
            np.cross(z_axis, x_axis), "Coordinate system y-axis is invalid."
        )
        return CoordinateSystem._build_local_to_global_from_axes(x_axis, y_axis, z_axis)

    @staticmethod
    def build_local_to_global_matrix(
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

        angle_xy, angle_yz, angle_xz = CoordinateSystem._deg_to_rad(
            np.array([angle_xy_deg, angle_yz_deg, angle_xz_deg], dtype=float)
        )

        c1 = snap_trig_value(math.cos(angle_xy))
        s1 = snap_trig_value(math.sin(angle_xy))
        c2 = snap_trig_value(math.cos(angle_yz))
        s2 = snap_trig_value(math.sin(angle_yz))
        c3 = snap_trig_value(math.cos(angle_xz))
        s3 = snap_trig_value(math.sin(angle_xz))

        local_to_global = np.zeros((3, 3), dtype=float)
        local_to_global[0, 0] = c1 * c3 - s1 * s2 * s3
        local_to_global[0, 1] = -s1 * c2
        local_to_global[0, 2] = c1 * s3 + s1 * s2 * c3
        local_to_global[1, 0] = s1 * c3 + c1 * s2 * s3
        local_to_global[1, 1] = c1 * c2
        local_to_global[1, 2] = s1 * s3 - c1 * s2 * c3
        local_to_global[2, 0] = -c2 * s3
        local_to_global[2, 1] = s2
        local_to_global[2, 2] = c2 * c3
        return local_to_global

    @staticmethod
    def build_global_to_local_matrix(
        angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float
    ) -> np.ndarray:
        return CoordinateSystem.build_local_to_global_matrix(
            angle_xy_deg, angle_yz_deg, angle_xz_deg
        ).T

    @staticmethod
    def build_global_to_local(
        angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float
    ) -> np.ndarray:
        return CoordinateSystem.build_global_to_local_matrix(
            angle_xy_deg, angle_yz_deg, angle_xz_deg
        )

    @classmethod
    def from_angles(
        cls: Type[TCoordinateSystem],
        origin: Sequence[float],
        angles_deg: Sequence[float],
    ) -> TCoordinateSystem:
        thxy_deg, thyz_deg, thzx_deg = angles_deg
        ctor = cast(Any, cls)
        return ctor(origin, thxy_deg=thxy_deg, thyz_deg=thyz_deg, thzx_deg=thzx_deg)

    @classmethod
    def from_points(
        cls: Type[TCoordinateSystem],
        origin: Sequence[float],
        nxax: Sequence[float],
        nxypl: Sequence[float],
    ) -> TCoordinateSystem:
        ctor = cast(Any, cls)
        return ctor(origin, nxax=nxax, nxypl=nxypl)

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
            rotation = CoordinateSystem.build_local_to_global_matrix(
                thxy_deg, thyz_deg, thzx_deg
            )
        else:
            if nxax is None:
                nxax = (1.0, 0.0, 0.0)
            if nxypl is None:
                nxypl = (0.0, 1.0, 0.0)
            rotation = CoordinateSystem._build_local_to_global_from_points(
                np.asarray(origin, dtype=float),
                np.asarray(nxax, dtype=float),
                np.asarray(nxypl, dtype=float),
            )
        super().__init__(origin, rotation, "cartesian")

    def vector_components_to_local(
        self, global_vectors: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        return super().vector_components_to_local(global_vectors, global_positions)


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
            rotation = CoordinateSystem.build_local_to_global_matrix(
                thxy_deg, thyz_deg, thzx_deg
            )
        else:
            if nxax is None:
                nxax = (1.0, 0.0, 0.0)
            if nxypl is None:
                nxypl = (0.0, 0.0, -1.0)
            rotation = CoordinateSystem._build_local_to_global_from_points(
                np.asarray(origin, dtype=float),
                np.asarray(nxax, dtype=float),
                np.asarray(nxypl, dtype=float),
            )
        super().__init__(origin, rotation, "cylindrical")

    def vector_components_to_local(
        self, global_vectors: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        local_vectors = super().vector_components_to_local(
            global_vectors, global_positions
        )
        local_positions = super().to_local(global_positions)
        x = local_positions[:, 0]
        y = local_positions[:, 1]
        r = np.sqrt(x * x + y * y)

        e_r = np.zeros_like(local_positions)
        e_theta = np.zeros_like(local_positions)
        e_z = np.zeros_like(local_positions)
        e_r[:, 0] = 1.0
        e_theta[:, 1] = 1.0
        e_z[:, 2] = 1.0

        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            e_r[nonzero, 0] = x[nonzero] * inv_r
            e_r[nonzero, 1] = y[nonzero] * inv_r
            e_theta[nonzero, 0] = -y[nonzero] * inv_r
            e_theta[nonzero, 1] = x[nonzero] * inv_r

        comp_r = np.einsum("ij,ij->i", e_r, local_vectors)
        comp_theta = np.einsum("ij,ij->i", e_theta, local_vectors)
        comp_z = local_vectors[:, 2]
        return np.stack((comp_r, comp_theta, comp_z), axis=1)

    def vector_components_to_global(
        self, local_components: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        local_positions = super().to_local(global_positions)
        x = local_positions[:, 0]
        y = local_positions[:, 1]
        r = np.sqrt(x * x + y * y)

        e_r = np.zeros_like(local_positions)
        e_theta = np.zeros_like(local_positions)
        e_z = np.zeros_like(local_positions)
        e_r[:, 0] = 1.0
        e_theta[:, 1] = 1.0
        e_z[:, 2] = 1.0

        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            e_r[nonzero, 0] = x[nonzero] * inv_r
            e_r[nonzero, 1] = y[nonzero] * inv_r
            e_theta[nonzero, 0] = -y[nonzero] * inv_r
            e_theta[nonzero, 1] = x[nonzero] * inv_r

        local_components = np.asarray(local_components, dtype=float)
        comp_r = local_components[:, 0]
        comp_theta = local_components[:, 1]
        comp_z = local_components[:, 2]
        local_vectors = (
            e_r * comp_r[:, np.newaxis]
            + e_theta * comp_theta[:, np.newaxis]
            + e_z * comp_z[:, np.newaxis]
        )
        return local_vectors @ self._local_to_global_cart.T

    def to_local(self, global_points: np.ndarray) -> np.ndarray:
        local = super().to_local(global_points)
        x = local[:, 0]
        y = local[:, 1]
        r = np.sqrt(x * x + y * y)
        theta = np.zeros_like(r)
        nonzero = r > 0.0
        theta[nonzero] = np.arctan2(y[nonzero], x[nonzero])
        theta_deg = CoordinateSystem._rad_to_deg(theta)
        return np.stack((r, theta_deg, local[:, 2]), axis=1)

    def to_global(self, local_points: np.ndarray) -> np.ndarray:
        r = local_points[:, 0]
        theta = CoordinateSystem._deg_to_rad(local_points[:, 1])
        z = local_points[:, 2]
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        local_cart = np.stack((x, y, z), axis=1)
        return self._origin + local_cart @ self._local_to_global_cart.T


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
            rotation = CoordinateSystem.build_local_to_global_matrix(
                thxy_deg, thyz_deg, thzx_deg
            )
        else:
            if nxax is None:
                nxax = (1.0, 0.0, 0.0)
            if nxypl is None:
                nxypl = (0.0, 0.0, -1.0)
            rotation = CoordinateSystem._build_local_to_global_from_points(
                np.asarray(origin, dtype=float),
                np.asarray(nxax, dtype=float),
                np.asarray(nxypl, dtype=float),
            )
        super().__init__(origin, rotation, "spherical")

    def vector_components_to_local(
        self, global_vectors: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        local_vectors = super().vector_components_to_local(
            global_vectors, global_positions
        )
        local_positions = super().to_local(global_positions)
        x = local_positions[:, 0]
        y = local_positions[:, 1]
        z = local_positions[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        rho = np.sqrt(x * x + y * y)

        e_r = np.zeros_like(local_positions)
        e_theta = np.zeros_like(local_positions)
        e_phi = np.zeros_like(local_positions)
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
                inv_rho = 1.0 / rho[nonzero_rho]
                cos_theta = z[nonzero_rho] * inv_r[nonzero_rho]
                sin_theta = rho[nonzero_rho] * inv_r[nonzero_rho]
                cos_phi = x[nonzero_rho] * inv_rho
                sin_phi = y[nonzero_rho] * inv_rho
                e_theta[nonzero_rho, 0] = cos_theta * cos_phi
                e_theta[nonzero_rho, 1] = cos_theta * sin_phi
                e_theta[nonzero_rho, 2] = -sin_theta
                e_phi[nonzero_rho, 0] = -sin_phi
                e_phi[nonzero_rho, 1] = cos_phi
                e_phi[nonzero_rho, 2] = 0.0

        comp_r = np.einsum("ij,ij->i", e_r, local_vectors)
        comp_theta = np.einsum("ij,ij->i", e_theta, local_vectors)
        comp_phi = np.einsum("ij,ij->i", e_phi, local_vectors)
        return np.stack((comp_r, comp_theta, comp_phi), axis=1)

    def vector_components_to_global(
        self, local_components: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        local_positions = super().to_local(global_positions)
        x = local_positions[:, 0]
        y = local_positions[:, 1]
        z = local_positions[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        rho = np.sqrt(x * x + y * y)

        e_r = np.zeros_like(local_positions)
        e_theta = np.zeros_like(local_positions)
        e_phi = np.zeros_like(local_positions)
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
                inv_rho = 1.0 / rho[nonzero_rho]
                cos_theta = z[nonzero_rho] * inv_r[nonzero_rho]
                sin_theta = rho[nonzero_rho] * inv_r[nonzero_rho]
                cos_phi = x[nonzero_rho] * inv_rho
                sin_phi = y[nonzero_rho] * inv_rho
                e_theta[nonzero_rho, 0] = cos_theta * cos_phi
                e_theta[nonzero_rho, 1] = cos_theta * sin_phi
                e_theta[nonzero_rho, 2] = -sin_theta
                e_phi[nonzero_rho, 0] = -sin_phi
                e_phi[nonzero_rho, 1] = cos_phi
                e_phi[nonzero_rho, 2] = 0.0

        local_components = np.asarray(local_components, dtype=float)
        comp_r = local_components[:, 0]
        comp_theta = local_components[:, 1]
        comp_phi = local_components[:, 2]
        local_vectors = (
            e_r * comp_r[:, np.newaxis]
            + e_theta * comp_theta[:, np.newaxis]
            + e_phi * comp_phi[:, np.newaxis]
        )
        return local_vectors @ self._local_to_global_cart.T

    def to_local(self, global_points: np.ndarray) -> np.ndarray:
        local = super().to_local(global_points)
        x = local[:, 0]
        y = local[:, 1]
        z = local[:, 2]
        r = np.sqrt(x * x + y * y + z * z)
        theta = np.zeros_like(r)
        phi = np.zeros_like(r)
        nonzero = r > 0.0
        theta[nonzero] = np.arccos(z[nonzero] / r[nonzero])
        phi[nonzero] = np.arctan2(y[nonzero], x[nonzero])
        theta_deg = CoordinateSystem._rad_to_deg(theta)
        phi_deg = CoordinateSystem._rad_to_deg(phi)
        return np.stack((r, theta_deg, phi_deg), axis=1)

    def to_global(self, local_points: np.ndarray) -> np.ndarray:
        r = local_points[:, 0]
        theta = CoordinateSystem._deg_to_rad(local_points[:, 1])
        phi = CoordinateSystem._deg_to_rad(local_points[:, 2])
        sin_theta = np.sin(theta)
        x = r * sin_theta * np.cos(phi)
        y = r * sin_theta * np.sin(phi)
        z = r * np.cos(theta)
        local_cart = np.stack((x, y, z), axis=1)
        return self._origin + local_cart @ self._local_to_global_cart.T

from abc import ABC, abstractmethod
from typing import Any, Mapping, Optional, Sequence, Type, TypeVar, cast
from pytec3d.logging_utils import get_logger

import numpy as np


LOGGER = get_logger("pytec3d.coordinate_system")

X_AXIS_CYLINDRICAL_ANGLES = (0.0, 0.0, 90.0)
Y_AXIS_CYLINDRICAL_ANGLES = (0.0, -90.0, 0.0)
Z_AXIS_CYLINDRICAL_ANGLES = (0.0, 0.0, 0.0)

TCoordinateSystem = TypeVar("TCoordinateSystem", bound="CoordinateSystem")

class CoordinateSystem(ABC):
    @abstractmethod
    def __init__(
        self, origin: Sequence[float], local_to_global_cart: np.ndarray, cs_type: str
    ):
        """Store validated origin/rotation and coordinate-system type.

        Args:
            origin: Coordinate-system origin.
            local_to_global_cart: 3x3 local-to-global rotation matrix.
            cs_type: Coordinate-system type name.

        Returns:
            None
        """
        self._origin = CoordinateSystem._as_vector3(origin, "origin")
        self._local_to_global_cart = CoordinateSystem._as_rotation_matrix(
            local_to_global_cart,
            "local_to_global_cart",
        )
        if not isinstance(cs_type, str) or not cs_type.strip():
            raise ValueError("cs_type must be a non-empty string.")
        self._type = cs_type
        LOGGER.debug("Initialized %s coordinate system.", cs_type)

    @staticmethod
    def _as_vector3(values: Sequence[float], arg_name: str) -> np.ndarray:
        """Convert an input sequence to a finite shape-(3,) NumPy array.

        Args:
            values: Candidate vector.
            arg_name: Parameter name used in error messages.

        Returns:
            np.ndarray: Validated shape-(3,) vector.
        """
        arr = np.asarray(values, dtype=float)
        if arr.shape != (3,):
            raise ValueError(f"{arg_name} must have shape (3,), got {arr.shape}.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{arg_name} must contain only finite values.")
        return arr

    @staticmethod
    def _as_points3(values: np.ndarray, arg_name: str) -> np.ndarray:
        """Convert points/vectors to a finite shape-(N, 3) NumPy array.

        Args:
            values: Candidate points/vectors array.
            arg_name: Parameter name used in error messages.

        Returns:
            np.ndarray: Validated shape-(N, 3) array.
        """
        arr = np.asarray(values, dtype=float)
        if arr.ndim == 1:
            if arr.shape[0] != 3:
                raise ValueError(
                    f"{arg_name} must have shape (N, 3) or (3,), got {arr.shape}."
                )
            arr = arr.reshape(1, 3)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(
                f"{arg_name} must have shape (N, 3) or (3,), got {arr.shape}."
            )
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{arg_name} must contain only finite values.")
        return arr

    @staticmethod
    def _as_rotation_matrix(values: np.ndarray, arg_name: str) -> np.ndarray:
        """Validate a finite orthonormal 3x3 rotation matrix.

        Args:
            values: Candidate rotation matrix.
            arg_name: Parameter name used in error messages.

        Returns:
            np.ndarray: Validated shape-(3, 3) rotation matrix.
        """
        arr = np.asarray(values, dtype=float)
        if arr.shape != (3, 3):
            raise ValueError(f"{arg_name} must have shape (3, 3), got {arr.shape}.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{arg_name} must contain only finite values.")
        gram = arr.T @ arr
        if not np.allclose(gram, np.eye(3), atol=1e-8):
            raise ValueError(
                f"{arg_name} must be orthonormal (R^T R == I within tolerance)."
            )
        return arr

    @staticmethod
    def _validate_vector_mapping_inputs(
        vectors: np.ndarray,
        positions: np.ndarray,
        vectors_name: str,
        positions_name: str,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Validate vector/position arrays and ensure matching row counts.

        Args:
            vectors: Vector-components array.
            positions: Associated position array.
            vectors_name: Vector argument name.
            positions_name: Position argument name.

        Returns:
            tuple[np.ndarray, np.ndarray]: Validated vectors and positions.
        """
        vectors_arr = CoordinateSystem._as_points3(vectors, vectors_name)
        positions_arr = CoordinateSystem._as_points3(positions, positions_name)
        if vectors_arr.shape[0] != positions_arr.shape[0]:
            raise ValueError(
                f"{vectors_name} and {positions_name} must have the same number of rows."
            )
        return vectors_arr, positions_arr

    @property
    def type(self) -> str:
        """Return coordinate-system type identifier.

        Returns:
            str: Coordinate-system type.
        """
        return self._type

    @property
    def origin(self) -> np.ndarray:
        """Return coordinate-system origin.

        Returns:
            np.ndarray: Origin vector.
        """
        return self._origin

    @property
    def local_to_global_cart(self) -> np.ndarray:
        """Return the 3x3 local-to-global Cartesian rotation matrix.

        Returns:
            np.ndarray: Rotation matrix (a copy).
        """
        return self._local_to_global_cart.copy()

    def axis_direction(self) -> np.ndarray:
        """Return the global direction of the local z-axis.

        Returns:
            np.ndarray: Direction vector.
        """
        return self._local_to_global_cart[:, 2].copy()

    def local_to_global_matrices(self, global_positions: np.ndarray) -> np.ndarray:
        """Build local-to-global rotation matrices evaluated at positions.

        Args:
            global_positions: Global Cartesian positions with shape (N, 3) or (3,).

        Returns:
            np.ndarray: Rotation matrices with shape (N, 3, 3).
        """
        positions = CoordinateSystem._as_points3(global_positions, "global_positions")
        n = positions.shape[0]
        mats = np.zeros((n, 3, 3), dtype=float)
        basis = np.eye(3, dtype=float)
        for axis in range(3):
            local_axis = np.repeat(basis[axis : axis + 1, :], n, axis=0)
            global_axis = self.vector_components_to_global(local_axis, positions)
            mats[:, :, axis] = global_axis
        return mats

    def local_to_global_matrix(self, global_position: np.ndarray) -> np.ndarray:
        """Build one local-to-global rotation matrix at a single position.

        Args:
            global_position: Global Cartesian position with shape (3,).

        Returns:
            np.ndarray: Rotation matrix with shape (3, 3).
        """
        return self.local_to_global_matrices(global_position)[0]

    def rotate_vector_to_global(
        self, local_components: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        """Rotate local vector components into global Cartesian components.

        Args:
            local_components: Local components with shape (N, 3) or (3,).
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Rotated global Cartesian vectors.
        """
        return self.vector_components_to_global(local_components, global_positions)

    def rotate_tensor_to_global(
        self, local_tensors: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        """Rotate local tensors into global Cartesian tensors.

        Args:
            local_tensors: Tensor(s) with shape (3, 3) or (N, 3, 3).
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Rotated tensor(s), shape (3, 3) or (N, 3, 3).
        """
        positions = CoordinateSystem._as_points3(global_positions, "global_positions")
        mats = self.local_to_global_matrices(positions)

        tensors = np.asarray(local_tensors, dtype=float)
        single = False
        if tensors.shape == (3, 3):
            tensors = tensors[np.newaxis, :, :]
            single = True
        if tensors.ndim != 3 or tensors.shape[1:] != (3, 3):
            raise ValueError("local_tensors must have shape (3, 3) or (N, 3, 3).")
        if tensors.shape[0] != positions.shape[0]:
            if tensors.shape[0] == 1:
                tensors = np.repeat(tensors, positions.shape[0], axis=0)
            else:
                raise ValueError("local_tensors and global_positions must have matching row counts.")

        rotated = np.einsum("nij,njk,nlk->nil", mats, tensors, mats)
        return rotated[0] if single else rotated

    def to_local(self, global_points: np.ndarray) -> np.ndarray:
        """Map global Cartesian points into local Cartesian coordinates.

        Args:
            global_points: Global Cartesian points.

        Returns:
            np.ndarray: Local Cartesian points.
        """
        points = CoordinateSystem._as_points3(global_points, "global_points")
        return (points - self._origin) @ self._local_to_global_cart

    def to_global(self, local_points: np.ndarray) -> np.ndarray:
        """Map local Cartesian points into global Cartesian coordinates.

        Args:
            local_points: Local Cartesian points.

        Returns:
            np.ndarray: Global Cartesian points.
        """
        points = CoordinateSystem._as_points3(local_points, "local_points")
        return self._origin + points @ self._local_to_global_cart.T

    def vector_components_to_local(
        self, global_vectors: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        """Map vector components from global Cartesian basis to local Cartesian basis.

        Args:
            global_vectors: Global Cartesian vector components.
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Local Cartesian vector components.
        """
        vectors, _ = CoordinateSystem._validate_vector_mapping_inputs(
            global_vectors,
            global_positions,
            "global_vectors",
            "global_positions",
        )
        return vectors @ self._local_to_global_cart

    def vector_components_to_global(
        self, local_components: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        """Map vector components from local Cartesian basis to global Cartesian basis.

        Args:
            local_components: Local Cartesian vector components.
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Global Cartesian vector components.
        """
        components, _ = CoordinateSystem._validate_vector_mapping_inputs(
            local_components,
            global_positions,
            "local_components",
            "global_positions",
        )
        return components @ self._local_to_global_cart.T

    @staticmethod
    def _deg_to_rad(values: np.ndarray) -> np.ndarray:
        """Convert degrees to radians.

        Args:
            values: Values in degrees.

        Returns:
            np.ndarray: Values in radians.
        """
        return np.deg2rad(values)

    @staticmethod
    def _rad_to_deg(values: np.ndarray) -> np.ndarray:
        """Convert radians to degrees.

        Args:
            values: Values in radians.

        Returns:
            np.ndarray: Values in degrees.
        """
        return np.rad2deg(values)

    @staticmethod
    def _normalize_or_error(vec: np.ndarray, error_message: str) -> np.ndarray:
        """Normalize a vector or raise ValueError when its norm is zero.

        Args:
            vec: Candidate vector.
            error_message: Error message for invalid vectors.

        Returns:
            np.ndarray: Unit-normalized vector.
        """
        norm = np.linalg.norm(vec)
        if norm <= 0.0:
            raise ValueError(error_message)
        return vec / norm

    @staticmethod
    def _build_local_to_global_from_axes(
        x_axis: np.ndarray, y_axis: np.ndarray, z_axis: np.ndarray
    ) -> np.ndarray:
        """Build a rotation matrix from local basis vectors expressed in global coordinates.

        Args:
            x_axis: Local x-axis in global coordinates.
            y_axis: Local y-axis in global coordinates.
            z_axis: Local z-axis in global coordinates.

        Returns:
            np.ndarray: Rotation matrix.
        """
        rotation = np.zeros((3, 3), dtype=float)
        rotation[:, 0] = x_axis
        rotation[:, 1] = y_axis
        rotation[:, 2] = z_axis
        return rotation

    @staticmethod
    def _build_local_to_global_from_points(
        origin: np.ndarray, nxax: np.ndarray, nxypl: np.ndarray
    ) -> np.ndarray:
        """Build a right-handed rotation matrix from origin, x-axis point, and xy-plane point.

        Args:
            origin: Origin point.
            nxax: Point defining local x-axis direction.
            nxypl: Point defining local xy-plane orientation.

        Returns:
            np.ndarray: Rotation matrix.
        """
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
        """Build one local-to-global rotation matrix from Euler-like angle triplet in degrees.

        Args:
            angle_xy_deg: Angle in XY plane.
            angle_yz_deg: Angle in YZ plane.
            angle_xz_deg: Angle in XZ plane.

        Returns:
            np.ndarray: Rotation matrix.
        """
        for name, value in (
            ("angle_xy_deg", angle_xy_deg),
            ("angle_yz_deg", angle_yz_deg),
            ("angle_xz_deg", angle_xz_deg),
        ):
            if not np.isfinite(float(value)):
                raise ValueError(f"{name} must be finite.")
        angles = np.array([[angle_xy_deg, angle_yz_deg, angle_xz_deg]], dtype=float)
        coeffs = CoordinateSystem.build_local_to_global_coefficients(angles)
        return coeffs[0].reshape(3, 3)

    @staticmethod
    def build_local_to_global_coefficients(angles_deg: np.ndarray) -> np.ndarray:
        """Build flattened local-to-global rotation coefficients for an array of angle triplets.

        Args:
            angles_deg: Array with shape (N, 3).

        Returns:
            np.ndarray: Flattened rotation coefficients with shape (N, 9).
        """
        angles = np.asarray(angles_deg, dtype=float)
        if angles.ndim != 2 or angles.shape[1] != 3:
            raise ValueError("angles_deg must have shape (N, 3).")
        if not np.all(np.isfinite(angles)):
            raise ValueError("angles_deg must contain only finite values.")

        eps = 1e-12

        def snap_trig_values(values: np.ndarray) -> np.ndarray:
            out = values.copy()
            out[np.abs(out) < eps] = 0.0
            out[np.abs(out - 1.0) < eps] = 1.0
            out[np.abs(out + 1.0) < eps] = -1.0
            return out

        angles_rad = CoordinateSystem._deg_to_rad(angles)
        angle_xy = angles_rad[:, 0]
        angle_yz = angles_rad[:, 1]
        angle_xz = angles_rad[:, 2]

        c1 = snap_trig_values(np.cos(angle_xy))
        s1 = snap_trig_values(np.sin(angle_xy))
        c2 = snap_trig_values(np.cos(angle_yz))
        s2 = snap_trig_values(np.sin(angle_yz))
        c3 = snap_trig_values(np.cos(angle_xz))
        s3 = snap_trig_values(np.sin(angle_xz))

        n_rows = angles.shape[0]
        coeffs = np.empty((n_rows, 9), dtype=float)
        coeffs[:, 0] = c1 * c3 - s1 * s2 * s3
        coeffs[:, 1] = -s1 * c2
        coeffs[:, 2] = c1 * s3 + s1 * s2 * c3
        coeffs[:, 3] = s1 * c3 + c1 * s2 * s3
        coeffs[:, 4] = c1 * c2
        coeffs[:, 5] = s1 * s3 - c1 * s2 * c3
        coeffs[:, 6] = -c2 * s3
        coeffs[:, 7] = s2
        coeffs[:, 8] = c2 * c3
        return coeffs

    @staticmethod
    def build_global_to_local_matrix(
        angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float
    ) -> np.ndarray:
        """Build one global-to-local rotation matrix from angle triplet in degrees.

        Args:
            angle_xy_deg: Angle in XY plane.
            angle_yz_deg: Angle in YZ plane.
            angle_xz_deg: Angle in XZ plane.

        Returns:
            np.ndarray: Rotation matrix.
        """
        return CoordinateSystem.build_local_to_global_matrix(
            angle_xy_deg, angle_yz_deg, angle_xz_deg
        ).T

    @staticmethod
    def build_global_to_local(
        angle_xy_deg: float, angle_yz_deg: float, angle_xz_deg: float
    ) -> np.ndarray:
        """Backward-compatible alias for build_global_to_local_matrix.

        Args:
            angle_xy_deg: Angle in XY plane.
            angle_yz_deg: Angle in YZ plane.
            angle_xz_deg: Angle in XZ plane.

        Returns:
            np.ndarray: Rotation matrix.
        """
        return CoordinateSystem.build_global_to_local_matrix(
            angle_xy_deg, angle_yz_deg, angle_xz_deg
        )

    @classmethod
    def from_angles(
        cls: Type[TCoordinateSystem],
        origin: Sequence[float],
        angles_deg: Sequence[float],
    ) -> TCoordinateSystem:
        """Construct a coordinate system from origin and three orientation angles.

        Args:
            origin: Origin point.
            angles_deg: Sequence of three orientation angles in degrees.

        Returns:
            TCoordinateSystem: Constructed coordinate-system instance.
        """
        if not isinstance(angles_deg, Sequence) or len(angles_deg) != 3:
            raise ValueError("angles_deg must be a sequence with exactly 3 values.")
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
        """Construct a coordinate system from origin and two reference points.

        Args:
            origin: Origin point.
            nxax: Point defining local x-axis direction.
            nxypl: Point defining local xy-plane orientation.

        Returns:
            TCoordinateSystem: Constructed coordinate-system instance.
        """
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
        local_to_global_cart: Optional[np.ndarray] = None,
    ):
        """Create a Cartesian coordinate system from either angles or defining points.

        Args:
            origin: Origin point.
            nxax: Optional point defining local x-axis direction.
            nxypl: Optional point defining local xy-plane orientation.
            thxy_deg: Optional XY angle.
            thyz_deg: Optional YZ angle.
            thzx_deg: Optional ZX angle.

        Returns:
            None
        """
        angles = (thxy_deg, thyz_deg, thzx_deg)
        if any(angle is not None for angle in angles) and not all(
            angle is not None for angle in angles
        ):
            raise ValueError(
                "thxy_deg, thyz_deg and thzx_deg must all be provided together."
            )

        if local_to_global_cart is not None:
            if nxax is not None or nxypl is not None or any(angle is not None for angle in angles):
                raise ValueError(
                    "local_to_global_cart cannot be combined with nxax/nxypl/angles."
                )
            rotation = CoordinateSystem._as_rotation_matrix(
                local_to_global_cart, "local_to_global_cart"
            )
        elif thxy_deg is not None and thyz_deg is not None and thzx_deg is not None:
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
        """Convert global vector components to local Cartesian components.

        Args:
            global_vectors: Global Cartesian vectors.
            global_positions: Associated positions.

        Returns:
            np.ndarray: Local Cartesian vectors.
        """
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
        local_to_global_cart: Optional[np.ndarray] = None,
    ):
        """Create a cylindrical coordinate system from either angles or defining points.

        Args:
            origin: Origin point.
            nxax: Optional point defining local x-axis direction.
            nxypl: Optional point defining local xy-plane orientation.
            thxy_deg: Optional XY angle.
            thyz_deg: Optional YZ angle.
            thzx_deg: Optional ZX angle.

        Returns:
            None
        """
        angles = (thxy_deg, thyz_deg, thzx_deg)
        if any(angle is not None for angle in angles) and not all(
            angle is not None for angle in angles
        ):
            raise ValueError(
                "thxy_deg, thyz_deg and thzx_deg must all be provided together."
            )

        if local_to_global_cart is not None:
            if nxax is not None or nxypl is not None or any(angle is not None for angle in angles):
                raise ValueError(
                    "local_to_global_cart cannot be combined with nxax/nxypl/angles."
                )
            rotation = CoordinateSystem._as_rotation_matrix(
                local_to_global_cart, "local_to_global_cart"
            )
        elif thxy_deg is not None and thyz_deg is not None and thzx_deg is not None:
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
        """Convert global Cartesian vector components into (r, theta, z) components.

        Args:
            global_vectors: Global Cartesian vectors.
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Cylindrical vector components.
        """
        local_vectors = super().vector_components_to_local(
            global_vectors, global_positions
        )
        local_positions = super().to_local(global_positions)
        x = local_positions[:, 0]
        y = local_positions[:, 1]
        r = np.sqrt(x * x + y * y)

        local_components = np.empty_like(local_vectors)
        local_components[:, 0] = local_vectors[:, 0]
        local_components[:, 1] = local_vectors[:, 1]
        local_components[:, 2] = local_vectors[:, 2]
        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            vx = local_vectors[nonzero, 0]
            vy = local_vectors[nonzero, 1]
            local_components[nonzero, 0] = (x[nonzero] * vx + y[nonzero] * vy) * inv_r
            local_components[nonzero, 1] = (-y[nonzero] * vx + x[nonzero] * vy) * inv_r
        return local_components

    def vector_components_to_global(
        self, local_components: np.ndarray, global_positions: np.ndarray
    ) -> np.ndarray:
        """Convert (r, theta, z) components into global Cartesian components.

        Args:
            local_components: Cylindrical vector components.
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Global Cartesian vectors.
        """
        local_positions = super().to_local(global_positions)
        x = local_positions[:, 0]
        y = local_positions[:, 1]
        r = np.sqrt(x * x + y * y)

        local_components = np.asarray(local_components, dtype=float)
        comp_r = local_components[:, 0]
        comp_theta = local_components[:, 1]
        comp_z = local_components[:, 2]
        local_vectors = np.empty_like(local_components)
        local_vectors[:, 0] = comp_r
        local_vectors[:, 1] = comp_theta
        local_vectors[:, 2] = comp_z
        nonzero = r > 0.0
        if np.any(nonzero):
            inv_r = 1.0 / r[nonzero]
            xr = x[nonzero] * inv_r
            yr = y[nonzero] * inv_r
            local_vectors[nonzero, 0] = xr * comp_r[nonzero] - yr * comp_theta[nonzero]
            local_vectors[nonzero, 1] = yr * comp_r[nonzero] + xr * comp_theta[nonzero]
        return local_vectors @ self._local_to_global_cart.T

    def to_local(self, global_points: np.ndarray) -> np.ndarray:
        """Map global points to (r, theta_deg, z) cylindrical coordinates.

        Args:
            global_points: Global Cartesian points.

        Returns:
            np.ndarray: Cylindrical coordinates.
        """
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
        """Map cylindrical coordinates (r, theta_deg, z) back to global Cartesian points.

        Args:
            local_points: Cylindrical coordinates.

        Returns:
            np.ndarray: Global Cartesian points.
        """
        local_points = CoordinateSystem._as_points3(local_points, "local_points")
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
        local_to_global_cart: Optional[np.ndarray] = None,
    ):
        """Create a spherical coordinate system from either angles or defining points.

        Args:
            origin: Origin point.
            nxax: Optional point defining local x-axis direction.
            nxypl: Optional point defining local xy-plane orientation.
            thxy_deg: Optional XY angle.
            thyz_deg: Optional YZ angle.
            thzx_deg: Optional ZX angle.

        Returns:
            None
        """
        angles = (thxy_deg, thyz_deg, thzx_deg)
        if any(angle is not None for angle in angles) and not all(
            angle is not None for angle in angles
        ):
            raise ValueError(
                "thxy_deg, thyz_deg and thzx_deg must all be provided together."
            )

        if local_to_global_cart is not None:
            if nxax is not None or nxypl is not None or any(angle is not None for angle in angles):
                raise ValueError(
                    "local_to_global_cart cannot be combined with nxax/nxypl/angles."
                )
            rotation = CoordinateSystem._as_rotation_matrix(
                local_to_global_cart, "local_to_global_cart"
            )
        elif thxy_deg is not None and thyz_deg is not None and thzx_deg is not None:
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
        """Convert global Cartesian vector components into (r, theta, phi) components.

        Args:
            global_vectors: Global Cartesian vectors.
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Spherical vector components.
        """
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
        """Convert (r, theta, phi) components into global Cartesian components.

        Args:
            local_components: Spherical vector components.
            global_positions: Associated global positions.

        Returns:
            np.ndarray: Global Cartesian vectors.
        """
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
        """Map global points to (r, theta_deg, phi_deg) spherical coordinates.

        Args:
            global_points: Global Cartesian points.

        Returns:
            np.ndarray: Spherical coordinates.
        """
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
        """Map spherical coordinates (r, theta_deg, phi_deg) back to global Cartesian points.

        Args:
            local_points: Spherical coordinates.

        Returns:
            np.ndarray: Global Cartesian points.
        """
        local_points = CoordinateSystem._as_points3(local_points, "local_points")
        r = local_points[:, 0]
        theta = CoordinateSystem._deg_to_rad(local_points[:, 1])
        phi = CoordinateSystem._deg_to_rad(local_points[:, 2])
        sin_theta = np.sin(theta)
        x = r * sin_theta * np.cos(phi)
        y = r * sin_theta * np.sin(phi)
        z = r * np.cos(theta)
        local_cart = np.stack((x, y, z), axis=1)
        return self._origin + local_cart @ self._local_to_global_cart.T


class XAxisCylindricalCoordinateSystem(CylindricalCoordinateSystem):
    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        local_to_global_cart: Optional[np.ndarray] = None,
    ):
        if local_to_global_cart is not None:
            super().__init__(origin=origin, local_to_global_cart=local_to_global_cart)
            return
        super().__init__(
            origin=origin,
            thxy_deg=X_AXIS_CYLINDRICAL_ANGLES[0],
            thyz_deg=X_AXIS_CYLINDRICAL_ANGLES[1],
            thzx_deg=X_AXIS_CYLINDRICAL_ANGLES[2],
        )


class YAxisCylindricalCoordinateSystem(CylindricalCoordinateSystem):
    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        local_to_global_cart: Optional[np.ndarray] = None,
    ):
        if local_to_global_cart is not None:
            super().__init__(origin=origin, local_to_global_cart=local_to_global_cart)
            return
        super().__init__(
            origin=origin,
            thxy_deg=Y_AXIS_CYLINDRICAL_ANGLES[0],
            thyz_deg=Y_AXIS_CYLINDRICAL_ANGLES[1],
            thzx_deg=Y_AXIS_CYLINDRICAL_ANGLES[2],
        )


class ZAxisCylindricalCoordinateSystem(CylindricalCoordinateSystem):
    def __init__(
        self,
        origin: Sequence[float] = (0.0, 0.0, 0.0),
        local_to_global_cart: Optional[np.ndarray] = None,
    ):
        if local_to_global_cart is not None:
            super().__init__(origin=origin, local_to_global_cart=local_to_global_cart)
            return
        super().__init__(
            origin=origin,
            thxy_deg=Z_AXIS_CYLINDRICAL_ANGLES[0],
            thyz_deg=Z_AXIS_CYLINDRICAL_ANGLES[1],
            thzx_deg=Z_AXIS_CYLINDRICAL_ANGLES[2],
        )


#: Names of the "type" attribute of a coordinate system in the model data JSON.
COORDINATE_SYSTEM_TYPES = {
    "cartesian": CartesianCoordinateSystem,
    "cylindrical": CylindricalCoordinateSystem,
    "x_axis_cylindrical": XAxisCylindricalCoordinateSystem,
    "y_axis_cylindrical": YAxisCylindricalCoordinateSystem,
    "z_axis_cylindrical": ZAxisCylindricalCoordinateSystem,
    "spherical": SphericalCoordinateSystem,
}

_AXIS_CYLINDRICAL_TYPES = (
    XAxisCylindricalCoordinateSystem,
    YAxisCylindricalCoordinateSystem,
    ZAxisCylindricalCoordinateSystem,
)


def coordinate_system_from_dict(spec: Mapping[str, Any]) -> CoordinateSystem:
    """Build a coordinate system from its attributes in the model data JSON.

    Attributes:
        type: One of ``COORDINATE_SYSTEM_TYPES`` (required).
        origin: Origin point, default (0, 0, 0).
        matrix: 3x3 local-to-global rotation matrix, or
        angles: [thxy, thyz, thzx] in degrees, or
        points: {"nxax": [...], "nxypl": [...]} defining the local x axis and xy plane.
            At most one of the three; none gives the default orientation of the type.
            The x/y/z axis cylindrical types only take ``matrix``.

    Args:
        spec: Mapping with the attributes above.

    Returns:
        CoordinateSystem: Instance of the class named by ``type``.

    Raises:
        ValueError: On an unknown type or attribute, or an invalid combination.
    """
    if not isinstance(spec, Mapping):
        raise ValueError(f"A coordinate system must be an object, got {type(spec).__name__}.")
    unknown = set(spec) - {"type", "origin", "matrix", "angles", "points"}
    if unknown:
        raise ValueError(f"Unknown coordinate system attributes: {sorted(unknown)}.")

    type_name = spec.get("type")
    cls = COORDINATE_SYSTEM_TYPES.get(type_name)
    if cls is None:
        raise ValueError(
            f"Coordinate system type must be one of {sorted(COORDINATE_SYSTEM_TYPES)}, got {type_name!r}."
        )

    orientations = [key for key in ("matrix", "angles", "points") if spec.get(key) is not None]
    if len(orientations) > 1:
        raise ValueError(f"Give one of matrix, angles or points, not {orientations}.")
    if cls in _AXIS_CYLINDRICAL_TYPES and orientations not in ([], ["matrix"]):
        raise ValueError(f"The type {type_name!r} has a fixed orientation; only 'matrix' is accepted.")

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
    """Return the model data JSON attributes (type, origin, matrix) of a coordinate system."""
    return {
        "type": coordinate_system.type,
        "origin": coordinate_system.origin.tolist(),
        "matrix": coordinate_system.local_to_global_cart.tolist(),
    }

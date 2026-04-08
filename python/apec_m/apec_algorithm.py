from __future__ import annotations

import copy
import os
from typing import Any, Callable, Dict, Optional, TYPE_CHECKING, Tuple

import numpy as np
import pandas as pd
import vtk
from .selector import Selector
from .model_edit import ModelEdit
from .constants import (
    DEFAULT_APEC_ALGORITHM_INPUT_PARAMETERS,
    DEFAULT_APEC_ALGORITHM_STATUS,
)

import vtk.util.numpy_support as vtknp

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    from .apec_m import APEC_M

class ApecAlgorithm:
    """
    Python scaffold of the legacy APEC equilibrium/evolution analysis objects.

    Notes:
        - This class mirrors the old C++ analysis API shape.
        - Methods that previously depended on libMesh internals are intentionally
          left as VTK-grid placeholders and can be refined incrementally.
        - Only time-step logic is supported (no substeps).
    """

    STEP_FIND_SEPARATRIX = "find_separatrix"
    STEP_COMPUTE_JPLASMA = "compute_jplasma_gradshafranov"
    STEP_CONTROL_CURRENTS = "control_currents"
    STEP_CHECK_CONVERGENCE = "check_convergence"

    def __init__(
        self,
        apec: "APEC_M",
        input_parameters: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.apec: "APEC_M" = apec
        self.input_parameters: Dict[str, Any] = copy.deepcopy(
            DEFAULT_APEC_ALGORITHM_INPUT_PARAMETERS
        )
        if input_parameters:
            self.input_parameters.update(input_parameters)

        self.status: Dict[str, Any] = copy.deepcopy(DEFAULT_APEC_ALGORITHM_STATUS)

        self._first_absolute_iteration: bool = True
        self._previous_flux_vector: Optional[np.ndarray] = None
        self._old_control_current_density_base: Optional[np.ndarray] = None
        self._diagnostics_frames: Dict[str, pd.DataFrame] = {}
        self._step_overrides: Dict[str, Callable[..., Any]] = {}
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)
        self._supported_step_overrides = {
            self.STEP_FIND_SEPARATRIX,
            self.STEP_COMPUTE_JPLASMA,
            self.STEP_CONTROL_CURRENTS,
            self.STEP_CHECK_CONVERGENCE,
        }

    def reset(self) -> None:
        self.status = copy.deepcopy(DEFAULT_APEC_ALGORITHM_STATUS)
        self._first_absolute_iteration = True
        self._previous_flux_vector = None
        self._old_control_current_density_base = None
        self._diagnostics_frames.clear()
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)

    @staticmethod
    def _vector_component_indices(mesh_dimension: int) -> Tuple[int, int, int]:
        if mesh_dimension == 2:
            # Legacy C++ labels for 2D: ["_z", "_x", "_y"]
            return 2, 0, 1
        # Legacy C++ labels for 3D: ["_y", "_x", "_z"]
        return 1, 0, 2

    def set_toroidal_slices(self) -> None:
        nodes = np.asarray(
            self.input_parameters["toroidal_slices_nodes"], dtype=np.int64
        )
        elems = np.asarray(
            self.input_parameters["toroidal_slices_elems"], dtype=np.int64
        )
        if nodes.ndim != 2 or nodes.shape[1] < 2:
            raise ValueError(
                "input_parameters['toroidal_slices_nodes'] must be a 2D array "
                "with at least 2 columns (source ids + >=1 slice)."
            )
        if elems.ndim != 2 or elems.shape[1] < 2:
            raise ValueError(
                "input_parameters['toroidal_slices_elems'] must be a 2D array "
                "with at least 2 columns (source ids + >=1 slice)."
            )
        self._toroidal_slices_nodes = nodes
        self._toroidal_slices_elems = elems

    def _require_toroidal_slices(self) -> Tuple[np.ndarray, np.ndarray]:
        if (
            self._toroidal_slices_nodes.size == 0
            or self._toroidal_slices_elems.size == 0
        ):
            self.set_toroidal_slices()
        return self._toroidal_slices_nodes, self._toroidal_slices_elems

    @staticmethod
    def _average_from_slice_targets(
        source_values: np.ndarray,
        source_ids: np.ndarray,
        target_ids_matrix: np.ndarray,
        source_size: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        valid_source = (source_ids >= 0) & (source_ids < source_size)
        if not np.any(valid_source):
            return (
                np.empty((0,), dtype=np.int64),
                np.empty((0,), dtype=float),
            )

        src = source_ids[valid_source]
        dst = target_ids_matrix[valid_source]
        valid_dst = (dst >= 0) & (dst < source_values.shape[0])

        clipped_dst = np.where(valid_dst, dst, 0)
        gathered = source_values[clipped_dst]
        gathered = np.where(valid_dst, gathered, 0.0)

        counts = np.sum(valid_dst, axis=1)
        rows = counts > 0
        if not np.any(rows):
            return (
                np.empty((0,), dtype=np.int64),
                np.empty((0,), dtype=float),
            )

        means = np.sum(gathered, axis=1)[rows] / counts[rows]
        return src[rows], np.asarray(means, dtype=float)

    @staticmethod
    def _broadcast_rows_to_slice_targets(
        source_values: np.ndarray,
        source_ids: np.ndarray,
        target_ids_matrix: np.ndarray,
        target_size: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        valid_source = (source_ids >= 0) & (source_ids < source_values.shape[0])
        if not np.any(valid_source):
            return (
                np.empty((0,), dtype=np.int64),
                np.empty((0,), dtype=float),
            )

        src = source_ids[valid_source]
        dst = target_ids_matrix[valid_source]

        repeated_src_values = np.repeat(source_values[src], dst.shape[1])
        dst_flat = dst.reshape(-1)

        valid_dst = (dst_flat >= 0) & (dst_flat < target_size)
        return (
            dst_flat[valid_dst].astype(np.int64, copy=False),
            repeated_src_values[valid_dst].astype(float, copy=False),
        )

    def _get_mesh(self, mesh_name: Optional[str] = None) -> "vtk.vtkUnstructuredGrid":
        if mesh_name is None:
            return self.apec.mesh
        mesh = self.apec.get_mesh(mesh_name)
        if mesh is None:
            raise ValueError(f"Mesh '{mesh_name}' does not exist.")
        return mesh

    def _diagnostics_mesh_key(self, mesh_name: Optional[str] = None) -> str:
        if mesh_name is not None:
            return str(mesh_name)
        if self.apec.active_mesh_name is not None:
            return str(self.apec.active_mesh_name)
        return "__active__"

    def _diagnostics_csv_path(self, mesh_name: Optional[str] = None) -> str:
        _ = mesh_name
        diagnostics_dir = os.path.join(os.getcwd(), "__apec__")
        os.makedirs(diagnostics_dir, exist_ok=True)
        return os.path.join(diagnostics_dir, "apec_iteration_diagnostics.csv")

    def _build_diagnostics_row(self, converged: bool) -> Dict[str, Any]:
        row: Dict[str, Any] = {
            "time": float(self.apec.time_step["time"]),
            "picard_iteration": int(self.status["picard_iteration"]),
            "flux_metric": float(self.status["flux_norm"]),
            "flux_threshold": float(self.status["flux_tolerance"]),
            "converged": bool(converged),
        }

        for key, value in self.status.items():
            if key in (
                "picard_iteration",
                "converged",
                "flux_norm",
                "flux_tolerance",
            ):
                continue
            if key == "magnetic_axis":
                row["magnetic_axis_x"] = float(value[0])
                row["magnetic_axis_y"] = float(value[1])
                row["magnetic_axis_z"] = float(value[2])
            elif key == "limiter_point":
                row["limiter_x"] = float(value[0])
                row["limiter_y"] = float(value[1])
                row["limiter_z"] = float(value[2])
            elif key == "x_point":
                row["x_point_x"] = float(value[0])
                row["x_point_y"] = float(value[1])
                row["x_point_z"] = float(value[2])
            elif key == "b_target_axis":
                row["b_target_axis_x"] = float(value[0])
                row["b_target_axis_y"] = float(value[1])
            else:
                row[key] = value
        return row

    def set_step_override(
        self,
        step: str,
        override: Optional[Callable[..., Any]],
    ) -> None:
        if step not in self._supported_step_overrides:
            raise ValueError(
                f"Unsupported step '{step}'. "
                f"Supported steps: {sorted(self._supported_step_overrides)}."
            )
        if override is None:
            self._step_overrides.pop(step, None)
            return
        if not callable(override):
            raise TypeError("override must be callable or None.")
        self._step_overrides[step] = override

    def clear_step_overrides(self) -> None:
        self._step_overrides.clear()

    def _dispatch_step(
        self,
        step: str,
        default_impl: Callable[..., Any],
        mesh_name: Optional[str] = None,
    ) -> Any:
        override = self._step_overrides[step] if step in self._step_overrides else None
        if override is None:
            return default_impl(mesh_name=mesh_name)
        return override(self, mesh_name=mesh_name)

    def set_psi_magnetic_axis(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)

        psi_array = mesh.GetPointData().GetArray("PSI")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsle(mesh)
        selection_array = mesh.GetPointData().GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        psi = np.asarray(vtknp.vtk_to_numpy(psi_array), dtype=float).reshape(-1)
        psi_selected = psi[selected_ids]
        axis_local_index = int(np.argmax(psi_selected))
        axis_index = int(selected_ids[axis_local_index])
        psi_axis = float(psi[axis_index])

        points = mesh.GetPoints()
        xyz = vtknp.vtk_to_numpy(points.GetData())
        axis_point = tuple(float(v) for v in xyz[axis_index])

        self.status["psi_magnetic_axis"] = psi_axis
        self.status["psi_max"] = psi_axis
        self.status["magnetic_axis"] = axis_point

    def reset_plasma_components(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)
        cell_data = mesh.GetCellData()

        js_array = cell_data.GetArray("Js")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.set_value_on_selection(mesh, "PLASMA", 1, type="elem")

        selection_array = cell_data.GetArray("selection_cell")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        js_np = vtknp.vtk_to_numpy(js_array)
        if js_np.ndim == 1:
            js_np[selected_ids] = 0.0
        else:
            js_np[selected_ids, :] = 0.0
        js_array.Modified()

    def identify_plasma_region(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)

        psi_min = float(self.status["psi_min"])
        if not np.isfinite(psi_min):
            raise ValueError("status['psi_min'] is not set; cannot identify plasma region.")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")

        # Strict PSI < psi_min using an upper bound just below psi_min.
        upper_strict = np.nextafter(psi_min, -np.inf)
        Selector.select(
            mesh,
            "PSI",
            (-np.inf, upper_strict),
            type="node",
            mode="r",
        )

        Selector.esln(mesh, all_nodes=False, mode="r")
        Selector.set_value_on_selection(mesh, "PLASMA", 0, type="elem")

    def remove_plasma_out_of_X_point(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)
        cell_data = mesh.GetCellData()

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsel_ext(mesh)
        selection_array = cell_data.GetArray("selection_cell")
        vacuum_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
        vacuum_ids = np.flatnonzero(vacuum_mask)

        centroids_array = cell_data.GetArray("cell_centroids")
        centroids = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
        x_centroid = np.asarray(centroids[:, 0], dtype=float)
        y_centroid = np.asarray(centroids[:, 1], dtype=float)

        magnetic_axis = self.status["magnetic_axis"]
        axis_x = float(magnetic_axis[0])
        axis_y = float(magnetic_axis[1])

        plasma_current = float(self.input_parameters["plasma_current"])
        current_sign = 1.0 if plasma_current >= 0.0 else -1.0

        b_cell_array = cell_data.GetArray("Be")
        b_cell = np.asarray(vtknp.vtk_to_numpy(b_cell_array), dtype=float)
        bx = b_cell[:, 0]
        by = b_cell[:, 1]

        b_rot = current_sign * (
            by[vacuum_ids] * (x_centroid[vacuum_ids] - axis_x)
            - bx[vacuum_ids] * (y_centroid[vacuum_ids] - axis_y)
        )
        excluded_local_mask = b_rot >= 0.0

        excluded_cell_ids = vacuum_ids[excluded_local_mask].astype(np.int64)
        Selector.select(
            mesh,
            "ID",
            set(excluded_cell_ids.tolist()),
            type="elem",
            mode="s",
        )
        Selector.set_value_on_selection(mesh, "PLASMA", 0, type="elem")
        Selector.nsle(mesh, mode="u")

    def check_limiter(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)

        psi_array = mesh.GetPointData().GetArray("PSI")

        selection_array = mesh.GetPointData().GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        psi = np.asarray(vtknp.vtk_to_numpy(psi_array), dtype=float).reshape(-1)
        psi_selected = psi[selected_ids]
        limiter_local_index = int(np.argmax(psi_selected))
        limiter_index = int(selected_ids[limiter_local_index])
        psi_limiter = float(psi[limiter_index])

        points = mesh.GetPoints()
        xyz = vtknp.vtk_to_numpy(points.GetData())
        limiter_point = tuple(float(v) for v in xyz[limiter_index])

        self.status["last_psi_limiter"] = psi_limiter
        self.status["limiter_point"] = limiter_point

    def check_X_point(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)
        point_data = mesh.GetPointData()

        psi_array = point_data.GetArray("PSI")

        b_array = point_data.GetArray("B")
        b = np.asarray(vtknp.vtk_to_numpy(b_array), dtype=float)
        bx = b[:, 0]
        by = b[:, 1]

        axis = self.status["magnetic_axis"]
        axis_x = float(axis[0])
        axis_y = float(axis[1])
        min_distance = float(self.input_parameters["min_distance"])

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsle(mesh)
        selection_array = point_data.GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        points = mesh.GetPoints()
        xyz = vtknp.vtk_to_numpy(points.GetData())
        dx = xyz[selected_ids, 0] - axis_x
        dy = xyz[selected_ids, 1] - axis_y
        dist = np.sqrt(dx * dx + dy * dy)
        candidate_ids = selected_ids[dist > min_distance]

        b_norm = np.sqrt(np.square(bx[candidate_ids]) + np.square(by[candidate_ids]))
        min_local_idx = int(np.argmin(b_norm))
        x_point_node_id = int(candidate_ids[min_local_idx])

        psi = np.asarray(vtknp.vtk_to_numpy(psi_array), dtype=float).reshape(-1)
        psi_x = float(psi[x_point_node_id])

        self.status["last_psi_x_point"] = psi_x
        self.status["x_point"] = tuple(float(v) for v in xyz[x_point_node_id])

    def set_b_target_axis(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)

        point_data = mesh.GetPointData()
        b_array = point_data.GetArray("B")
        b = np.asarray(vtknp.vtk_to_numpy(b_array), dtype=float)

        target_point = (
            float(self.input_parameters["target_central_axis_x"]),
            float(self.input_parameters["target_central_axis_y"]),
            float(self.input_parameters["target_central_axis_z"]),
        )

        locator = vtk.vtkCellLocator()
        locator.SetDataSet(mesh)
        locator.BuildLocator()
        cell_id = int(locator.FindCell(target_point))
        if cell_id < 0:
            raise ValueError(
                "Target magnetic axis point is not inside the mesh: "
                f"{target_point}."
            )

        cell = mesh.GetCell(cell_id)
        n_points = int(cell.GetNumberOfPoints())
        if n_points <= 0:
            raise ValueError(f"Cell {cell_id} has no points; cannot interpolate B.")

        closest_point = [0.0, 0.0, 0.0]
        pcoords = [0.0, 0.0, 0.0]
        sub_id = vtk.reference(0)
        dist2 = vtk.reference(0.0)
        weights = [0.0] * n_points
        inside = int(
            cell.EvaluatePosition(
                target_point,
                closest_point,
                sub_id,
                pcoords,
                dist2,
                weights,
            )
        )
        if inside < 0:
            raise ValueError(
                "Could not evaluate interpolation weights for target magnetic axis "
                f"point {target_point} in cell {cell_id}."
            )

        point_ids = cell.GetPointIds()
        b_target_x = 0.0
        b_target_y = 0.0
        for local_idx, weight in enumerate(weights):
            pid = int(point_ids.GetId(local_idx))
            b_target_x += float(weight) * float(b[pid, 0])
            b_target_y += float(weight) * float(b[pid, 1])

        target = (float(b_target_x), float(b_target_y))
        self.status["b_target_axis"] = target

    def compute_correction_currents(self) -> None:
        bradx = float(self.input_parameters["bradx_value"])
        brady = float(self.input_parameters["brady_value"])
        bverx = float(self.input_parameters["bverx_value"])
        bvery = float(self.input_parameters["bvery_value"])
        cpgain = float(self.input_parameters["cp_gain"])
        btx, bty = self.status["b_target_axis"]

        det = bradx * bvery - bverx * brady
        if abs(det) <= np.finfo(float).eps:
            correction = (0.0, 0.0)
        else:
            dira = -(bverx * bty - bvery * btx) / det
            dive = -(brady * btx - bradx * bty) / det
            correction = (-cpgain * dira, -cpgain * dive)

        self.status["radial_correction_current"] = float(correction[0])
        self.status["vertical_correction_current"] = float(correction[1])

    def apply_control_currents(
        self,
        mesh_name: Optional[str] = None,
    ) -> None:
        mesh = self._get_mesh(mesh_name)

        cell_data = mesh.GetCellData()
        js_array = cell_data.GetArray("Js")
        js_np = np.asarray(vtknp.vtk_to_numpy(js_array), dtype=float)

        pfc2_array = cell_data.GetArray("PFC_2")
        pfc3_array = cell_data.GetArray("PFC_3")
        pfc4_array = cell_data.GetArray("PFC_4")
        pfc5_array = cell_data.GetArray("PFC_5")

        pfc2_mask = np.asarray(vtknp.vtk_to_numpy(pfc2_array) == 1, dtype=bool)
        pfc3_mask = np.asarray(vtknp.vtk_to_numpy(pfc3_array) == 1, dtype=bool)
        pfc4_mask = np.asarray(vtknp.vtk_to_numpy(pfc4_array) == 1, dtype=bool)
        pfc5_mask = np.asarray(vtknp.vtk_to_numpy(pfc5_array) == 1, dtype=bool)

        selected_mask = pfc2_mask | pfc3_mask | pfc4_mask | pfc5_mask
        selected_ids = np.flatnonzero(selected_mask)
        radial_current = float(self.status["radial_correction_current"])
        vertical_current = float(self.status["vertical_correction_current"])

        real_constant_array = cell_data.GetArray("CONSTANT")
        real_constant_ids = np.asarray(vtknp.vtk_to_numpy(real_constant_array), dtype=int)
        constants_df = self.apec.constants
        id_col = "CONSTANTS_ID"
        area_col = "CS"

        # Keep old semantics for duplicate IDs: the last occurrence wins.
        area_by_constant = constants_df.groupby(id_col, sort=False)[area_col].last()
        selected_rc_ids = pd.Series(real_constant_ids[selected_ids], index=selected_ids)
        selected_areas = selected_rc_ids.map(area_by_constant)
        if selected_areas.isna().any():
            missing_ids = np.unique(
                selected_rc_ids[selected_areas.isna()].to_numpy(dtype=int)
            )
            raise ValueError(
                f"Real constant IDs not found in constants table: {missing_ids.tolist()}."
            )
        selected_areas_np = selected_areas.to_numpy(dtype=float)
        if np.any(np.abs(selected_areas_np) <= np.finfo(float).eps):
            bad_ids = np.unique(
                selected_rc_ids[np.abs(selected_areas_np) <= np.finfo(float).eps]
                .to_numpy(dtype=int)
            )
            raise ValueError(
                f"Real constant IDs with invalid area (CS=0): {bad_ids.tolist()}."
            )

        base_js_z = self._old_control_current_density_base
        if base_js_z is None or base_js_z.shape[0] != js_np.shape[0]:
            base_js_z = np.asarray(js_np[:, 2], dtype=float).copy()
            self._old_control_current_density_base = base_js_z

        corr_current = np.zeros(selected_ids.shape[0], dtype=float)
        sel_pfc2 = pfc2_mask[selected_ids]
        sel_pfc3 = pfc3_mask[selected_ids]
        sel_pfc4 = pfc4_mask[selected_ids]
        sel_pfc5 = pfc5_mask[selected_ids]
        # JPF2=(xirad+IPF2)/APF2
        corr_current[sel_pfc2] = radial_current
        # JPF3=(xiver+IPF3)/APF3
        corr_current[sel_pfc3] = vertical_current
        # JPF4=(xiver+IPF4)/APF4
        corr_current[sel_pfc4] = vertical_current
        # JPF5=(-xirad+IPF5)/APF5
        corr_current[sel_pfc5] = -radial_current

        correction_current_density = corr_current / selected_areas_np
        js_np[selected_ids, 2] = base_js_z[selected_ids] + correction_current_density

        js_array.Modified()

    def find_separatrix(self, mesh_name: Optional[str] = None) -> None:
        self._dispatch_step(
            self.STEP_FIND_SEPARATRIX,
            self._find_separatrix_default,
            mesh_name=mesh_name,
        )

    def _find_separatrix_default(self, mesh_name: Optional[str] = None) -> None:
        self.reset_plasma_components(mesh_name=mesh_name)
        self.set_psi_magnetic_axis(mesh_name=mesh_name)
        self.remove_plasma_out_of_X_point(mesh_name=mesh_name)
        self.check_limiter(mesh_name=mesh_name)
        self.check_X_point(mesh_name=mesh_name)
        psi_limiter = float(self.status["last_psi_limiter"])
        psi_x_point = float(self.status["last_psi_x_point"])
        if psi_limiter > psi_x_point:
            self.status["psi_min"] = psi_limiter
            self.status["separatrix_configuration"] = "limiter"
        else:
            self.status["psi_min"] = psi_x_point
            self.status["separatrix_configuration"] = "X-point"
        self.identify_plasma_region(mesh_name=mesh_name)

    def compute_jplasma_gradshafranov(
        self, mesh_name: Optional[str] = None
    ) -> None:
        self._dispatch_step(
            self.STEP_COMPUTE_JPLASMA,
            self._compute_jplasma_gradshafranov_default,
            mesh_name=mesh_name,
        )

    def _compute_jplasma_gradshafranov_default(
        self, mesh_name: Optional[str] = None
    ) -> None:
        mesh = self._get_mesh(mesh_name)

        cell_data = mesh.GetCellData()

        js_array = cell_data.GetArray("Js")
        js_np = np.asarray(vtknp.vtk_to_numpy(js_array), dtype=float)

        psie_array = cell_data.GetArray("PSIe")
        psi_elem = np.asarray(vtknp.vtk_to_numpy(psie_array), dtype=float).reshape(-1)

        be_array = cell_data.GetArray("Be")
        be_np = np.asarray(vtknp.vtk_to_numpy(be_array), dtype=float)
        bx = np.asarray(be_np[:, 0], dtype=float)
        by = np.asarray(be_np[:, 1], dtype=float)

        volume_array = cell_data.GetArray("cell_volume")
        cell_volume = np.asarray(vtknp.vtk_to_numpy(volume_array), dtype=float).reshape(-1)

        centroids_array = cell_data.GetArray("cell_centroids")
        centroids = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
        elem_centroid_rad = np.asarray(centroids[:, 0], dtype=float)

        psi_magnetic_axis = float(self.status["psi_magnetic_axis"])
        psi_min = float(self.status["psi_min"])
        psi_den = psi_min - psi_magnetic_axis
        if abs(psi_den) <= np.finfo(float).eps:
            self.status["total_plasma_current_before_scaling"] = 0.0
            self.status["plasma_current_scaling_lambda"] = 0.0
            self.status["plasma_element_count"] = 0
            return

        beta_profile = float(self.input_parameters["beta_p"])
        alpha_mp = float(self.input_parameters["alpha_mp"])
        alpha_np = float(self.input_parameters["alpha_np"])
        r0 = float(self.input_parameters["R0"])
        af0 = float(self.input_parameters["AF0"])
        target_plasma_current = float(self.input_parameters["plasma_current"])

        Selector.select(mesh, "PLASMA", 1, type="elem", mode="s")
        selection_array = cell_data.GetArray("selection_cell")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        psi_segn = (psi_elem[selected_ids] - psi_magnetic_axis) / psi_den
        one_minus = 1.0 - np.power(psi_segn, alpha_mp)
        pprime_h2 = np.power(one_minus, alpha_np)
        grad_shafranov_current_density = pprime_h2 * (
            beta_profile * elem_centroid_rad[selected_ids] / r0
            + (1.0 - beta_profile) * (r0 / elem_centroid_rad[selected_ids])
        )

        total_plasma_current = float(
            np.sum(grad_shafranov_current_density * cell_volume[selected_ids])
        )
        lambda_scale = 0.0
        if abs(total_plasma_current) > np.finfo(float).eps:
            lambda_scale = target_plasma_current / total_plasma_current

        js_toroidal = grad_shafranov_current_density * lambda_scale

        mu_0 = 4.0 * np.pi * 1e-7
        ffprime_h2 = np.power(one_minus, alpha_np + 1.0)
        c_const = (
            2.0
            * mu_0
            * r0
            * (1.0 - beta_profile)
            * (psi_magnetic_axis - psi_min)
            * lambda_scale
            / (alpha_np + 1.0)
        )
        f_value = np.sqrt(np.maximum(c_const * ffprime_h2 + af0 * af0, 0.0))
        ffprime = pprime_h2 * (1.0 - beta_profile) * r0
        fprime = (ffprime * lambda_scale) / f_value

        js_radial = fprime * bx[selected_ids]
        js_vertical = fprime * by[selected_ids]

        js_np[selected_ids, 0] = js_radial
        js_np[selected_ids, 1] = js_vertical
        js_np[selected_ids, 2] = js_toroidal
        js_array.Modified()

        self.status["total_plasma_current_before_scaling"] = total_plasma_current
        self.status["plasma_current_scaling_lambda"] = float(lambda_scale)
        self.status["plasma_element_count"] = int(selected_ids.size)

    def control_currents(self, mesh_name: Optional[str] = None) -> None:
        self._dispatch_step(
            self.STEP_CONTROL_CURRENTS,
            self._control_currents_default,
            mesh_name=mesh_name,
        )

    def _control_currents_default(self, mesh_name: Optional[str] = None) -> None:
        self.set_b_target_axis(mesh_name=mesh_name)
        self.compute_correction_currents()
        self.apply_control_currents(mesh_name=mesh_name)

    def check_convergence(self, mesh_name: Optional[str] = None) -> bool:
        result = self._dispatch_step(
            self.STEP_CHECK_CONVERGENCE,
            self._check_convergence_default,
            mesh_name=mesh_name,
        )
        return bool(result)

    def _check_convergence_default(
        self, mesh_name: Optional[str] = None
    ) -> bool:
        mesh = self._get_mesh(mesh_name)

        point_data = mesh.GetPointData()
        psi_array = point_data.GetArray("PSI")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsle(mesh)
        selection_array = point_data.GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        psi = np.asarray(vtknp.vtk_to_numpy(psi_array), dtype=float).reshape(-1)
        current_flux_vector = psi[selected_ids].copy()

        tolerance = float(self.input_parameters["tolerance_value"])
        self.status["flux_tolerance"] = tolerance

        if self._first_absolute_iteration or self._previous_flux_vector is None:
            self._previous_flux_vector = current_flux_vector.copy()
            self.status["flux_norm"] = np.nan
            self.status["converged"] = False
            return False

        if self._previous_flux_vector.size != current_flux_vector.size:
            self._previous_flux_vector = current_flux_vector.copy()
            self.status["flux_norm"] = np.nan
            self.status["converged"] = False
            return False

        flux_norm = float(np.max(np.abs(current_flux_vector - self._previous_flux_vector)))
        self._previous_flux_vector = current_flux_vector.copy()

        converged = bool(flux_norm <= tolerance)
        self.status["flux_norm"] = flux_norm
        self.status["converged"] = converged
        return converged

    # Backward-compatible alias.
    def check_poloidal_flux_convergence(
        self, mesh_name: Optional[str] = None
    ) -> bool:
        return self.check_convergence(mesh_name=mesh_name)

    def toroidal_average(
        self,
        equilibrium_mesh_name: Optional[str] = None,
        evolution_mesh_name: Optional[str] = None,
    ) -> None:
        from .coordinate_system import CylindricalCoordinateSystem

        equilibrium_name = str(
            equilibrium_mesh_name
            or self.input_parameters["equilibrium_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        evolution_name = str(
            evolution_mesh_name
            or self.input_parameters["evolution_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        if not equilibrium_name or not evolution_name:
            raise ValueError(
                "toroidal_average requires equilibrium and evolution mesh names."
            )

        equilibrium_mesh = self._get_mesh(equilibrium_name)
        evolution_mesh = self._get_mesh(evolution_name)
        node_slices, elem_slices = self._require_toroidal_slices()

        eq_point_data = equilibrium_mesh.GetPointData()
        evol_point_data = evolution_mesh.GetPointData()
        eq_cell_data = equilibrium_mesh.GetCellData()
        evol_cell_data = evolution_mesh.GetCellData()

        eq_psi_array = eq_point_data.GetArray("PSI")
        evol_psi_array = evol_point_data.GetArray("PSI")
        eq_be_array = eq_cell_data.GetArray("Be")
        evol_be_array = evol_cell_data.GetArray("Be")
        eq_psie_array = eq_cell_data.GetArray("PSIe")
        evol_psie_array = evol_cell_data.GetArray("PSIe")
        eq_js_array = eq_cell_data.GetArray("Js")
        evol_jt_array = evol_cell_data.GetArray("Jt")
        if (
            eq_psi_array is None
            or evol_psi_array is None
            or eq_be_array is None
            or evol_be_array is None
            or eq_psie_array is None
            or evol_psie_array is None
            or eq_js_array is None
            or evol_jt_array is None
        ):
            raise ValueError(
                "Required arrays for toroidal_average are missing "
                "(PSI, Be, PSIe, Js, Jt)."
            )

        eq_psi = np.asarray(vtknp.vtk_to_numpy(eq_psi_array), dtype=float).reshape(-1)
        evol_psi = np.asarray(vtknp.vtk_to_numpy(evol_psi_array), dtype=float).reshape(-1)

        src_nodes = np.asarray(node_slices[:, 0], dtype=np.int64)
        dst_nodes = np.asarray(node_slices[:, 1:], dtype=np.int64)
        node_ids, node_values = self._average_from_slice_targets(
            source_values=evol_psi,
            source_ids=src_nodes,
            target_ids_matrix=dst_nodes,
            source_size=eq_psi.shape[0],
        )
        eq_psi[node_ids] = node_values
        eq_psi_array.Modified()

        cylindrical = CylindricalCoordinateSystem()
        ModelEdit.map_vector_field_to_local_components(
            evolution_mesh,
            "Jt",
            cylindrical,
            location="elem",
        )
        ModelEdit.map_vector_field_to_local_components(
            equilibrium_mesh,
            "Be",
            cylindrical,
            location="elem",
        )
        ModelEdit.map_vector_field_to_local_components(
            equilibrium_mesh,
            "Js",
            cylindrical,
            location="elem",
        )

        eq_be = np.asarray(vtknp.vtk_to_numpy(eq_be_array), dtype=float)
        evol_be = np.asarray(vtknp.vtk_to_numpy(evol_be_array), dtype=float)
        eq_psie = np.asarray(vtknp.vtk_to_numpy(eq_psie_array), dtype=float).reshape(-1)
        evol_psie = np.asarray(vtknp.vtk_to_numpy(evol_psie_array), dtype=float).reshape(-1)
        eq_js = np.asarray(vtknp.vtk_to_numpy(eq_js_array), dtype=float)
        evol_jt = np.asarray(vtknp.vtk_to_numpy(evol_jt_array), dtype=float)

        eq_dim = int(equilibrium_mesh.GetMaxSpatialDimension())
        evol_dim = int(evolution_mesh.GetMaxSpatialDimension())
        eq_main, eq_rad, eq_ver = self._vector_component_indices(eq_dim)
        evol_main, evol_rad, evol_ver = self._vector_component_indices(evol_dim)

        src_elems = np.asarray(elem_slices[:, 0], dtype=np.int64)
        dst_elems = np.asarray(elem_slices[:, 1:], dtype=np.int64)

        src_ids, rad_values = self._average_from_slice_targets(
            source_values=evol_be[:, evol_rad],
            source_ids=src_elems,
            target_ids_matrix=dst_elems,
            source_size=eq_be.shape[0],
        )
        eq_be[src_ids, eq_rad] = rad_values

        src_ids, ver_values = self._average_from_slice_targets(
            source_values=evol_be[:, evol_ver],
            source_ids=src_elems,
            target_ids_matrix=dst_elems,
            source_size=eq_be.shape[0],
        )
        eq_be[src_ids, eq_ver] = ver_values

        src_ids, psie_values = self._average_from_slice_targets(
            source_values=evol_psie,
            source_ids=src_elems,
            target_ids_matrix=dst_elems,
            source_size=eq_psie.shape[0],
        )
        eq_psie[src_ids] = psie_values

        coils_array = eq_cell_data.GetArray("COILS")
        if coils_array is not None:
            coils_mask = np.asarray(vtknp.vtk_to_numpy(coils_array) == 1, dtype=bool)
            src_ids, js_main_values = self._average_from_slice_targets(
                source_values=evol_jt[:, evol_main],
                source_ids=src_elems,
                target_ids_matrix=dst_elems,
                source_size=eq_js.shape[0],
            )
            if src_ids.size > 0:
                valid_coils = coils_mask[src_ids]
                eq_js[src_ids[valid_coils], eq_main] = js_main_values[valid_coils]

        eq_be_array.Modified()
        eq_psie_array.Modified()
        eq_js_array.Modified()

        ModelEdit.map_vector_field_to_global_components(
            evolution_mesh,
            "Jt",
            cylindrical,
            location="elem",
        )
        ModelEdit.map_vector_field_to_global_components(
            equilibrium_mesh,
            "Be",
            cylindrical,
            location="elem",
        )
        ModelEdit.map_vector_field_to_global_components(
            equilibrium_mesh,
            "Js",
            cylindrical,
            location="elem",
        )

    def toroidal_average_psi(
        self,
        equilibrium_mesh_name: Optional[str] = None,
        static_mesh_name: Optional[str] = None,
    ) -> None:
        equilibrium_name = str(
            equilibrium_mesh_name
            or self.input_parameters["equilibrium_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        static_name = str(
            static_mesh_name
            or self.input_parameters["static_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        if not equilibrium_name or not static_name:
            raise ValueError(
                "toroidal_average_psi requires equilibrium and static mesh names."
            )

        equilibrium_mesh = self._get_mesh(equilibrium_name)
        static_mesh = self._get_mesh(static_name)
        node_slices, _ = self._require_toroidal_slices()

        eq_psi_array = equilibrium_mesh.GetPointData().GetArray("PSI")
        static_psi_array = static_mesh.GetPointData().GetArray("PSI")
        if eq_psi_array is None or static_psi_array is None:
            raise ValueError("PointData array 'PSI' not found on equilibrium/static mesh.")

        eq_psi = np.asarray(vtknp.vtk_to_numpy(eq_psi_array), dtype=float).reshape(-1)
        static_psi = np.asarray(vtknp.vtk_to_numpy(static_psi_array), dtype=float).reshape(-1)

        src_nodes = np.asarray(node_slices[:, 0], dtype=np.int64)
        dst_nodes = np.asarray(node_slices[:, 1:], dtype=np.int64)
        node_ids, node_values = self._average_from_slice_targets(
            source_values=static_psi,
            source_ids=src_nodes,
            target_ids_matrix=dst_nodes,
            source_size=eq_psi.shape[0],
        )
        eq_psi[node_ids] = node_values
        eq_psi_array.Modified()

    def toroidal_projection_plasma_current(
        self,
        static_mesh_name: Optional[str] = None,
        equilibrium_mesh_name: Optional[str] = None,
    ) -> None:
        from .coordinate_system import CylindricalCoordinateSystem

        static_name = str(
            static_mesh_name
            or self.input_parameters["static_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        equilibrium_name = str(
            equilibrium_mesh_name
            or self.input_parameters["equilibrium_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        if not static_name or not equilibrium_name:
            raise ValueError(
                "toroidal_projection_plasma_current requires static and equilibrium mesh names."
            )

        static_mesh = self._get_mesh(static_name)
        equilibrium_mesh = self._get_mesh(equilibrium_name)
        _, elem_slices = self._require_toroidal_slices()

        static_js_array = static_mesh.GetCellData().GetArray("Js")
        equilibrium_js_array = equilibrium_mesh.GetCellData().GetArray("Js")
        if static_js_array is None or equilibrium_js_array is None:
            raise ValueError("CellData array 'Js' not found on static/equilibrium mesh.")

        cylindrical = CylindricalCoordinateSystem()
        ModelEdit.map_vector_field_to_local_components(
            static_mesh,
            "Js",
            cylindrical,
            location="elem",
        )

        static_js = np.asarray(vtknp.vtk_to_numpy(static_js_array), dtype=float)
        equilibrium_js = np.asarray(vtknp.vtk_to_numpy(equilibrium_js_array), dtype=float)

        static_dim = int(static_mesh.GetMaxSpatialDimension())
        equilibrium_dim = int(equilibrium_mesh.GetMaxSpatialDimension())
        static_main, static_c1, static_c2 = self._vector_component_indices(static_dim)
        equilibrium_main, equilibrium_c1, equilibrium_c2 = self._vector_component_indices(
            equilibrium_dim
        )

        src_elems = np.asarray(elem_slices[:, 0], dtype=np.int64)
        dst_elems = np.asarray(elem_slices[:, 1:], dtype=np.int64)

        dst_ids, main_values = self._broadcast_rows_to_slice_targets(
            source_values=equilibrium_js[:, equilibrium_main],
            source_ids=src_elems,
            target_ids_matrix=dst_elems,
            target_size=static_js.shape[0],
        )
        static_js[dst_ids, static_main] = main_values

        if static_dim == 3:
            dst_ids, c1_values = self._broadcast_rows_to_slice_targets(
                source_values=equilibrium_js[:, equilibrium_c1],
                source_ids=src_elems,
                target_ids_matrix=dst_elems,
                target_size=static_js.shape[0],
            )
            static_js[dst_ids, static_c1] = c1_values

            dst_ids, c2_values = self._broadcast_rows_to_slice_targets(
                source_values=equilibrium_js[:, equilibrium_c2],
                source_ids=src_elems,
                target_ids_matrix=dst_elems,
                target_size=static_js.shape[0],
            )
            static_js[dst_ids, static_c2] = c2_values

        static_js_array.Modified()
        ModelEdit.map_vector_field_to_global_components(
            static_mesh,
            "Js",
            cylindrical,
            location="elem",
        )

    def impressed_eddy_currents_coils_currents(
        self,
        evolution_mesh_name: Optional[str] = None,
        static_mesh_name: Optional[str] = None,
    ) -> None:
        evolution_name = str(
            evolution_mesh_name
            or self.input_parameters["evolution_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        static_name = str(
            static_mesh_name
            or self.input_parameters["static_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        if not evolution_name or not static_name:
            raise ValueError(
                "impressed_eddy_currents_coils_currents requires evolution and static mesh names."
            )

        evolution_mesh = self._get_mesh(evolution_name)
        static_mesh = self._get_mesh(static_name)

        evol_cell_data = evolution_mesh.GetCellData()
        static_cell_data = static_mesh.GetCellData()
        evol_jt_array = evol_cell_data.GetArray("Jt")
        static_js_array = static_cell_data.GetArray("Js")
        if evol_jt_array is None or static_js_array is None:
            raise ValueError("CellData array 'Jt'/'Js' not found on evolution/static mesh.")

        evol_jt = np.asarray(vtknp.vtk_to_numpy(evol_jt_array), dtype=float)
        static_js = np.asarray(vtknp.vtk_to_numpy(static_js_array), dtype=float)

        evol_passive = evol_cell_data.GetArray("PASSIVE_STRUCTURES")
        evol_coils = evol_cell_data.GetArray("COILS")
        static_passive = static_cell_data.GetArray("PASSIVE_STRUCTURES")
        static_coils = static_cell_data.GetArray("COILS")

        evol_mask = np.asarray(vtknp.vtk_to_numpy(evol_passive) == 1, dtype=bool) | np.asarray(
            vtknp.vtk_to_numpy(evol_coils) == 1, dtype=bool
        )
        static_mask = np.asarray(vtknp.vtk_to_numpy(static_passive) == 1, dtype=bool) | np.asarray(
            vtknp.vtk_to_numpy(static_coils) == 1, dtype=bool
        )
        shared_ids = np.flatnonzero(evol_mask & static_mask)

        static_js[shared_ids, 2] = evol_jt[shared_ids, 2]

        if int(evolution_mesh.GetMaxSpatialDimension()) == 3:
            static_js[shared_ids, 0] = evol_jt[shared_ids, 0]
            static_js[shared_ids, 1] = evol_jt[shared_ids, 1]

        static_js_array.Modified()

    def set_plasma_current(
        self,
        static_mesh_name: Optional[str] = None,
        evolution_mesh_name: Optional[str] = None,
    ) -> None:
        static_name = str(
            static_mesh_name
            or self.input_parameters["static_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        evolution_name = str(
            evolution_mesh_name
            or self.input_parameters["evolution_mesh"]
            or self.apec.active_mesh_name
            or ""
        )
        if not static_name or not evolution_name:
            raise ValueError(
                "set_plasma_current requires static and evolution mesh names."
            )

        static_mesh = self._get_mesh(static_name)
        evolution_mesh = self._get_mesh(evolution_name)

        static_cell_data = static_mesh.GetCellData()
        evolution_cell_data = evolution_mesh.GetCellData()
        static_js_array = static_cell_data.GetArray("Js")
        evolution_js_array = evolution_cell_data.GetArray("Js")
        vacuum_static = static_cell_data.GetArray("VACUUM_CHAMBER")
        vacuum_evolution = evolution_cell_data.GetArray("VACUUM_CHAMBER")

        static_js = np.asarray(vtknp.vtk_to_numpy(static_js_array), dtype=float)
        evolution_js = np.asarray(vtknp.vtk_to_numpy(evolution_js_array), dtype=float)

        static_mask = np.asarray(vtknp.vtk_to_numpy(vacuum_static) == 1, dtype=bool)
        evolution_mask = np.asarray(vtknp.vtk_to_numpy(vacuum_evolution) == 1, dtype=bool)
        shared_ids = np.flatnonzero(static_mask & evolution_mask)

        evolution_js[shared_ids, 2] = static_js[shared_ids, 2]

        if int(evolution_mesh.GetMaxSpatialDimension()) == 3:
            evolution_js[shared_ids, 0] = static_js[shared_ids, 0]
            evolution_js[shared_ids, 1] = static_js[shared_ids, 1]

        evolution_js_array.Modified()

    def set_poloidal_currents(
        self,
        static_mesh_name: Optional[str] = None,
    ) -> None:
        # Placeholder for old static-mesh poloidal current setup.
        _ = static_mesh_name

    def append_iteration_diagnostics(
        self,
        converged: bool,
        mesh_name: Optional[str] = None,
    ) -> None:
        self.status["converged"] = bool(converged)

        if self.apec._world_rank() != 0:
            return

        mesh_key = self._diagnostics_mesh_key(mesh_name)
        row = self._build_diagnostics_row(converged)
        row_df = pd.DataFrame([row])

        if mesh_key in self._diagnostics_frames:
            diagnostics_df = pd.concat(
                [self._diagnostics_frames[mesh_key], row_df],
                ignore_index=True,
            )
        else:
            diagnostics_df = row_df

        self._diagnostics_frames[mesh_key] = diagnostics_df
        diagnostics_df.to_csv(self._diagnostics_csv_path(mesh_name), index=False)

    def calculate_poloidal_flux(self, mesh_name: Optional[str] = None) -> None:
        from .coordinate_system import CylindricalCoordinateSystem
        from .cyclic import Cyclic

        mesh = self._get_mesh(mesh_name)
        ModelEdit.element_solution(mesh, "A")
        Cyclic.poloidal_flux(
            mesh=mesh,
            coordinate_system=CylindricalCoordinateSystem(),
        )

    def equilibrium_solve(
        self,
        nodal_solution: bool = True,
    ) -> bool:
        mesh_name = str(
            self.input_parameters["equilibrium_mesh"] or self.apec.active_mesh_name or ""
        )
        if not mesh_name:
            raise RuntimeError(
                "equilibrium_solve requires an active mesh or input_parameters['equilibrium_mesh']."
            )

        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )
        self.status["converged"] = False

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1

            if not self._first_absolute_iteration:
                self.find_separatrix(mesh_name=mesh_name)
                self.compute_jplasma_gradshafranov(mesh_name=mesh_name)
                self.control_currents(mesh_name=mesh_name)

            self.apec.set_active_mesh(mesh_name)
            self.apec.solve(nodal_solution=nodal_solution)
            self.calculate_poloidal_flux(mesh_name=mesh_name)

            converged = self.check_convergence(mesh_name=mesh_name)
            self.append_iteration_diagnostics(converged, mesh_name=mesh_name)
            if converged:
                self.status["converged"] = True
                return True
            # self.apec.write_output()
            self._first_absolute_iteration = False

        self.status["converged"] = False
        return False

    def evolution_solve(
        self,
        nodal_solution: bool = False,
    ) -> bool:
        equilibrium_mesh = str(
            self.input_parameters["equilibrium_mesh"] or self.apec.active_mesh_name or ""
        )
        evolution_mesh = str(
            self.input_parameters["evolution_mesh"] or self.apec.active_mesh_name or ""
        )
        static_mesh = str(self.input_parameters["static_mesh"] or equilibrium_mesh)

        if not evolution_mesh:
            raise RuntimeError(
                "evolution_solve requires an active mesh or input_parameters['evolution_mesh']."
            )

        if self._first_absolute_iteration:
            self.set_toroidal_slices()

        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1

            if not self._first_absolute_iteration:
                self.find_separatrix(mesh_name=equilibrium_mesh)
                self.compute_jplasma_gradshafranov(mesh_name=equilibrium_mesh)

                self.apec.set_active_mesh(static_mesh)
                self.toroidal_projection_plasma_current(
                    static_mesh_name=static_mesh,
                    equilibrium_mesh_name=equilibrium_mesh,
                )
                self.impressed_eddy_currents_coils_currents(
                    evolution_mesh_name=evolution_mesh,
                    static_mesh_name=static_mesh,
                )

                self.apec.set_active_mesh(static_mesh)
                self.apec.time_step["transient_effects"] = False
                self.apec.solve(nodal_solution=nodal_solution)
                self.calculate_poloidal_flux(mesh_name=static_mesh)
                # self.apec.write_output()

                self.toroidal_average_psi(
                    equilibrium_mesh_name=equilibrium_mesh,
                    static_mesh_name=static_mesh,
                )

                converged = self.check_convergence(mesh_name=equilibrium_mesh)
                self.append_iteration_diagnostics(
                    converged,
                    mesh_name=equilibrium_mesh,
                )
                if converged:
                    self.status["converged"] = True
                    return True

                self.set_plasma_current(
                    static_mesh_name=static_mesh,
                    evolution_mesh_name=evolution_mesh,
                )
                self.apec.time_step["transient_effects"] = True


            self.apec.set_active_mesh(evolution_mesh)
            self.apec.solve(nodal_solution=nodal_solution)
            self.calculate_poloidal_flux(mesh_name=evolution_mesh)
            # self.apec.write_output()
            self.toroidal_average(
                equilibrium_mesh_name=equilibrium_mesh,
                evolution_mesh_name=evolution_mesh,
            )

            if self._first_absolute_iteration:
                self.check_convergence(mesh_name=equilibrium_mesh)
                self._first_absolute_iteration = False
                self.status["converged"] = True
                return True

        self.status["converged"] = False
        return False

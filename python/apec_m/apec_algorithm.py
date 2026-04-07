from __future__ import annotations

import copy
import os
from typing import Any, Callable, Dict, Optional, TYPE_CHECKING, Tuple

import numpy as np
import pandas as pd
import vtk
from .selector import Selector
from .model_edit import ModelEdit

import vtk.util.numpy_support as vtknp

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    from .apec_m import APEC_M


DEFAULT_APEC_ALGORITHM_INPUT_PARAMETERS: Dict[str, Any] = {
    "target_central_axis_x": 0.0,
    "target_central_axis_y": 0.0,
    "target_central_axis_z": 0.0,
    "R0": 0.0,
    "max_picard_iterations": 1,
    "min_distance": 0.3,
    "cp_gain": 0.0,
    "tolerance_value": 1e-6,
    "bradx_value": 0.0,
    "brady_value": 0.0,
    "bverx_value": 0.0,
    "bvery_value": 0.0,
    "alpha_mp": 1.0001,
    "alpha_np": 0.83197,
    "plasma_current": 0.0,
    "beta_p": 0.0,
    "AF0": 5.3 * 6.2,
    "print_full_output": 0,
    "equilibrium_mesh": "",
    "static_mesh": "",
    "evolution_mesh": "",
    "n_slices": 1,
    "toroidal_slices_nodes": np.zeros((0, 0), dtype=np.int64),
    "toroidal_slices_elems": np.zeros((0, 0), dtype=np.int64),
}


DEFAULT_APEC_ALGORITHM_STATUS: Dict[str, Any] = {
    "psi_magnetic_axis": 0.0,
    "psi_max": 0.0,
    "psi_min": 0.0,
    "magnetic_axis": (0.0, 0.0, 0.0),
    "b_target_axis": (0.0, 0.0),
    "last_psi_limiter": 0.0,
    "last_psi_x_point": 0.0,
    "separatrix_configuration": "uninitialized",
    "limiter_point": (np.nan, np.nan, np.nan),
    "x_point": (np.nan, np.nan, np.nan),
    "flux_norm": np.nan,
    "flux_tolerance": np.nan,
    "radial_correction_current": 0.0,
    "vertical_correction_current": 0.0,
    "total_plasma_current_before_scaling": 0.0,
    "plasma_current_scaling_lambda": 0.0,
    "plasma_element_count": 0,
    "picard_iteration": 0,
    "converged": False,
}


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
        self._old_control_current_density_base: Dict[str, np.ndarray] = {}
        self._diagnostics_frames: Dict[str, pd.DataFrame] = {}
        self._step_overrides: Dict[str, Callable[..., Any]] = {}
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
        self._old_control_current_density_base.clear()
        self._diagnostics_frames.clear()

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
        if psi_array is None:
            raise ValueError("PointData array 'PSI' not found on the active mesh.")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsle(mesh)
        selection_array = mesh.GetPointData().GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        if selected_ids.size == 0:
            raise ValueError(
                "No nodes selected with VACUUM_CHAMBER == 1; cannot find magnetic axis."
            )

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

        if cell_data.GetArray("VACUUM_CHAMBER") is None:
            raise ValueError("CellData array 'VACUUM_CHAMBER' not found on the mesh.")
        if cell_data.GetArray("PLASMA") is None:
            raise ValueError("CellData array 'PLASMA' not found on the mesh.")

        js_array = cell_data.GetArray("Js")
        if js_array is None:
            raise ValueError("CellData array 'Js' not found on the mesh.")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.set_value_on_selection(mesh, "PLASMA", 1, type="elem")

        selection_array = cell_data.GetArray("selection_cell")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        if selected_ids.size == 0:
            return

        js_np = vtknp.vtk_to_numpy(js_array)
        if js_np.ndim == 1:
            js_np[selected_ids] = 0.0
        else:
            js_np[selected_ids, :] = 0.0
        js_array.Modified()

    def identify_plasma_region(self, mesh_name: Optional[str] = None) -> None:
        mesh = self._get_mesh(mesh_name)

        psi_array = mesh.GetPointData().GetArray("PSI")
        if psi_array is None:
            raise ValueError("PointData array 'PSI' not found on the active mesh.")

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

        vacuum_array = cell_data.GetArray("VACUUM_CHAMBER")
        if vacuum_array is None:
            raise ValueError("CellData array 'VACUUM_CHAMBER' not found on the mesh.")
        if cell_data.GetArray("PLASMA") is None:
            raise ValueError("CellData array 'PLASMA' not found on the mesh.")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsel_ext(mesh)
        selection_array = cell_data.GetArray("selection_cell")
        vacuum_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
        if not np.any(vacuum_mask):
            return
        vacuum_ids = np.flatnonzero(vacuum_mask)

        centroids_array = cell_data.GetArray("cell_centroids")
        if centroids_array is None:
            raise ValueError(
                "CellData array 'cell_centroids' not found on the mesh."
            )
        centroids = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
        if centroids.ndim != 2 or centroids.shape[1] < 2:
            raise ValueError(
                "CellData array 'cell_centroids' must have at least 2 components."
            )
        x_centroid = np.asarray(centroids[:, 0], dtype=float)
        y_centroid = np.asarray(centroids[:, 1], dtype=float)

        magnetic_axis = self.status["magnetic_axis"]
        axis_x = float(magnetic_axis[0])
        axis_y = float(magnetic_axis[1])

        plasma_current = float(self.input_parameters["plasma_current"])
        current_sign = 1.0 if plasma_current >= 0.0 else -1.0

        b_cell_array = cell_data.GetArray("Be")
        if b_cell_array is None:
            raise ValueError("CellData array 'Be' not found on the mesh.")
        b_cell = np.asarray(vtknp.vtk_to_numpy(b_cell_array), dtype=float)
        if b_cell.ndim != 2 or b_cell.shape[1] < 2:
            raise ValueError("CellData array 'Be' must have at least 2 components.")
        bx = b_cell[:, 0]
        by = b_cell[:, 1]

        b_rot = current_sign * (
            by[vacuum_ids] * (x_centroid[vacuum_ids] - axis_x)
            - bx[vacuum_ids] * (y_centroid[vacuum_ids] - axis_y)
        )
        excluded_local_mask = b_rot >= 0.0
        if not np.any(excluded_local_mask):
            return

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
        if psi_array is None:
            raise ValueError("PointData array 'PSI' not found on the active mesh.")

        # Assume here to have calle before remove_plasma_out_of_X_point, so the limiter is the node with max PSI among those connected to the vacuum chamber.

        selection_array = mesh.GetPointData().GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        if selected_ids.size == 0:
            raise ValueError(
                "No nodes selected with NEW_VACUUM_CHAMBER_BOUNDARY == 1; cannot find limiter."
            )

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
        if psi_array is None:
            raise ValueError("PointData array 'PSI' not found on the active mesh.")

        b_array = point_data.GetArray("B")
        if b_array is None:
            raise ValueError("PointData array 'B' not found on the active mesh.")
        b = np.asarray(vtknp.vtk_to_numpy(b_array), dtype=float)
        if b.ndim != 2 or b.shape[1] < 2:
            raise ValueError("PointData array 'B' must have at least 2 components.")
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
        if selected_ids.size == 0:
            raise ValueError(
                "No nodes selected with VACUUM_CHAMBER == 1; cannot evaluate X-point."
            )

        points = mesh.GetPoints()
        xyz = vtknp.vtk_to_numpy(points.GetData())
        dx = xyz[selected_ids, 0] - axis_x
        dy = xyz[selected_ids, 1] - axis_y
        dist = np.sqrt(dx * dx + dy * dy)
        candidate_ids = selected_ids[dist > min_distance]
        if candidate_ids.size == 0:
            raise ValueError(
                "No VACUUM_CHAMBER nodes found outside min_distance from magnetic axis."
            )

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
        if b_array is None:
            raise ValueError("PointData array 'B' not found on the active mesh.")

        b = np.asarray(vtknp.vtk_to_numpy(b_array), dtype=float)
        if b.ndim != 2 or b.shape[1] < 2:
            raise ValueError("PointData array 'B' must have at least 2 components.")

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
        if mesh_name is None:
            mesh = self.apec.mesh
            mesh_key = str(self.apec.active_mesh_name or "__active__")
        else:
            mesh = self.apec.get_mesh(mesh_name)
            if mesh is None:
                raise ValueError(f"Mesh '{mesh_name}' does not exist.")
            mesh_key = str(mesh_name)

        cell_data = mesh.GetCellData()
        js_array = cell_data.GetArray("Js")
        if js_array is None:
            raise ValueError("CellData array 'Js' not found on the mesh.")
        js_np = np.asarray(vtknp.vtk_to_numpy(js_array), dtype=float)
        if js_np.ndim != 2 or js_np.shape[1] < 3:
            raise ValueError("CellData array 'Js' must have at least 3 components.")

        pfc2_array = cell_data.GetArray("PFC_2")
        pfc3_array = cell_data.GetArray("PFC_3")
        pfc4_array = cell_data.GetArray("PFC_4")
        pfc5_array = cell_data.GetArray("PFC_5")
        if (
            pfc2_array is None
            or pfc3_array is None
            or pfc4_array is None
            or pfc5_array is None
        ):
            raise ValueError(
                "CellData arrays 'PFC_2', 'PFC_3', 'PFC_4', and 'PFC_5' are required."
            )

        pfc2_mask = np.asarray(vtknp.vtk_to_numpy(pfc2_array) == 1, dtype=bool)
        pfc3_mask = np.asarray(vtknp.vtk_to_numpy(pfc3_array) == 1, dtype=bool)
        pfc4_mask = np.asarray(vtknp.vtk_to_numpy(pfc4_array) == 1, dtype=bool)
        pfc5_mask = np.asarray(vtknp.vtk_to_numpy(pfc5_array) == 1, dtype=bool)

        selected_mask = pfc2_mask | pfc3_mask | pfc4_mask | pfc5_mask
        selected_ids = np.flatnonzero(selected_mask)
        radial_current = float(self.status["radial_correction_current"])
        vertical_current = float(self.status["vertical_correction_current"])
        if selected_ids.size == 0:
            return

        real_constant_array = cell_data.GetArray("CONSTANT")
        if real_constant_array is None:
            raise ValueError(
                "CellData array 'CONSTANT' not found on the mesh."
            )
        real_constant_ids = np.asarray(vtknp.vtk_to_numpy(real_constant_array), dtype=int)

        constants_df = self.apec.constants
        if constants_df is None or constants_df.empty:
            raise ValueError("APEC constants table is empty; cannot apply control currents.")

        id_col = "CONSTANTS_ID"
        area_col = "CS"
        if id_col not in constants_df.columns or area_col not in constants_df.columns:
            raise ValueError(
                "Constants table must contain ID column "
                "'CONSTANTS_ID' and area column 'CS'."
            )

        constant_lookup: Dict[int, float] = {}
        for _, row in constants_df.iterrows():
            cid = int(row[id_col])
            area = float(row[area_col])
            constant_lookup[cid] = area

        base_js_z = (
            self._old_control_current_density_base[mesh_key]
            if mesh_key in self._old_control_current_density_base
            else None
        )
        if base_js_z is None or base_js_z.shape[0] != js_np.shape[0]:
            base_js_z = np.asarray(js_np[:, 2], dtype=float).copy()
            self._old_control_current_density_base[mesh_key] = base_js_z

        for cell_id in selected_ids:
            if pfc2_mask[cell_id]:
                # JPF2=(xirad+IPF2)/APF2
                correction_current = radial_current
            elif pfc3_mask[cell_id]:
                # JPF3=(xiver+IPF3)/APF3
                correction_current = vertical_current
            elif pfc4_mask[cell_id]:
                # JPF4=(xiver+IPF4)/APF4
                correction_current = vertical_current
            elif pfc5_mask[cell_id]:
                # JPF5=(-xirad+IPF5)/APF5
                correction_current = -radial_current
            else:
                continue

            rc_id = int(real_constant_ids[cell_id])
            if rc_id not in constant_lookup:
                raise ValueError(
                    f"Real constant ID {rc_id} not found in constants table."
                )

            coil_area = constant_lookup[rc_id]
            if abs(coil_area) <= np.finfo(float).eps:
                raise ValueError(
                    f"Real constant ID {rc_id} has invalid coil area ({coil_area})."
                )

            correction_current_density = correction_current / coil_area
            js_np[cell_id, 2] = base_js_z[cell_id] + correction_current_density

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
        from .coordinate_system import CylindricalCoordinateSystem

        mesh = self._get_mesh(mesh_name)

        cell_data = mesh.GetCellData()

        plasma_array = cell_data.GetArray("PLASMA")
        if plasma_array is None:
            raise ValueError("CellData array 'PLASMA' not found on the mesh.")

        js_array = cell_data.GetArray("Js")
        if js_array is None:
            raise ValueError("CellData array 'Js' not found on the mesh.")
        js_np = np.asarray(vtknp.vtk_to_numpy(js_array), dtype=float)
        if js_np.ndim != 2 or js_np.shape[1] < 3:
            raise ValueError("CellData array 'Js' must have at least 3 components.")

        psie_array = cell_data.GetArray("PSIe")
        psi_elem = np.asarray(vtknp.vtk_to_numpy(psie_array), dtype=float).reshape(-1)

        be_array = cell_data.GetArray("Be")
        if be_array is None:
            raise ValueError("CellData array 'Be' not found on the mesh.")
        be_np = np.asarray(vtknp.vtk_to_numpy(be_array), dtype=float)
        if be_np.ndim != 2 or be_np.shape[1] < 2:
            raise ValueError("CellData array 'Be' must have at least 2 components.")
        bx = np.asarray(be_np[:, 0], dtype=float)
        by = np.asarray(be_np[:, 1], dtype=float)

        volume_array = cell_data.GetArray("cell_volume")
        cell_volume = np.asarray(vtknp.vtk_to_numpy(volume_array), dtype=float).reshape(-1)

        centroids_array = cell_data.GetArray("cell_centroids")
        centroids = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
        if centroids.ndim != 2 or centroids.shape[1] < 1:
            raise ValueError(
                "CellData array 'cell_centroids' must have at least 1 component."
            )
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
        if selected_ids.size == 0:
            self.status["total_plasma_current_before_scaling"] = 0.0
            self.status["plasma_current_scaling_lambda"] = 0.0
            self.status["plasma_element_count"] = 0
            return

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

        # cylindrical_coordinate_system = CylindricalCoordinateSystem()
        # ModelEdit.map_vector_field_to_global_components(
        #     mesh,
        #     array_name="Js",
        #     coordinate_system=cylindrical_coordinate_system,
        #     location="elem",
        #     selected_only=True,
        # )

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
        if psi_array is None:
            raise ValueError("PointData array 'PSI' not found on the mesh.")

        if point_data.GetArray("VACUUM_CHAMBER_NODES") is None:
            raise ValueError(
                "PointData array 'VACUUM_CHAMBER_NODES' not found on the mesh."
            )

        Selector.select(mesh, "VACUUM_CHAMBER_NODES", 1, type="node", mode="s")
        selection_array = point_data.GetArray("selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        if selected_ids.size == 0:
            raise ValueError(
                "No nodes selected with VACUUM_CHAMBER_NODES == 1; "
                "cannot check poloidal flux convergence."
            )

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
        # Placeholder for old evolution->equilibrium averaging.
        _ = equilibrium_mesh_name
        _ = evolution_mesh_name

    def toroidal_average_psi(
        self,
        equilibrium_mesh_name: Optional[str] = None,
        static_mesh_name: Optional[str] = None,
    ) -> None:
        # Placeholder for static->equilibrium PSI averaging.
        _ = equilibrium_mesh_name
        _ = static_mesh_name

    def toroidal_projection_plasma_current(
        self,
        static_mesh_name: Optional[str] = None,
        equilibrium_mesh_name: Optional[str] = None,
    ) -> None:
        # Placeholder for equilibrium->static current projection.
        _ = static_mesh_name
        _ = equilibrium_mesh_name

    def impressed_eddy_currents_coils_currents(
        self,
        evolution_mesh_name: Optional[str] = None,
        static_mesh_name: Optional[str] = None,
    ) -> None:
        # Placeholder for evolution->static current copy.
        _ = evolution_mesh_name
        _ = static_mesh_name

    def set_plasma_current(
        self,
        static_mesh_name: Optional[str] = None,
        evolution_mesh_name: Optional[str] = None,
    ) -> None:
        # Placeholder for static->evolution plasma current transfer.
        _ = static_mesh_name
        _ = evolution_mesh_name

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
            self.apec.write_output()
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

        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1

            if not self._first_absolute_iteration:
                self.find_separatrix(mesh_name=equilibrium_mesh)
                self.compute_jplasma_gradshafranov(mesh_name=equilibrium_mesh)

                self.toroidal_projection_plasma_current(
                    static_mesh_name=static_mesh,
                    equilibrium_mesh_name=equilibrium_mesh,
                )
                self.impressed_eddy_currents_coils_currents(
                    evolution_mesh_name=evolution_mesh,
                    static_mesh_name=static_mesh,
                )

                self.apec.set_active_mesh(static_mesh)
                self.apec.solve(nodal_solution=nodal_solution)
                self.calculate_poloidal_flux(mesh_name=static_mesh)

                self.toroidal_average_psi(
                    equilibrium_mesh_name=equilibrium_mesh,
                    static_mesh_name=static_mesh,
                )

                if self.check_convergence(mesh_name=equilibrium_mesh):
                    self.status["converged"] = True
                    return True

                self.set_plasma_current(
                    static_mesh_name=static_mesh,
                    evolution_mesh_name=evolution_mesh,
                )

            self.apec.set_active_mesh(evolution_mesh)
            self.apec.solve(nodal_solution=nodal_solution)
            self.calculate_poloidal_flux(mesh_name=evolution_mesh)
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

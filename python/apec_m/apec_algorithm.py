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
        self._b_target_axis_cache: Optional[
            Tuple[Tuple[Any, ...], np.ndarray, np.ndarray]
        ] = None
        self._step_overrides: Dict[str, Callable[..., Any]] = {}
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)
        self._supported_step_overrides = {
            self.STEP_FIND_SEPARATRIX,
            self.STEP_COMPUTE_JPLASMA,
            self.STEP_CONTROL_CURRENTS,
            self.STEP_CHECK_CONVERGENCE,
        }
        self._static_solver: Optional[Any] = None
        self._equilibrium_mesh: Optional[vtk.vtkUnstructuredGrid] = None

    def reset(self) -> None:
        self.status = copy.deepcopy(DEFAULT_APEC_ALGORITHM_STATUS)
        self._first_absolute_iteration = True
        self._previous_flux_vector = None
        self._old_control_current_density_base = None
        self._b_target_axis_cache = None
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)

    def set_toroidal_slices(self) -> None:
        if (
            self._toroidal_slices_nodes.size != 0
            and self._toroidal_slices_elems.size != 0
        ):
            return

        n_slices = int(self.input_parameters["n_slices"])
        elem_n_slices = n_slices

        if self._equilibrium_mesh is None:
            raise ValueError(
                "Equilibrium mesh is not set; cannot build toroidal slices."
            )

        if self.apec.mesh.GetMaxSpatialDimension() == 3:
            elem_n_slices = n_slices - 1
        if self.apec.mesh.GetMaxSpatialDimension() == 2 and (elem_n_slices != 1 or n_slices != 1):
            raise ValueError(
                "For 2D meshes, n_slices and elem_n_slices must both be 1."
            )
            
        from .coordinate_system import CylindricalCoordinateSystem
        from .cyclic import Cyclic

        tolerance = float(self.input_parameters["tolerance"])
        coordinate_system = CylindricalCoordinateSystem()
        self._toroidal_slices_nodes = Cyclic.build_toroidal_slices(
            self._equilibrium_mesh,
            self.apec.mesh,
            location="node",
            tolerance=tolerance,
            n_slices=n_slices,
            coordinate_system=coordinate_system,
        )
        self._toroidal_slices_elems = Cyclic.build_toroidal_slices(
            self._equilibrium_mesh,
            self.apec.mesh,
            location="elem",
            tolerance=tolerance,
            n_slices=elem_n_slices,
            coordinate_system=coordinate_system,
        )

    def _diagnostics_csv_path(self) -> str:
        diagnostics_dir = os.path.join(os.getcwd(), "__apec__")
        os.makedirs(diagnostics_dir, exist_ok=True)
        return os.path.join(diagnostics_dir, "apec_iteration_diagnostics.csv")

    def _clear_diagnostics_file(self) -> None:
        csv_path = self._diagnostics_csv_path()
        if os.path.exists(csv_path):
            os.remove(csv_path)

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
    ) -> Any:
        override = self._step_overrides[step] if step in self._step_overrides else None
        if override is None:
            return default_impl()
        return override(self)

    def set_psi_magnetic_axis(self) -> None:
        mesh = self._equilibrium_mesh

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

    def reset_plasma_components(self) -> None:
        mesh = self._equilibrium_mesh
        cell_data = mesh.GetCellData()

        js_array = cell_data.GetArray("Js")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.set_value_on_selection(mesh, "PLASMA", 1, type="elem")

        selection_array = cell_data.GetArray("selection_cell")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        js_np = vtknp.vtk_to_numpy(js_array)

        js_np[selected_ids, :] = 0.0
        js_array.Modified()

    def identify_plasma_region(self) -> None:
        mesh = self._equilibrium_mesh

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

    def remove_plasma_out_of_X_point(self) -> None:
        mesh = self._equilibrium_mesh
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

    def check_limiter(self) -> None:
        mesh = self._equilibrium_mesh

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

    def check_X_point(self) -> None:
        mesh = self._equilibrium_mesh
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

    def set_b_target_axis(self) -> None:
        mesh = self._equilibrium_mesh

        point_data = mesh.GetPointData()
        b_array = point_data.GetArray("B")
        b = np.asarray(vtknp.vtk_to_numpy(b_array), dtype=float)

        target_point = (
            float(self.input_parameters["target_central_axis_x"]),
            float(self.input_parameters["target_central_axis_y"]),
            float(self.input_parameters["target_central_axis_z"]),
        )

        cache_key = (
            id(mesh),
            mesh.GetNumberOfPoints(),
            mesh.GetNumberOfCells(),
            target_point,
        )
        if (
            self._b_target_axis_cache is None
            or self._b_target_axis_cache[0] != cache_key
        ):
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
            point_ids_np = np.fromiter(
                (int(point_ids.GetId(local_idx)) for local_idx in range(n_points)),
                dtype=np.int64,
                count=n_points,
            )
            weights_np = np.asarray(weights, dtype=float)
            self._b_target_axis_cache = (cache_key, point_ids_np, weights_np)
        else:
            _, point_ids_np, weights_np = self._b_target_axis_cache

        b_target = weights_np @ b[point_ids_np, :2]
        target = (float(b_target[0]), float(b_target[1]))
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

    def apply_control_currents(self) -> None:
        mesh = self._equilibrium_mesh

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

    def find_separatrix(self) -> None:
        self._dispatch_step(
            self.STEP_FIND_SEPARATRIX,
            self._find_separatrix_default,
        )

    def _find_separatrix_default(self) -> None:
        self.reset_plasma_components()
        self.set_psi_magnetic_axis()
        self.remove_plasma_out_of_X_point()
        self.check_limiter()
        self.check_X_point()
        psi_limiter = float(self.status["last_psi_limiter"])
        psi_x_point = float(self.status["last_psi_x_point"])
        if psi_limiter > psi_x_point:
            self.status["psi_min"] = psi_limiter
            self.status["separatrix_configuration"] = "limiter"
        else:
            self.status["psi_min"] = psi_x_point
            self.status["separatrix_configuration"] = "X-point"
        self.identify_plasma_region()
    def compute_jplasma_gradshafranov(self) -> None:
        self._dispatch_step(
            self.STEP_COMPUTE_JPLASMA,
            self._compute_jplasma_gradshafranov_default,
        )

    def _compute_jplasma_gradshafranov_default(self) -> None:
        mesh = self._equilibrium_mesh

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
        psi_segn = np.maximum(psi_segn, 0.0)
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

    def control_currents(self) -> None:
        self._dispatch_step(
            self.STEP_CONTROL_CURRENTS,
            self._control_currents_default,
        )

    def _control_currents_default(self) -> None:
        self.set_b_target_axis()
        self.compute_correction_currents()
        self.apply_control_currents()

    def check_convergence(self) -> bool:
        result = self._dispatch_step(
            self.STEP_CHECK_CONVERGENCE,
            self._check_convergence_default,
        )
        return bool(result)

    def _check_convergence_default(self) -> bool:
        mesh = self._equilibrium_mesh

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
    def check_poloidal_flux_convergence(self) -> bool:
        return self.check_convergence()

    def toroidal_average(self, psi_only: bool = False) -> None:
        equilibrium_mesh = self._equilibrium_mesh
        evolution_mesh = self.apec.mesh
        self.set_toroidal_slices()
        node_slices = self._toroidal_slices_nodes
        elem_slices = self._toroidal_slices_elems

        eq_point_data = equilibrium_mesh.GetPointData()
        evol_point_data = evolution_mesh.GetPointData()

        eq_psi_array = eq_point_data.GetArray("PSI")
        evol_psi_array = evol_point_data.GetArray("PSI")
        if eq_psi_array is None or evol_psi_array is None:
            raise ValueError(
                "Required arrays for toroidal_average are missing (PSI)."
            )

        eq_psi = np.asarray(vtknp.vtk_to_numpy(eq_psi_array), dtype=float).reshape(-1)
        evol_psi = np.asarray(vtknp.vtk_to_numpy(evol_psi_array), dtype=float).reshape(-1)

        src_nodes = np.asarray(node_slices[:, 0], dtype=np.int64)
        dst_nodes = np.asarray(node_slices[:, 1:], dtype=np.int64)

        psi_values = np.mean(evol_psi[dst_nodes], axis=1)
        eq_psi[src_nodes] = psi_values
        eq_psi_array.Modified()

        if psi_only:
            return

        eq_cell_data = equilibrium_mesh.GetCellData()
        evol_cell_data = evolution_mesh.GetCellData()

        eq_b_array = eq_point_data.GetArray("B")
        evol_b_array = evol_point_data.GetArray("B")
        eq_be_array = eq_cell_data.GetArray("Be")
        evol_be_array = evol_cell_data.GetArray("Be")
        eq_psie_array = eq_cell_data.GetArray("PSIe")
        evol_psie_array = evol_cell_data.GetArray("PSIe")
        if (
            eq_b_array is None
            or evol_b_array is None
            or eq_be_array is None
            or evol_be_array is None
            or eq_psie_array is None
            or evol_psie_array is None
        ):
            raise ValueError(
                "Required arrays for toroidal_average are missing "
                "(PSI, B, Be, PSIe)."
            )

        eq_b = np.asarray(vtknp.vtk_to_numpy(eq_b_array), dtype=float)
        eq_be = np.asarray(vtknp.vtk_to_numpy(eq_be_array), dtype=float)
        evol_b = np.asarray(vtknp.vtk_to_numpy(evol_b_array), dtype=float)
        evol_be = np.asarray(vtknp.vtk_to_numpy(evol_be_array), dtype=float)
        eq_psie = np.asarray(vtknp.vtk_to_numpy(eq_psie_array), dtype=float).reshape(-1)
        evol_psie = np.asarray(vtknp.vtk_to_numpy(evol_psie_array), dtype=float).reshape(-1)

        from .coordinate_system import CylindricalCoordinateSystem

        cylindrical = CylindricalCoordinateSystem()

        src_elems = np.asarray(elem_slices[:, 0], dtype=np.int64)
        dst_elems = np.asarray(elem_slices[:, 1:], dtype=np.int64)

        converted_vectors_to_cylindrical = evolution_mesh.GetMaxSpatialDimension() > 2
        if converted_vectors_to_cylindrical:
            flat_dst_nodes = dst_nodes.reshape(-1)
            flat_dst_elems = dst_elems.reshape(-1)
            evol_node_positions = ModelEdit._get_positions(evolution_mesh, "node")
            evol_elem_positions = ModelEdit._get_positions(evolution_mesh, "elem")
            evol_b_values = cylindrical.vector_components_to_local(
                evol_b[flat_dst_nodes],
                evol_node_positions[flat_dst_nodes],
            ).reshape((*dst_nodes.shape, 3))
            evol_be_values = cylindrical.vector_components_to_local(
                evol_be[flat_dst_elems],
                evol_elem_positions[flat_dst_elems],
            ).reshape((*dst_elems.shape, 3))
        else:
            evol_b_values = evol_b[dst_nodes]
            evol_be_values = evol_be[dst_elems]

        eq_b[src_nodes] = np.mean(evol_b_values, axis=1)
        eq_be[src_elems] = np.mean(evol_be_values, axis=1)

        if converted_vectors_to_cylindrical:
            eq_node_positions = ModelEdit._get_positions(equilibrium_mesh, "node")
            eq_elem_positions = ModelEdit._get_positions(equilibrium_mesh, "elem")
            eq_b[src_nodes] = cylindrical.vector_components_to_global(
                eq_b[src_nodes],
                eq_node_positions[src_nodes],
            )
            eq_be[src_elems] = cylindrical.vector_components_to_global(
                eq_be[src_elems],
                eq_elem_positions[src_elems],
            )

        eq_b_array.Modified()
        eq_be_array.Modified()

        psie_values = np.mean(evol_psie[dst_elems], axis=1)
        eq_psie[src_elems] = psie_values
        eq_psie_array.Modified()

    def toroidal_projection_current(self) -> None:
        from .coordinate_system import CylindricalCoordinateSystem

        static_mesh = self.apec.mesh
        equilibrium_mesh = self._equilibrium_mesh
        self.set_toroidal_slices()
        elem_slices = self._toroidal_slices_elems

        static_js_array = static_mesh.GetCellData().GetArray("Js")
        equilibrium_js_array = equilibrium_mesh.GetCellData().GetArray("Js")

        cylindrical = CylindricalCoordinateSystem()
        static_js = np.asarray(vtknp.vtk_to_numpy(static_js_array), dtype=float)
        equilibrium_js = np.asarray(vtknp.vtk_to_numpy(equilibrium_js_array), dtype=float)

        src_elems = np.asarray(elem_slices[:, 0], dtype=np.int64)
        dst_elems = np.asarray(elem_slices[:, 1:], dtype=np.int64)
        projected_js = equilibrium_js[src_elems]
        static_mesh_is_3d = static_mesh.GetMaxSpatialDimension() > 2
        if static_mesh_is_3d:
            projected_js = projected_js[:, [0, 2, 1]]

        if static_mesh_is_3d:
            flat_dst_elems = dst_elems.reshape(-1)
            local_js = np.broadcast_to(
                projected_js[:, np.newaxis, :],
                (*dst_elems.shape, 3),
            ).reshape(-1, 3)
            elem_positions = ModelEdit._get_positions(static_mesh, "elem")
            static_js[flat_dst_elems] = cylindrical.vector_components_to_global(
                local_js,
                elem_positions[flat_dst_elems],
            )
            selection_mask = np.zeros(static_mesh.GetNumberOfCells(), dtype=bool)
            selection_mask[flat_dst_elems] = True
            Selector.apply_boolean_selection_mask(
                selection_mask,
                static_mesh.GetCellData(),
                "elem",
                "s",
            )
        else:
            static_js[dst_elems] = projected_js[:, np.newaxis, :]

        static_js_array.Modified()

    def append_iteration_diagnostics(
        self,
        converged: bool,
    ) -> None:
        self.status["converged"] = bool(converged)
        if self.apec._world_rank() != 0:
            return

        row = self._build_diagnostics_row(converged)
        row_df = pd.DataFrame([row])
        csv_path = self._diagnostics_csv_path()
        write_header = (not os.path.exists(csv_path)) or os.path.getsize(csv_path) == 0
        row_df.to_csv(csv_path, mode="a", header=write_header, index=False)

    def calculate_poloidal_flux(self) -> None:
        from .coordinate_system import CylindricalCoordinateSystem
        from .cyclic import Cyclic

        mesh = self.apec.mesh
        ModelEdit.element_solution(mesh, "A")
        Cyclic.poloidal_flux(
            mesh=mesh,
            coordinate_system=CylindricalCoordinateSystem(),
        )

    def equilibrium_solve(self) -> bool:
        if self._first_absolute_iteration:
            self._clear_diagnostics_file()
        self._equilibrium_mesh = self.apec.mesh
        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )
        self.status["converged"] = False

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1

            if not self._first_absolute_iteration:
                self.find_separatrix()
                self.compute_jplasma_gradshafranov()
                self.control_currents()

            self.apec.solve()
            self.calculate_poloidal_flux()

            converged = self.check_convergence()
            # self.apec.write_output()
            self.append_iteration_diagnostics(converged)
            if converged:
                self.status["converged"] = True
                return True
            self._first_absolute_iteration = False

        self.status["converged"] = False
        return False

    def evolution_solve(self) -> bool:
        if self._first_absolute_iteration:
            self._clear_diagnostics_file()
            Selector.select(self._equilibrium_mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
            Selector.nsle(self._equilibrium_mesh, mode="s")
            Selector.select_all(self.apec.mesh)
            self.set_toroidal_slices()

        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1

            if not self._first_absolute_iteration:
                self.toroidal_average()
                self.find_separatrix()
                self.compute_jplasma_gradshafranov()

                self.toroidal_projection_current()
                self.apec.solve()
                self.calculate_poloidal_flux()
                # self.apec.write_output()

                self.toroidal_average(psi_only=True)

                converged = self.check_convergence()
                self.append_iteration_diagnostics(
                    converged,
                )
                if converged:
                    self.status["converged"] = True
                    return True
            else:
                self.apec.solve()
                self.calculate_poloidal_flux()
                # self.apec.write_output()
                self.check_convergence()
                self._first_absolute_iteration = False
                self.status["converged"] = True
                return True

        self.status["converged"] = False
        return False

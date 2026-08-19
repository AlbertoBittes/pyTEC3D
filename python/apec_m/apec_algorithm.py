from __future__ import annotations

import copy
import os
from typing import Any, Callable, Dict, List, Optional, TYPE_CHECKING, Tuple, cast

import numpy as np
import pandas as pd
import vtk
try:
    from vmtk.vtkvmtkDifferentialGeometryPython import vtkvmtkUnstructuredGridGradientFilter
except ImportError:
    vtkvmtkUnstructuredGridGradientFilter = None

from .selector import Selector
from .model_edit import ModelEdit
from .coordinate_system import CylindricalCoordinateSystem, YAxisCylindricalCoordinateSystem
from .logging_utils import get_logger
from .validation_utils import ensure_vtk_unstructured_grid
from .constants import (
    DEFAULT_APEC_ALGORITHM_INPUT_PARAMETERS,
    DEFAULT_APEC_ALGORITHM_STATUS,
)

import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    from .apec_m import APEC_M


LOGGER = get_logger("apec_m.algorithm")

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
        """Initialize algorithm state and optional parameter overrides.

        Args:
            apec: Owning APEC_M object.
            input_parameters: Optional input-parameter overrides.

        Returns:
            None
        """
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
        self._halo_node_groups: Tuple[np.ndarray, ...] = ()
        self._halo_boundary_graph_cache: Optional[
            Tuple[np.ndarray, Dict[int, np.ndarray]]
        ] = None
        self._supported_step_overrides = {
            self.STEP_FIND_SEPARATRIX,
            self.STEP_COMPUTE_JPLASMA,
            self.STEP_CONTROL_CURRENTS,
            self.STEP_CHECK_CONVERGENCE,
        }
        self._static_solver: Optional[Any] = None
        self.equilibrium_mesh: Optional[vtk.vtkUnstructuredGrid] = None
        LOGGER.debug(
            "ApecAlgorithm initialized%s.",
            " with parameter overrides" if input_parameters else "",
        )

    def reset(self) -> None:
        """Reset iterative state and cached data to default values.

        Returns:
            None
        """
        self.status = copy.deepcopy(DEFAULT_APEC_ALGORITHM_STATUS)
        self._first_absolute_iteration = True
        self._previous_flux_vector = None
        self._old_control_current_density_base = None
        self._b_target_axis_cache = None
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)
        self._halo_node_groups = ()
        self._halo_boundary_graph_cache = None
        LOGGER.debug("ApecAlgorithm state reset.")

    def _require_equilibrium_mesh(self) -> vtk.vtkUnstructuredGrid:
        """Return a validated equilibrium mesh required by internal algorithm steps.

        Returns:
            vtk.vtkUnstructuredGrid: Validated equilibrium mesh.
        """
        if self.equilibrium_mesh is None:
            raise ValueError("equilibrium_mesh is not set.")
        return ensure_vtk_unstructured_grid(
            self.equilibrium_mesh,
            arg_name="equilibrium_mesh",
        )

    @staticmethod
    def _require_array(
        data: "vtk.vtkDataSetAttributes",
        name: str,
    ) -> "vtk.vtkDataArray":
        """Return a named VTK array or raise a descriptive error if missing.

        Args:
            data: VTK data attributes object (CellData or PointData).
            name: Array name to look up.

        Returns:
            vtk.vtkDataArray: Found array.

        Raises:
            ValueError: If the named array is not present in the data object.
        """
        arr = data.GetArray(name)
        if arr is None:
            raise ValueError(
                f"Required array '{name}' not found. "
                "Ensure the mesh is fully initialized before running the algorithm."
            )
        return arr

    @staticmethod
    def build_toroidal_slices(
        input_mesh: "vtk.vtkUnstructuredGrid",
        target_mesh: "vtk.vtkUnstructuredGrid",
        location: str = "node",
        tolerance: float = 1e-6,
        extrusion_angle: float = 360.0,
        n_divisions: Optional[int] = None,
        coordinate_system: Optional["CylindricalCoordinateSystem"] = None,
        n_slices: Optional[int] = None,
    ) -> np.ndarray:
        """
        Build a toroidal projection map from a 2D input mesh to a target mesh.

        The first column stores source IDs; remaining columns store matched target IDs.
        For 3D targets this returns one column per toroidal station.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            target_mesh: Target unstructured grid (typically 3D).
            location: "node" to map point IDs or "elem" to map cell IDs.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            extrusion_angle: Total toroidal extrusion angle in degrees.
            n_divisions: Number of toroidal element divisions.
            coordinate_system: Cylindrical coordinate system used for local coordinates.
            n_slices: Backward-compatible alias for n_divisions.

        Returns:
            np.ndarray: Integer map with shape (num_input, 1 + num_matches).
        """
        ensure_vtk_unstructured_grid(input_mesh, arg_name="input_mesh")
        ensure_vtk_unstructured_grid(target_mesh, arg_name="target_mesh")
        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")
        if tolerance <= 0:
            raise ValueError("tolerance must be positive.")
        if extrusion_angle <= 0:
            raise ValueError("extrusion_angle must be positive.")

        if input_mesh.GetNumberOfCells() > 0 and input_mesh.GetMaxSpatialDimension() > 2:
            raise ValueError("Input mesh must be 2D.")

        if n_divisions is None and n_slices is not None:
            n_divisions = n_slices

        if coordinate_system is None:
            coordinate_system = YAxisCylindricalCoordinateSystem()
        assert coordinate_system is not None
        if coordinate_system.type != "cylindrical":
            raise ValueError(
                "coordinate_system must be cylindrical for toroidal projection."
            )

        target_is_3d = target_mesh.GetMaxSpatialDimension() > 2
        if target_is_3d:
            if n_divisions is None:
                raise ValueError("n_divisions must be provided for 3D target meshes.")
            try:
                n_divisions = int(n_divisions)
            except (TypeError, ValueError) as exc:
                raise ValueError("n_divisions must be a positive integer.") from exc
            if n_divisions <= 0:
                raise ValueError("n_divisions must be a positive integer.")
            num_matches = n_divisions + 1 if location == "node" else n_divisions
            if location == "node":
                expected_theta = np.linspace(
                    0.0,
                    float(extrusion_angle),
                    n_divisions + 1,
                    dtype=float,
                )
            else:
                step = float(extrusion_angle) / float(n_divisions)
                expected_theta = (np.arange(n_divisions, dtype=float) + 0.5) * step
        else:
            if input_mesh.GetMaxSpatialDimension() > 2:
                raise ValueError("2D target meshes can only be mapped from 2D input meshes.")
            num_matches = 1
            expected_theta = np.zeros((1,), dtype=float)

        def angular_distance(theta: np.ndarray, expected: float) -> np.ndarray:
            return np.abs((theta - expected + 180.0) % 360.0 - 180.0)

        input_local, input_ids = ModelEdit.to_local_coordinates(
            input_mesh,
            coordinate_system,
            location=location,
            selected_only=True,
        )
        target_local, target_ids = ModelEdit.to_local_coordinates(
            target_mesh,
            coordinate_system,
            location=location,
            selected_only=True,
        )

        if input_local.size == 0:
            return np.empty((0, 1 + num_matches), dtype=np.int64)

        input_radius = input_local[:, 0]
        input_axial = input_local[:, 2]
        target_radius = target_local[:, 0]
        target_theta = target_local[:, 1]
        target_axial = target_local[:, 2]

        target_count = target_radius.shape[0]
        if target_count == 0:
            raise ValueError("Target mesh has no points/cells to map.")
        if num_matches > target_count:
            raise ValueError(
                "Number of requested toroidal matches cannot exceed the number of "
                "available target points/cells."
            )

        num_input = input_radius.shape[0]
        tolerance_current = tolerance
        counts = np.empty(num_input, dtype=int)

        axis_mask = input_radius <= tolerance
        required_counts = np.where(axis_mask, 1, num_matches)

        max_iterations = 50
        for _ in range(max_iterations):
            for input_id in range(num_input):
                dr = np.abs(target_radius - input_radius[input_id])
                dz = np.abs(target_axial - input_axial[input_id])
                mask = (dr <= tolerance_current) & (dz <= tolerance_current)
                counts[input_id] = int(np.count_nonzero(mask))
            if np.all(counts >= required_counts):
                break
            tolerance_current *= 2.0
        else:
            raise ValueError(
                "Toroidal projection map did not converge within the "
                "tolerance expansion limit."
            )

        mapping = np.full((num_input, 1 + num_matches), -1, dtype=np.int64)
        mapping[:, 0] = input_ids.astype(np.int64, copy=False)
        for row_idx in range(num_input):
            dr = np.abs(target_radius - input_radius[row_idx])
            dz = np.abs(target_axial - input_axial[row_idx])
            mask = (dr <= tolerance_current) & (dz <= tolerance_current)
            candidate_idx = np.where(mask)[0]
            if candidate_idx.size == 0:
                continue
            candidate_ids = target_ids[candidate_idx]
            if axis_mask[row_idx]:
                trimmed = np.full(
                    (num_matches,),
                    candidate_ids[0],
                    dtype=candidate_ids.dtype,
                )
            else:
                trimmed = np.empty((num_matches,), dtype=candidate_ids.dtype)
                used = np.zeros(candidate_idx.size, dtype=bool)
                poloidal_dist = np.hypot(dr[candidate_idx], dz[candidate_idx])
                for theta_idx, theta_value in enumerate(expected_theta):
                    theta_dist = angular_distance(
                        target_theta[candidate_idx],
                        theta_value,
                    )
                    score = theta_dist + poloidal_dist / max(
                        tolerance_current,
                        np.finfo(float).eps,
                    )
                    if np.any(~used):
                        score = np.where(used, np.inf, score)
                    pick = int(np.argmin(score))
                    trimmed[theta_idx] = candidate_ids[pick]
                    used[pick] = True
            mapping[row_idx, 1 : 1 + trimmed.size] = trimmed[:num_matches].astype(
                np.int64,
                copy=False,
            )

        return mapping

    def _set_toroidal_slices(self) -> None:
        """Build toroidal node/element mapping slices when they are not cached yet.

        Returns:
            None
        """
        if (
            self._toroidal_slices_nodes.size != 0
            and self._toroidal_slices_elems.size != 0
        ):
            LOGGER.debug("Toroidal slices already cached; skipping rebuild.")
            return
        LOGGER.debug("Building toroidal node/element slices.")

        equilibrium_mesh = self._require_equilibrium_mesh()

        target_dimension = self.apec.mesh.GetMaxSpatialDimension()
        extrusion_angle = float(self.input_parameters.get("extrusion_angle", 360.0))
        n_toroidal_divisions = int(
            self.input_parameters.get("n_toroidal_divisions", -1)
        )
        if target_dimension > 2 and n_toroidal_divisions <= 0:
            raise ValueError(
                "n_toroidal_divisions must be a positive integer for 3D meshes."
            )
        if target_dimension <= 2 and n_toroidal_divisions <= 0:
            n_toroidal_divisions = 1

        tolerance = float(self.input_parameters["tolerance"])
        coordinate_system = YAxisCylindricalCoordinateSystem()
        self._toroidal_slices_nodes = self.build_toroidal_slices(
            equilibrium_mesh,
            self.apec.mesh,
            location="node",
            tolerance=tolerance,
            extrusion_angle=extrusion_angle,
            n_divisions=n_toroidal_divisions,
            coordinate_system=coordinate_system,
        )
        self._toroidal_slices_elems = self.build_toroidal_slices(
            equilibrium_mesh,
            self.apec.mesh,
            location="elem",
            tolerance=tolerance,
            extrusion_angle=extrusion_angle,
            n_divisions=n_toroidal_divisions,
            coordinate_system=coordinate_system,
        )

    def _diagnostics_csv_path(self) -> str:
        """Return diagnostics CSV output path under the __apec__ folder.

        Returns:
            str: Diagnostics CSV file path.
        """
        diagnostics_dir = os.path.join(os.getcwd(), "__apec__")
        os.makedirs(diagnostics_dir, exist_ok=True)
        return os.path.join(diagnostics_dir, "apec_iteration_diagnostics.csv")

    def _clear_diagnostics_file(self) -> None:
        """Remove diagnostics CSV on rank 0 before a new run.

        Returns:
            None
        """
        if self.apec._world_rank() != 0:
            return

        csv_path = self._diagnostics_csv_path()
        try:
            os.remove(csv_path)
        except FileNotFoundError:
            pass

    def _build_diagnostics_row(self, converged: bool) -> Dict[str, Any]:
        """Build one diagnostics record from current algorithm status.

        Args:
            converged: Current convergence flag.

        Returns:
            Dict[str, Any]: Row payload for CSV append.
        """
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
        """Register or clear a custom implementation for a supported algorithm step.

        Args:
            step: Step identifier.
            override: Callable override or None to clear.

        Returns:
            None
        """
        if step not in self._supported_step_overrides:
            raise ValueError(
                f"Unsupported step '{step}'. "
                f"Supported steps: {sorted(self._supported_step_overrides)}."
            )
        if override is None:
            self._step_overrides.pop(step, None)
            LOGGER.debug("Step override cleared for '%s'.", step)
            return
        if not callable(override):
            raise TypeError("override must be callable or None.")
        self._step_overrides[step] = override
        LOGGER.debug("Step override registered for '%s': %s.", step, override)

    def clear_step_overrides(self) -> None:
        """Remove all previously registered step overrides.

        Returns:
            None
        """
        self._step_overrides.clear()

    def _normalized_control_coils(self) -> List[Dict[str, Any]]:
        """Validate and normalize control-coil mapping from input parameters.

        Returns:
            List[Dict[str, Any]]: Normalized list with component/radial_gain/vertical_gain.
        """
        raw = self.input_parameters.get("control_coils", [])
        if raw is None:
            return []
        if not isinstance(raw, list):
            raise TypeError(
                "input_parameters['control_coils'] must be a list of dictionaries."
            )

        normalized: List[Dict[str, Any]] = []
        for idx, item in enumerate(raw, start=1):
            if not isinstance(item, dict):
                raise TypeError(
                    f"control_coils[{idx}] must be a dictionary, got {type(item).__name__}."
                )

            component = str(item.get("component", "")).strip()
            if not component:
                raise ValueError(
                    f"control_coils[{idx}] must define a non-empty 'component'."
                )

            if "radial_gain" not in item:
                raise ValueError(
                    f"control_coils[{idx}] must define 'radial_gain'."
                )
            if "vertical_gain" not in item:
                raise ValueError(
                    f"control_coils[{idx}] must define 'vertical_gain'."
                )

            radial_gain = item["radial_gain"]
            vertical_gain = item["vertical_gain"]
            normalized.append(
                {
                    "component": component,
                    "radial_gain": float(radial_gain),
                    "vertical_gain": float(vertical_gain),
                }
            )
        return normalized

    def _cache_halo_boundary_graph(self) -> None:
        """Cache the external boundary graph of selected vacuum-chamber cells.

        Returns:
            None
        """
        if self._halo_boundary_graph_cache is not None:
            LOGGER.debug("Halo boundary graph already cached; skipping rebuild.")
            return
        LOGGER.debug("Building halo boundary graph.")
        mesh = self._require_equilibrium_mesh()
        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        boundary_node_ids, boundary_adjacency = (
            Selector.get_selected_external_node_adjacency(mesh)
        )
        self._halo_boundary_graph_cache = (boundary_node_ids, boundary_adjacency)

    def _get_halo_boundary_graph(self) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
        """Return cached halo boundary graph, computing it on first use.

        Returns:
            Tuple[np.ndarray, Dict[int, np.ndarray]]: Boundary node ids and adjacency map.
        """
        if self._halo_boundary_graph_cache is None:
            self._cache_halo_boundary_graph()
        if self._halo_boundary_graph_cache is None:
            raise RuntimeError("Halo boundary graph cache is unexpectedly empty.")
        return self._halo_boundary_graph_cache

    @staticmethod
    def _connected_node_groups(
        node_ids: np.ndarray,
        boundary_adjacency: Dict[int, np.ndarray],
    ) -> Tuple[np.ndarray, ...]:
        """Split node ids into ordered edge-connected components on the boundary graph.

        Args:
            node_ids: Candidate node ids.
            boundary_adjacency: Boundary adjacency graph.

        Returns:
            Tuple[np.ndarray, ...]: Ordered connected node groups.
        """
        node_ids = np.unique(np.asarray(node_ids, dtype=np.int64).reshape(-1))
        if node_ids.size == 0:
            return ()
        if node_ids.size == 1:
            return (node_ids,)

        selected = set(int(node_id) for node_id in node_ids)
        visited: set[int] = set()
        groups: list[np.ndarray] = []

        for start_node in node_ids:
            start = int(start_node)
            if start in visited:
                continue

            stack = [start]
            component: list[int] = []
            visited.add(start)
            while stack:
                node = stack.pop()
                component.append(node)
                for neighbor in boundary_adjacency.get(node, ()):
                    neighbor = int(neighbor)
                    if neighbor in selected and neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)

            component_set = set(component)
            if len(component) <= 2:
                groups.append(np.asarray(component, dtype=np.int64))
                continue

            degrees = {
                node: sum(
                    int(neighbor) in component_set
                    for neighbor in boundary_adjacency.get(node, ())
                )
                for node in component
            }
            if any(degree > 2 for degree in degrees.values()):
                raise ValueError(
                    "Halo node group has branching edge connectivity; "
                    "cannot order it as a single boundary path."
                )

            endpoints = [node for node, degree in degrees.items() if degree == 1]
            if len(endpoints) != 2:
                raise ValueError("Halo node group is not an open edge-connected path.")

            ordered: list[int] = []
            previous = -1
            current = min(endpoints)
            ordered_set: set[int] = set()
            while current not in ordered_set:
                ordered.append(current)
                ordered_set.add(current)
                next_candidates = [
                    int(neighbor)
                    for neighbor in boundary_adjacency.get(current, ())
                    if (
                        int(neighbor) in component_set
                        and int(neighbor) != previous
                        and int(neighbor) not in ordered_set
                    )
                ]
                if not next_candidates:
                    break
                previous, current = current, next_candidates[0]

            if len(ordered) != len(component):
                raise ValueError("Halo node group is not a simple edge-connected path.")
            groups.append(np.asarray(ordered, dtype=np.int64))

        return tuple(sorted(groups, key=lambda group: int(group[0])))

    @staticmethod
    def _extract_main_plasma_region_cells(
        mesh: "vtk.vtkUnstructuredGrid",
        axis_node_id: int,
        psi_point: np.ndarray,
        psi_x: float,
        vac_cell_mask: np.ndarray,
    ) -> np.ndarray:
        """Traverses VACUUM_CHAMBER elements starting from the magnetic axis.
        
        Candidate elements must have ALL nodes satisfying PSI >= PSI_x.
        Adjacent elements MUST share >= 2 nodes to prevent BFS traversal 
        across a single-node X-point saddle contact.
        """
        n_cells = mesh.GetNumberOfCells()
        n_points = mesh.GetNumberOfPoints()

        if axis_node_id < 0 or axis_node_id >= n_points:
            return np.zeros(0, dtype=np.int64)

        psi_axis = psi_point[axis_node_id]
        if not np.isfinite(psi_axis):
            return np.zeros(0, dtype=np.int64)

        core_above = psi_axis >= psi_x

        cells = mesh.GetCells()
        offsets_array = cells.GetOffsetsArray()
        connectivity_array = cells.GetConnectivityArray()

        if offsets_array is None or connectivity_array is None:
            return np.zeros(0, dtype=np.int64)

        offsets = np.asarray(vtknp.vtk_to_numpy(offsets_array), dtype=np.int64)
        connectivity = np.asarray(vtknp.vtk_to_numpy(connectivity_array), dtype=np.int64)

        # 1. Filter VACUUM_CHAMBER cells where ALL nodes satisfy PSI condition
        candidate_cell_mask = np.zeros(n_cells, dtype=bool)
        vac_cell_ids = np.flatnonzero(vac_cell_mask)

        for cell_id in vac_cell_ids:
            start = offsets[cell_id]
            end = offsets[cell_id + 1]
            cell_nodes = connectivity[start:end]
            cell_psi = psi_point[cell_nodes]

            if not np.all(np.isfinite(cell_psi)):
                continue

            all_nodes_valid = np.all(cell_psi >= psi_x) if core_above else np.all(cell_psi <= psi_x)
            if all_nodes_valid:
                candidate_cell_mask[cell_id] = True

        # 2. Build reverse lookup: Node ID -> List of Candidate Cell IDs
        node_to_cand_cells = [[] for _ in range(n_points)]
        cand_cell_ids = np.flatnonzero(candidate_cell_mask)

        for cell_id in cand_cell_ids:
            start = offsets[cell_id]
            end = offsets[cell_id + 1]
            for node_id in connectivity[start:end]:
                node_to_cand_cells[node_id].append(cell_id)

        # 3. Seed BFS from candidate cells touching the magnetic axis node
        seed_cells = node_to_cand_cells[axis_node_id]
        if not seed_cells:
            return np.zeros(0, dtype=np.int64)

        # 4. Run cell-level BFS with strict edge-connectivity check (>= 2 common nodes)
        from collections import deque

        visited_cells = np.zeros(n_cells, dtype=bool)
        queue = deque(seed_cells)
        for cid in seed_cells:
            visited_cells[cid] = True

        while queue:
            curr_cid = queue.popleft()
            start = offsets[curr_cid]
            end = offsets[curr_cid + 1]
            curr_nodes = set(connectivity[start:end])

            # Get potential neighbor candidate cells sharing at least 1 node
            potential_neighbors = set()
            for node_id in curr_nodes:
                potential_neighbors.update(node_to_cand_cells[node_id])

            # Test edge-connectivity (>= 2 common nodes required)
            for nbr_cid in potential_neighbors:
                if visited_cells[nbr_cid]:
                    continue

                nbr_start = offsets[nbr_cid]
                nbr_end = offsets[nbr_cid + 1]
                nbr_nodes = connectivity[nbr_start:nbr_end]

                # Count overlapping nodes between current cell and neighbor
                common_nodes_count = sum(1 for node in nbr_nodes if node in curr_nodes)

                # If common_nodes_count == 1, it is a point-contact (X-point node) -> REJECT
                if common_nodes_count >= 2:
                    visited_cells[nbr_cid] = True
                    queue.append(nbr_cid)

        return np.flatnonzero(visited_cells).astype(np.int64)

    @staticmethod
    def _orient_halo_node_groups(
        groups: Tuple[np.ndarray, ...],
        points: np.ndarray,
        reference_point: Tuple[float, float, float],
        plasma_current: float,
    ) -> Tuple[np.ndarray, ...]:
        """Orient each halo node path consistently with the plasma-current direction.

        Args:
            groups: Ordered halo node groups.
            points: Mesh point coordinates.
            reference_point: Magnetic axis reference point.
            plasma_current: Plasma current sign/magnitude.

        Returns:
            Tuple[np.ndarray, ...]: Oriented node groups.
        """
        if not groups or plasma_current == 0.0:
            return groups

        reference_xy = np.asarray(reference_point[:2], dtype=float)
        target_clockwise = plasma_current > 0.0
        oriented_groups: list[np.ndarray] = []

        for group in groups:
            if group.size <= 1:
                oriented_groups.append(group)
                continue

            xy = np.asarray(points[group, :2], dtype=float) - reference_xy
            cross = xy[:-1, 0] * xy[1:, 1] - xy[:-1, 1] * xy[1:, 0]
            dot = np.einsum("ij,ij->i", xy[:-1], xy[1:])
            signed_angle = float(np.sum(np.arctan2(cross, dot)))
            if abs(signed_angle) <= np.finfo(float).eps:
                oriented_groups.append(group)
                continue

            is_clockwise = signed_angle < 0.0
            if is_clockwise != target_clockwise:
                group = group[::-1].copy()
            oriented_groups.append(group)

        return tuple(oriented_groups)

    def _dispatch_step(
        self,
        step: str,
        default_impl: Callable[..., Any],
    ) -> Any:
        """Dispatch to an override when provided, else execute the default implementation.

        Args:
            step: Step identifier.
            default_impl: Default callable implementation.

        Returns:
            Any: Step return value.
        """
        override = self._step_overrides[step] if step in self._step_overrides else None
        if override is None:
            return default_impl()
        return override(self)

    def _set_psi_magnetic_axis(self) -> None:
        """Compute and store magnetic-axis location and PSI from selected vacuum nodes.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()

        psi_array = ApecAlgorithm._require_array(mesh.GetPointData(), "PSI")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsle(mesh)
        selection_array = ApecAlgorithm._require_array(mesh.GetPointData(), "selection_point")
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
        self.status["magnetic_axis_node_id"] = axis_index
        self.status["magnetic_axis"] = axis_point
        LOGGER.debug(
            "Magnetic axis at (%.4g, %.4g, %.4g), PSI_axis=%.6g.",
            axis_point[0], axis_point[1], axis_point[2], psi_axis,
        )

    def _reset_plasma_components(self) -> None:
        """Reset PLASMA/AMPS and Js fields inside selected vacuum-chamber elements.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()
        cell_data = mesh.GetCellData()

        js_array = ApecAlgorithm._require_array(cell_data, "Js")

        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.set_value_on_selection(mesh, "PLASMA", 1, type="elem")
        Selector.nsle(mesh)
        Selector.set_value_on_selection(mesh, "AMPS", 0.0, type="node")

        selection_array = ApecAlgorithm._require_array(cell_data, "selection_cell")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        js_np = vtknp.vtk_to_numpy(js_array)

        js_np[selected_ids, :] = 0.0
        js_array.Modified()

    def _identify_plasma_region(self) -> None:
        """Identify plasma elements from PSI thresholds and optional halo settings.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()

        psi_2_min = float(self.status["psi_2_min"])
        psi_min = float(self.status["psi_min"])
        is_limiter = self.status["separatrix_configuration"] == "limiter"
        halo_delta1 = float(self.input_parameters.get("halo_delta1", 0.0))
        t0_halo = float(self.input_parameters.get("t0_halo", 0.0))
        halo_requested = (
            is_limiter
            and halo_delta1 > 0.0
            and self.apec.time_step["time"] >= t0_halo
        )

        point_data = mesh.GetPointData()
        cell_data = mesh.GetCellData()
        self._halo_node_groups = ()

        if halo_requested and psi_2_min >= psi_min:
            raise ValueError(
                "psi_2_min must be lower than psi_min when identifying plasma region with halo."
            )
        has_halo = halo_requested and psi_2_min < psi_min

        if has_halo:
            lower_strict = np.nextafter(psi_2_min, np.inf)
            upper_strict = np.nextafter(psi_min, -np.inf)
            Selector.select(
                mesh,
                "PSI",
                (lower_strict, upper_strict),
                type="node",
                mode="r",
            )
            selection_array = point_data.GetArray("selection_point")
            halo_node_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
            boundary_node_ids, boundary_adjacency = self._get_halo_boundary_graph()
            halo_node_ids = np.intersect1d(
                halo_node_ids,
                boundary_node_ids,
                assume_unique=True,
            )
            halo_node_groups = self._connected_node_groups(
                halo_node_ids,
                boundary_adjacency,
            )
            points = vtknp.vtk_to_numpy(mesh.GetPoints().GetData())
            self._halo_node_groups = self._orient_halo_node_groups(
                halo_node_groups,
                points,
                self.status["magnetic_axis"],
                float(self.input_parameters["plasma_current"]),
            )

        # Strict PSIe < psi_2_min using an upper bound just below psi_2_min.
        upper_strict = np.nextafter(psi_2_min, -np.inf)
        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.select(
            mesh,
            "PSIe",
            (-np.inf, upper_strict),
            type="elem",
            mode="r",
        )

        selection_array = cell_data.GetArray("selection_cell")
        selected_elem_count = int(np.count_nonzero(vtknp.vtk_to_numpy(selection_array)))

        if selected_elem_count > 0:
            Selector.set_value_on_selection(mesh, "PLASMA", 0, type="elem")

        # 1. Select VACUUM_CHAMBER elements
        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        vac_cell_sel_array = ApecAlgorithm._require_array(cell_data, "selection_cell")
        vac_cell_mask = np.asarray(vtknp.vtk_to_numpy(vac_cell_sel_array), dtype=bool)

        if np.any(vac_cell_mask):
            psi_x = self.status.get("last_psi_x_point")
            axis_node_id = int(self.status.get("magnetic_axis_node_id", -1))

            if psi_x is not None and axis_node_id >= 0:
                psi_x = float(psi_x)
                psi_point_array = ApecAlgorithm._require_array(point_data, "PSI")
                psi_point = np.asarray(vtknp.vtk_to_numpy(psi_point_array), dtype=float).reshape(-1)

                # 2, 3 & 4. Run cell-based BFS requiring ALL nodes >= PSI_x AND >= 2 shared nodes
                connected_cell_ids = self._extract_main_plasma_region_cells(
                    mesh,
                    axis_node_id,
                    psi_point,
                    psi_x,
                    vac_cell_mask,
                )

                connected_cell_mask = np.zeros(mesh.GetNumberOfCells(), dtype=bool)
                if connected_cell_ids.size > 0:
                    connected_cell_mask[connected_cell_ids] = True

                # 5. Unselect candidate cells that are topologically disconnected
                cells = mesh.GetCells()
                offsets = np.asarray(vtknp.vtk_to_numpy(cells.GetOffsetsArray()), dtype=np.int64)
                connectivity = np.asarray(vtknp.vtk_to_numpy(cells.GetConnectivityArray()), dtype=np.int64)

                psi_axis = psi_point[axis_node_id]
                core_above = psi_axis >= psi_x

                remove_mask = np.zeros(mesh.GetNumberOfCells(), dtype=bool)
                vac_cell_ids = np.flatnonzero(vac_cell_mask)

                for cell_id in vac_cell_ids:
                    start = offsets[cell_id]
                    end = offsets[cell_id + 1]
                    cell_nodes = connectivity[start:end]
                    cell_psi = psi_point[cell_nodes]

                    # Candidate condition: ALL nodes satisfy the PSI constraint
                    all_nodes_valid = np.all(cell_psi >= psi_x) if core_above else np.all(cell_psi <= psi_x)

                    # If element is a valid candidate BUT unreachable via edge-sharing BFS, unselect it
                    if all_nodes_valid and not connected_cell_mask[cell_id]:
                        remove_mask[cell_id] = True

                # 1. Any VACUUM_CHAMBER element not reached by the core BFS is non-plasma
                remove_mask = vac_cell_mask & (~connected_cell_mask)

                # 2. Reset PLASMA flag on all disconnected divertor/SOL elements
                if np.any(remove_mask):
                    Selector.apply_boolean_selection_mask(
                        remove_mask,
                        cell_data,
                        "elem",
                        "s",
                    )
                    Selector.set_value_on_selection(mesh, "PLASMA", 0, type="elem")

        halo_group_count = len(self._halo_node_groups)
        LOGGER.debug(
            "Plasma region identified: %d elements excluded, %d halo group(s)%s.",
            selected_elem_count,
            halo_group_count,
            " (halo active)" if halo_group_count > 0 else "",
        )

    def _check_limiter(self) -> None:
        """Locate limiter contact point and associated PSI among selected boundary nodes,
        excluding divertor/private-flux nodes with inverted poloidal field rotation.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()

        psi_array = ApecAlgorithm._require_array(mesh.GetPointData(), "PSI")
        psi = np.asarray(vtknp.vtk_to_numpy(psi_array), dtype=float).reshape(-1)

        # 1. Initial selection of vacuum chamber external boundary nodes attached to PLASMA
        Selector.select(mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsel_ext(mesh, mode="s")

        selection_array = ApecAlgorithm._require_array(mesh.GetPointData(), "selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        if selected_ids.size == 0:
            raise ValueError(
                "No limiter candidates found on the external vacuum-chamber "
                "boundary attached to PLASMA elements."
            )

        # 2. Compute PSI point gradient to evaluate field rotation
        grad_filter = vtk.vtkGradientFilter()
        grad_filter.SetInputData(mesh)
        grad_filter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "PSI")
        grad_filter.SetResultArrayName("__apec_psi_grad_limiter__")
        grad_filter.Update()

        grad_mesh = grad_filter.GetOutput()
        grad_array = grad_mesh.GetPointData().GetArray("__apec_psi_grad_limiter__")
        if grad_array is None:
            raise RuntimeError("Failed to compute PSI gradient for limiter candidate filtering.")
        grad = np.asarray(vtknp.vtk_to_numpy(grad_array), dtype=float)

        # 3. Filter out divertor nodes where (r - r_axis) . grad(PSI) has reversed sign
        axis_node_id = int(self.status.get("magnetic_axis_node_id", -1))
        psi_x = self.status.get("last_psi_x_point")

        points = mesh.GetPoints()
        xyz = vtknp.vtk_to_numpy(points.GetData())

        if axis_node_id >= 0 and axis_node_id < mesh.GetNumberOfPoints():
            axis_xyz = xyz[axis_node_id]
            psi_axis = psi[axis_node_id]

            reference_psi = float(psi_x) if psi_x is not None else float(psi.max())
            core_sign = 1.0 if reference_psi >= psi_axis else -1.0

            candidate_xyz = xyz[selected_ids]
            candidate_grad = grad[selected_ids]

            # 2D radial vector from magnetic axis in (R, Z)
            vec_r = candidate_xyz[:, :2] - axis_xyz[:2]
            dot_prod = vec_r[:, 0] * candidate_grad[:, 0] + vec_r[:, 1] * candidate_grad[:, 1]

            # Keep only nodes matching core outward gradient alignment
            valid_mask = (core_sign * dot_prod) > 0.0
            if np.any(valid_mask):
                selected_ids = selected_ids[valid_mask]

        # 4. Find maximum PSI among filtered main-chamber boundary nodes
        psi_selected = psi[selected_ids]
        limiter_local_index = int(np.argmax(psi_selected))
        limiter_index = int(selected_ids[limiter_local_index])
        psi_limiter = float(psi[limiter_index])

        limiter_point = tuple(float(v) for v in xyz[limiter_index])

        self.status["last_psi_limiter"] = psi_limiter
        self.status["limiter_node_id"] = limiter_index
        self.status["limiter_point"] = limiter_point
        LOGGER.debug(
            "Limiter contact at (%.4g, %.4g, %.4g), PSI_limiter=%.6g.",
            limiter_point[0], limiter_point[1], limiter_point[2], psi_limiter,
        )

    @staticmethod
    def _compute_hessian_vmtk(
        ugrid: "vtk.vtkUnstructuredGrid",
        scalar_array_name: str,
    ) -> "vtk.vtkUnstructuredGrid":
        """Compute a 9-component Hessian tensor via two VMTK gradient passes."""
        if vtkvmtkUnstructuredGridGradientFilter is None:
            raise RuntimeError("VMTK gradient filter is unavailable.")

        grad_filter_1 = vtkvmtkUnstructuredGridGradientFilter()
        grad_filter_1_dyn = cast(Any, grad_filter_1)
        grad_filter_1.SetInputData(ugrid)
        grad_filter_1_dyn.SetInputArrayName(scalar_array_name)
        grad_filter_1_dyn.SetGradientArrayName("__apec_gradient_1__")
        grad_filter_1_dyn.SetQuadratureOrder(3)
        grad_filter_1.Update()

        grad_filter_2 = vtkvmtkUnstructuredGridGradientFilter()
        grad_filter_2_dyn = cast(Any, grad_filter_2)
        grad_filter_2.SetInputConnection(grad_filter_1.GetOutputPort())
        grad_filter_2_dyn.SetInputArrayName("__apec_gradient_1__")
        grad_filter_2_dyn.SetGradientArrayName("__apec_psi_hessian__")
        grad_filter_2_dyn.SetQuadratureOrder(3)
        grad_filter_2.Update()

        return cast("vtk.vtkUnstructuredGrid", grad_filter_2.GetOutput())

    def _check_x_point(self) -> None:
        """Locate the X-point as a PSI saddle-point candidate.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()

        # Extract only VACUUM_CHAMBER elements before Hessian evaluation.
        work_mesh = vtk.vtkUnstructuredGrid()
        work_mesh.DeepCopy(mesh)
        work_point_data = work_mesh.GetPointData()

        original_point_ids_np = np.arange(work_mesh.GetNumberOfPoints(), dtype=np.int64)
        original_point_ids_vtk = vtknp.numpy_to_vtk(original_point_ids_np, deep=True)
        original_point_ids_vtk.SetName("__apec_original_point_id__")
        work_point_data.AddArray(original_point_ids_vtk)

        Selector.select(work_mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        selection_cell_array = ApecAlgorithm._require_array(work_mesh.GetCellData(), "selection_cell")
        selected_cell_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_cell_array) == 1).astype(np.int64)
        if selected_cell_ids.size == 0:
            raise ValueError("No VACUUM_CHAMBER elements found for X-point search.")

        selected_cell_ids_vtk = vtk.vtkIdTypeArray()
        selected_cell_ids_vtk.SetNumberOfValues(selected_cell_ids.size)
        for i, cell_id in enumerate(selected_cell_ids):
            selected_cell_ids_vtk.SetValue(i, int(cell_id))

        selection_node = vtk.vtkSelectionNode()
        selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
        selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
        selection_node.SetSelectionList(selected_cell_ids_vtk)

        selection = vtk.vtkSelection()
        selection.AddNode(selection_node)

        extract = vtk.vtkExtractSelection()
        extract.SetInputData(0, work_mesh)
        extract.SetInputData(1, selection)
        extract.Update()

        vacuum_mesh = vtk.vtkUnstructuredGrid.SafeDownCast(extract.GetOutput())
        if vacuum_mesh is None or vacuum_mesh.GetNumberOfPoints() == 0:
            raise ValueError("Failed to extract VACUUM_CHAMBER subgrid for X-point search.")

        vacuum_point_data = vacuum_mesh.GetPointData()
        psi_array = ApecAlgorithm._require_array(vacuum_point_data, "PSI")

        try:
            hessian_output = self._compute_hessian_vmtk(vacuum_mesh, "PSI")
        except RuntimeError:
            grad_filter = vtk.vtkGradientFilter()
            grad_filter.SetInputData(vacuum_mesh)
            grad_filter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "PSI")
            grad_filter.SetResultArrayName("__apec_psi_gradient__")
            if hasattr(grad_filter, "SetComputeGradient"):
                grad_filter.SetComputeGradient(True)
            grad_filter.Update()

            # Portable Hessian computation fallback: H(PSI) = grad(grad(PSI)).
            hessian_filter = vtk.vtkGradientFilter()
            hessian_filter.SetInputConnection(grad_filter.GetOutputPort())
            hessian_filter.SetInputScalars(
                vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,
                "__apec_psi_gradient__",
            )
            hessian_filter.SetResultArrayName("__apec_psi_hessian__")
            if hasattr(hessian_filter, "SetComputeGradient"):
                hessian_filter.SetComputeGradient(True)
            hessian_filter.Update()
            hessian_output = hessian_filter.GetOutput()

        hessian_point_data = hessian_output.GetPointData()
        hessian_array = hessian_point_data.GetArray("__apec_psi_hessian__")
        if hessian_array is None:
            raise ValueError("vtkGradientFilter did not produce a PSI Hessian (gradient-of-gradient) point array.")

        gradient_array = hessian_point_data.GetArray("__apec_gradient_1__")
        if gradient_array is None:
            gradient_array = hessian_point_data.GetArray("__apec_psi_gradient__")
        if gradient_array is None:
            gradient_array = hessian_point_data.GetArray("Gradients")
        if gradient_array is None:
            raise ValueError("Could not find PSI gradient array to enforce X-point stationarity criterion.")

        Selector.select_all(vacuum_mesh, type="elem")
        Selector.nsle(vacuum_mesh)
        # Exclude midside nod
        # es: their PSI values are interpolated and can
        # produce false minima near the X-point with second-order elements.
        Selector.select(vacuum_mesh, "MidsideNodes", 1, type="node", mode="u")
        selection_array = ApecAlgorithm._require_array(vacuum_point_data, "selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        # Exclude external boundary nodes; keep only interior candidates.
        # Selector.nsel_ext(vacuum_mesh, mode="s")
        # boundary_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        # selected_ids = np.setdiff1d(selected_ids, boundary_ids, assume_unique=False)

        if selected_ids.size == 0:
            raise ValueError("No interior nodes selected for X-point saddle search.")

        psi = np.asarray(vtknp.vtk_to_numpy(psi_array), dtype=float).reshape(-1)
        hessian = np.asarray(vtknp.vtk_to_numpy(hessian_array), dtype=float)
        if hessian.ndim != 2 or hessian.shape[1] != 9:
            raise ValueError(
                "Unexpected PSI Hessian shape: expected (n_points, 9) for scalar PSI on 2D mesh."
            )

        # VTK/VMTK flattened tensor layout: [xx, xy, xz, yx, yy, yz, zx, zy, zz].
        hxx = hessian[:, 0]
        hxy = hessian[:, 1]
        hyx = hessian[:, 3]
        hyy = hessian[:, 4]

        det_h = hxx * hyy - hxy * hyx
        finite_mask = np.isfinite(det_h)
        saddle_mask = finite_mask & (det_h < 0.0)
        saddle_ids = selected_ids[saddle_mask[selected_ids]]

        if saddle_ids.size == 0:
            raise ValueError("No PSI saddle-point candidates found for the X point.")

        gradient = np.asarray(vtknp.vtk_to_numpy(gradient_array), dtype=float)
        if gradient.ndim != 2 or gradient.shape[1] < 2:
            raise ValueError("Unexpected PSI gradient shape while enforcing X-point stationarity criterion.")
        grad_mag = np.sqrt(np.square(gradient[:, 0]) + np.square(gradient[:, 1]))
        finite_grad_ids = saddle_ids[np.isfinite(grad_mag[saddle_ids])]
        if finite_grad_ids.size == 0:
            raise ValueError("No finite gradient magnitudes found among PSI saddle-point candidates.")

        x_point_local_id = int(finite_grad_ids[int(np.argmin(grad_mag[finite_grad_ids]))])
        original_point_ids_array = ApecAlgorithm._require_array(
            hessian_point_data,
            "__apec_original_point_id__",
        )
        original_point_ids = np.asarray(vtknp.vtk_to_numpy(original_point_ids_array), dtype=np.int64).reshape(-1)
        x_point_node_id = int(original_point_ids[x_point_local_id])

        psi_x = float(psi[x_point_local_id])

        points = mesh.GetPoints()
        xyz = vtknp.vtk_to_numpy(points.GetData())
        self.status["last_psi_x_point"] = psi_x
        self.status["x_point_node_id"] = x_point_node_id
        x_point = tuple(float(v) for v in xyz[x_point_node_id])
        self.status["x_point"] = x_point
        LOGGER.debug(
            "X-point at (%.4g, %.4g, %.4g), PSI_x=%.6g.",
            x_point[0], x_point[1], x_point[2], psi_x,
        )

    def _set_b_target_axis(self) -> None:
        """Interpolate B at the requested target axis point and store it in status.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()

        point_data = mesh.GetPointData()
        b_array = point_data.GetArray("B")
        b = np.asarray(vtknp.vtk_to_numpy(b_array), dtype=float)

        target_point = (
            float(self.input_parameters["target_central_axis_x"]),
            float(self.input_parameters["target_central_axis_y"]),
            float(self.input_parameters["target_central_axis_z"]),
        )
        target_point_seq = [target_point[0], target_point[1], target_point[2]]

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
            cell_id = int(locator.FindCell(target_point_seq))
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
                    target_point_seq,
                    closest_point,
                    cast(Any, sub_id),
                    pcoords,
                    cast(Any, dist2),
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
        LOGGER.debug("B at target axis: Bx=%.6g, By=%.6g.", target[0], target[1])

    def _compute_correction_currents(self) -> None:
        """Compute radial/vertical control-current corrections from target-axis B.

        Returns:
            None
        """
        bradx = float(self.input_parameters["bradx_value"])
        brady = float(self.input_parameters["brady_value"])
        bverx = float(self.input_parameters["bverx_value"])
        bvery = float(self.input_parameters["bvery_value"])
        btx, bty = self.status["b_target_axis"]

        det = bradx * bvery - bverx * brady
        if abs(det) <= np.finfo(float).eps:
            correction = (0.0, 0.0)
        else:
            dira = -(bverx * bty - bvery * btx) / det
            dive = -(brady * btx - bradx * bty) / det
            correction = (-dira, -dive)

        self.status["radial_correction_current"] = float(correction[0])
        self.status["vertical_correction_current"] = float(correction[1])
        LOGGER.debug(
            "Correction currents: radial=%.6g A, vertical=%.6g A.",
            correction[0], correction[1],
        )

    def _apply_control_currents(self) -> None:
        """Apply computed control currents to configured coil regions via Js updates.

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()
        control_coils = self._normalized_control_coils()
        if not control_coils:
            LOGGER.debug("No control coils configured; skipping control-current update.")
            return

        cell_data = mesh.GetCellData()
        js_array = ApecAlgorithm._require_array(cell_data, "Js")
        js_np = np.asarray(vtknp.vtk_to_numpy(js_array), dtype=float)
        radial_current = float(self.status["radial_correction_current"])
        vertical_current = float(self.status["vertical_correction_current"])

        selected_mask = np.zeros(js_np.shape[0], dtype=bool)
        correction_current_all = np.zeros(js_np.shape[0], dtype=float)
        for coil in control_coils:
            component = coil["component"]
            component_array = ApecAlgorithm._require_array(cell_data, component)
            component_mask = np.asarray(vtknp.vtk_to_numpy(component_array) == 1, dtype=bool)
            selected_mask |= component_mask

            correction_current = (
                float(coil["radial_gain"]) * radial_current
                + float(coil["vertical_gain"]) * vertical_current
            )
            correction_current_all[component_mask] += correction_current

        selected_ids = np.flatnonzero(selected_mask)
        if selected_ids.size == 0:
            LOGGER.debug("Configured control coils selected zero elements; skipping update.")
            return

        real_constant_array = ApecAlgorithm._require_array(cell_data, "CONSTANT")
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

        correction_current_density = correction_current_all[selected_ids] / selected_areas_np
        js_np[selected_ids, 2] = base_js_z[selected_ids] + correction_current_density

        js_array.Modified()

    def find_separatrix(self) -> None:
        """Run the separatrix-detection step (or its registered override).

        Returns:
            None
        """
        LOGGER.info("Step: find_separatrix%s.", " (override)" if self.STEP_FIND_SEPARATRIX in self._step_overrides else "")
        self._dispatch_step(
            self.STEP_FIND_SEPARATRIX,
            self._find_separatrix_default,
        )

    def _find_separatrix_default(self) -> None:
        """Default separatrix detection and plasma-region classification workflow.

        Returns:
            None
        """
        self._reset_plasma_components()
        self._set_psi_magnetic_axis()
        self._check_x_point()
        self._check_limiter()
        psi_limiter = float(self.status["last_psi_limiter"])
        psi_x_point = float(self.status["last_psi_x_point"])
        if psi_limiter > psi_x_point:
            self.status["psi_min"] = psi_limiter
            self.status["separatrix_configuration"] = "limiter"
        else:
            self.status["psi_min"] = psi_x_point
            self.status["separatrix_configuration"] = "X-point"
        halo_delta1 = float(self.input_parameters.get("halo_delta1", 0.0))
        halo_delta2 = float(self.input_parameters.get("halo_delta2", 0.0))
        t0_halo = float(self.input_parameters.get("t0_halo", 0.0))
        if halo_delta1 < 0.0:
            raise ValueError("halo_delta1 must be non-negative.")
        if (
            halo_delta1 > 0.0
            and self.apec.time_step["time"] >= t0_halo
            and self.status["separatrix_configuration"] == "limiter"
        ):
            halo_scale = halo_delta1 * (
                1.0 + (self.apec.time_step["time"] - t0_halo) * halo_delta2
            )
            self.status["psi_2_min"] = (
                self.status["psi_min"]
                - halo_scale * (self.status["psi_max"] - self.status["psi_min"])
            )
        else:
            self.status["psi_2_min"] = self.status["psi_min"]
        self._identify_plasma_region()
        LOGGER.info(
            "Separatrix: config=%s, PSI_min=%.6g, PSI_2_min=%.6g, PSI_axis=%.6g.",
            self.status["separatrix_configuration"],
            self.status["psi_min"],
            self.status["psi_2_min"],
            self.status["psi_magnetic_axis"],
        )

    def compute_jplasma_gradshafranov(self) -> None:
        """Run Grad-Shafranov current computation (or its registered override).

        Returns:
            None
        """
        LOGGER.info("Step: compute_jplasma_gradshafranov%s.", " (override)" if self.STEP_COMPUTE_JPLASMA in self._step_overrides else "")
        self._dispatch_step(
            self.STEP_COMPUTE_JPLASMA,
            self._compute_jplasma_gradshafranov_default,
        )

    def _compute_jplasma_gradshafranov_default(self) -> None:
        """Compute Grad-Shafranov plasma current profile and update Js (and halo AMPS).

        Returns:
            None
        """
        mesh = self._require_equilibrium_mesh()

        point_data = mesh.GetPointData()
        cell_data = mesh.GetCellData()

        js_array = ApecAlgorithm._require_array(cell_data, "Js")
        js_np = np.asarray(vtknp.vtk_to_numpy(js_array), dtype=float)

        psi_point_array = ApecAlgorithm._require_array(point_data, "PSI")
        psi_point = np.asarray(vtknp.vtk_to_numpy(psi_point_array), dtype=float).reshape(-1)

        psie_array = ApecAlgorithm._require_array(cell_data, "PSIe")
        psi_elem = np.asarray(vtknp.vtk_to_numpy(psie_array), dtype=float).reshape(-1)

        be_array = ApecAlgorithm._require_array(cell_data, "Be")
        be_np = np.asarray(vtknp.vtk_to_numpy(be_array), dtype=float)
        bx = np.asarray(be_np[:, 0], dtype=float)
        by = np.asarray(be_np[:, 1], dtype=float)

        volume_array = ApecAlgorithm._require_array(cell_data, "cell_volume")
        cell_volume = np.asarray(vtknp.vtk_to_numpy(volume_array), dtype=float).reshape(-1)

        centroids_array = ApecAlgorithm._require_array(cell_data, "cell_centroids")
        centroids = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
        elem_centroid_rad = np.asarray(centroids[:, 0], dtype=float)

        psi_magnetic_axis = float(self.status["psi_magnetic_axis"])
        psi_min = float(self.status["psi_2_min"])
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

        amps_array = point_data.GetArray("AMPS")
        halo_currents_enabled = (
            self.apec.mesh.GetMaxSpatialDimension() > 2
            and self.status["separatrix_configuration"] == "limiter"
            and float(self.input_parameters.get("halo_delta1", 0.0)) > 0.0
            and self.apec.time_step["time"]
            >= float(self.input_parameters.get("t0_halo", 0.0))
            and float(self.status["psi_2_min"]) < float(self.status["psi_min"])
        )
        if (
            halo_currents_enabled
            and amps_array is not None
            and self._halo_node_groups
        ):
            amps_np = vtknp.vtk_to_numpy(amps_array)
            halo_angle = float(self.input_parameters.get("extrusion_angle", 360.0))
            halo_angle_factor = halo_angle / 360.0
            halo_flux_den = psi_min - psi_magnetic_axis

            def apply_halo_currents(node_ids: np.ndarray) -> None:
                if node_ids.size <= 1:
                    return

                psi_seg = (psi_point[node_ids] - psi_magnetic_axis) / halo_flux_den
                psi_seg = np.clip(psi_seg, 0.0, 1.0)

                integrand = np.power(
                    np.maximum(1.0 - np.power(psi_seg, alpha_mp), 0.0),
                    alpha_np,
                )
                dpsi_seg = np.diff(psi_seg)
                integral = np.empty_like(psi_seg)
                integral[0] = 0.0
                integral[1:] = np.cumsum(
                    0.5 * (integrand[:-1] + integrand[1:]) * dpsi_seg
                )

                f_const = 2.0 * mu_0 * r0 * (1.0 - beta_profile) * halo_flux_den * lambda_scale
                f_profile = np.sqrt(np.maximum(f_const * integral + af0 * af0, 0.0))
                segment_currents = (
                    (2.0 * np.pi / mu_0)
                    * (f_profile[:-1] - f_profile[1:])
                    * halo_angle_factor
                )
                nodal_segment_currents = 0.5 * segment_currents
                np.add.at(amps_np, node_ids[:-1], nodal_segment_currents)
                np.add.at(amps_np, node_ids[1:], nodal_segment_currents)

            for halo_node_group in self._halo_node_groups:
                apply_halo_currents(halo_node_group)
            amps_array.Modified()

        self.status["total_plasma_current_before_scaling"] = total_plasma_current
        self.status["plasma_current_scaling_lambda"] = float(lambda_scale)
        self.status["plasma_element_count"] = int(selected_ids.size)
        LOGGER.info(
            "Grad-Shafranov: %d plasma elements, I_plasma_raw=%.6g A, lambda=%.6g.",
            selected_ids.size, total_plasma_current, lambda_scale,
        )

    def control_currents(self) -> None:
        """Run control-current update step (or its registered override).

        Returns:
            None
        """
        LOGGER.info("Step: control_currents%s.", " (override)" if self.STEP_CONTROL_CURRENTS in self._step_overrides else "")
        self._dispatch_step(
            self.STEP_CONTROL_CURRENTS,
            self._control_currents_default,
        )

    def _control_currents_default(self) -> None:
        """Default control-current pipeline: sample B, compute correction, apply Js.

        Returns:
            None
        """
        self._set_b_target_axis()
        self._compute_correction_currents()
        self._apply_control_currents()
        LOGGER.info(
            "Control currents applied: radial=%.6g A, vertical=%.6g A.",
            self.status["radial_correction_current"],
            self.status["vertical_correction_current"],
        )

    def check_convergence(self) -> bool:
        """Run convergence check step (or its registered override).

        Returns:
            bool: True when convergence criterion is satisfied.
        """
        result = self._dispatch_step(
            self.STEP_CHECK_CONVERGENCE,
            self._check_convergence_default,
        )
        return bool(result)

    def _check_convergence_default(self) -> bool:
        """Check Picard convergence using max absolute PSI change on selected nodes.

        Returns:
            bool: True when convergence criterion is satisfied.
        """
        mesh = self._require_equilibrium_mesh()

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
        LOGGER.info(
            "Convergence check: flux_norm=%.6g, tolerance=%.6g -> %s.",
            flux_norm, tolerance, "CONVERGED" if converged else "not converged",
        )
        return converged

    # Backward-compatible alias.
    def check_poloidal_flux_convergence(self) -> bool:
        """Backward-compatible alias for :meth:`check_convergence`.

        Returns:
            bool: Convergence status.
        """
        return self.check_convergence()

    def _toroidal_average(self, psi_only: bool = False) -> None:
        """Average toroidally mapped values from evolution mesh into equilibrium mesh.

        Args:
            psi_only: When True, average only PSI values.

        Returns:
            None
        """
        LOGGER.debug("Toroidal average: psi_only=%s.", psi_only)
        equilibrium_mesh = self._require_equilibrium_mesh()
        evolution_mesh = self.apec.mesh
        self._set_toroidal_slices()
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

        eq_psie_array = eq_cell_data.GetArray("PSIe")
        evol_psie_array = evol_cell_data.GetArray("PSIe")
        if (
            eq_psie_array is None
            or evol_psie_array is None
        ):
            raise ValueError(
                "Required arrays for toroidal_average are missing "
                "(PSIe)."
            )

        eq_psie = np.asarray(vtknp.vtk_to_numpy(eq_psie_array), dtype=float).reshape(-1)
        evol_psie = np.asarray(vtknp.vtk_to_numpy(evol_psie_array), dtype=float).reshape(-1)

        src_elems = np.asarray(elem_slices[:, 0], dtype=np.int64)
        dst_elems = np.asarray(elem_slices[:, 1:], dtype=np.int64)

        psie_values = np.mean(evol_psie[dst_elems], axis=1)
        eq_psie[src_elems] = psie_values
        eq_psie_array.Modified()

    def _toroidal_projection_current(self) -> None:
        """Project equilibrium current state back onto the static/evolution mesh.

        Returns:
            None
        """
        from .coordinate_system import YAxisCylindricalCoordinateSystem

        static_mesh = self.apec.mesh
        equilibrium_mesh = self._require_equilibrium_mesh()
        self._set_toroidal_slices()
        node_slices = self._toroidal_slices_nodes
        elem_slices = self._toroidal_slices_elems

        static_amps_array = static_mesh.GetPointData().GetArray("AMPS")
        equilibrium_amps_array = equilibrium_mesh.GetPointData().GetArray("AMPS")
        static_js_array = static_mesh.GetCellData().GetArray("Js")
        equilibrium_js_array = equilibrium_mesh.GetCellData().GetArray("Js")

        cylindrical = YAxisCylindricalCoordinateSystem()
        static_amps = (
            np.asarray(vtknp.vtk_to_numpy(static_amps_array), dtype=float).reshape(-1)
            if static_amps_array is not None and equilibrium_amps_array is not None
            else None
        )
        equilibrium_amps = (
            np.asarray(vtknp.vtk_to_numpy(equilibrium_amps_array), dtype=float).reshape(-1)
            if static_amps_array is not None and equilibrium_amps_array is not None
            else None
        )
        static_js = np.asarray(vtknp.vtk_to_numpy(static_js_array), dtype=float)
        equilibrium_js = np.asarray(vtknp.vtk_to_numpy(equilibrium_js_array), dtype=float)

        src_nodes = np.asarray(node_slices[:, 0], dtype=np.int64)
        dst_nodes = np.asarray(node_slices[:, 1:], dtype=np.int64)
        src_elems = np.asarray(elem_slices[:, 0], dtype=np.int64)
        dst_elems = np.asarray(elem_slices[:, 1:], dtype=np.int64)
        if static_amps is not None and equilibrium_amps is not None:
            toroidal_node_count = dst_nodes.shape[1]
            if toroidal_node_count <= 0:
                raise ValueError("Toroidal node map has no target node columns.")

            flat_dst_nodes = dst_nodes.reshape(-1)
            valid_node_mask = flat_dst_nodes >= 0
            flat_dst_nodes = flat_dst_nodes[valid_node_mask]
            projected_amps = np.broadcast_to(
                equilibrium_amps[src_nodes, np.newaxis] / float(toroidal_node_count),
                dst_nodes.shape,
            ).reshape(-1)[valid_node_mask]

            static_amps[np.unique(flat_dst_nodes)] = 0.0
            np.add.at(static_amps, flat_dst_nodes, projected_amps)
            static_amps_array.Modified()

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

    def _append_iteration_diagnostics(
        self,
        converged: bool,
    ) -> None:
        """Append one solver-iteration diagnostics row to the CSV diagnostics file.

        Args:
            converged: Current convergence flag.

        Returns:
            None
        """
        if not isinstance(converged, (bool, np.bool_)):
            raise TypeError("converged must be a boolean value.")

        self.status["converged"] = bool(converged)
        if self.apec._world_rank() != 0:
            return

        row = self._build_diagnostics_row(converged)
        row_df = pd.DataFrame([row])
        csv_path = self._diagnostics_csv_path()
        write_header = (not os.path.exists(csv_path)) or os.path.getsize(csv_path) == 0
        row_df.to_csv(csv_path, mode="a", header=write_header, index=False)

    def _calculate_poloidal_flux(self) -> None:
        """Update Ae/PSI/PSIe fields from current vector-potential solution arrays.

        Returns:
            None
        """
        mesh = self.apec.mesh
        if mesh.GetPointData().GetArray("A") is None:
            raise ValueError("Required PointData array 'A' is missing before poloidal flux calculation.")

        ModelEdit.element_solution(mesh, "A")
        ModelEdit.poloidal_flux(
            mesh=mesh,
            coordinate_system=YAxisCylindricalCoordinateSystem(),
        )

    def equilibrium_solve(self) -> bool:
        """Run the equilibrium Picard loop until convergence or iteration limit.

        Returns:
            bool: True when convergence is reached.
        """
        LOGGER.info("Starting equilibrium solve.")
        if self.apec.solver is None:
            raise RuntimeError("No solver session available. Call load_mesh_into_solver() before equilibrium_solve().")

        if self._first_absolute_iteration:
            self._clear_diagnostics_file()
        self.equilibrium_mesh = self.apec.mesh
        if self._first_absolute_iteration:
            self._cache_halo_boundary_graph()
        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )
        self.status["converged"] = False

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1
            LOGGER.info(
                "Equilibrium Picard iteration %d/%d.",
                self.status["picard_iteration"],
                max_picard_iterations,
            )

            if not self._first_absolute_iteration:
                self.find_separatrix()
                self.compute_jplasma_gradshafranov()
                self.control_currents()

            self.apec.solve()
            self._calculate_poloidal_flux()

            converged = self.check_convergence()
            # self.apec.write_output()
            self._append_iteration_diagnostics(converged)
            if converged:
                self.status["converged"] = True
                LOGGER.info(
                    "Equilibrium solve converged at iteration %d (flux_norm=%s).",
                    self.status["picard_iteration"],
                    self.status["flux_norm"],
                )
                return True
            self._first_absolute_iteration = False

        self.status["converged"] = False
        LOGGER.error(
            "Equilibrium solve did not converge after %d iterations (last flux_norm=%s).",
            max_picard_iterations,
            self.status["flux_norm"],
        )
        raise RuntimeError(
            "Equilibrium solve did not converge after "
            f"{max_picard_iterations} Picard iterations "
            f"(last flux metric: {self.status['flux_norm']})."
        )

    def evolution_solve(self) -> bool:
        """Run the evolution Picard loop, including toroidal transfer operations.

        Returns:
            bool: True when evolution solve completes successfully.
        """
        LOGGER.info("Starting evolution solve.")
        if self.apec.solver is None:
            raise RuntimeError("No solver session available. Call load_mesh_into_solver() before evolution_solve().")

        equilibrium_mesh = self._require_equilibrium_mesh()

        if self._first_absolute_iteration:
            self._clear_diagnostics_file()
            Selector.select(equilibrium_mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
            Selector.nsle(equilibrium_mesh, mode="s")
            Selector.select_all(self.apec.mesh)
            self._set_toroidal_slices()
            self._cache_halo_boundary_graph()

        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1
            LOGGER.info(
                "Evolution Picard iteration %d/%d.",
                self.status["picard_iteration"],
                max_picard_iterations,
            )

            if not self._first_absolute_iteration:
                self._toroidal_average()
                self.find_separatrix()
                self.compute_jplasma_gradshafranov()

                self._toroidal_projection_current()
                self.apec.solve()
                self._calculate_poloidal_flux()
                # self.apec.write_output()

                self._toroidal_average(psi_only=True)

                converged = self.check_convergence()
                self._append_iteration_diagnostics(
                    converged,
                )
                if converged:
                    self.status["converged"] = True
                    LOGGER.info(
                        "Evolution solve converged at iteration %d (flux_norm=%s).",
                        self.status["picard_iteration"],
                        self.status["flux_norm"],
                    )
                    return True
            else:
                self._toroidal_projection_current()
                self.apec.solve()
                self._calculate_poloidal_flux()
                # self.apec.write_output()
                self.check_convergence()
                self._first_absolute_iteration = False
                self.status["converged"] = True
                LOGGER.info("Evolution solve completed initial absolute iteration.")
                return True

        self.status["converged"] = False
        LOGGER.error(
            "Evolution solve did not converge after %d iterations (last flux_norm=%s).",
            max_picard_iterations,
            self.status["flux_norm"],
        )
        raise RuntimeError(
            "Evolution solve did not converge after "
            f"{max_picard_iterations} Picard iterations "
            f"(last flux metric: {self.status['flux_norm']})."
        )

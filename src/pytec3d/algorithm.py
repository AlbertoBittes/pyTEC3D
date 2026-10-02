from __future__ import annotations

import copy
import os
from typing import Any, Dict, List, Optional, TYPE_CHECKING, Tuple, cast

import numpy as np
import pandas as pd
import vtk

from pytec3d.selector import Selector
from pytec3d.model_edit import ModelEdit
from pytec3d.coordinate_system import CylindricalCoordinateSystem, YAxisCylindricalCoordinateSystem
from pytec3d.logging_utils import get_logger
from pytec3d.constants import (
    DEFAULT_pyTEC3D_ALGORITHM_INPUT_PARAMETERS,
    DEFAULT_pyTEC3D_ALGORITHM_STATUS,
)

import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    from pytec3d_solver.solver import pyTEC3DSolver


LOGGER = get_logger("pytec3d.algorithm")

class pyTEC3DAlgorithm:
    """
    Python scaffold of the legacy pyTEC3D equilibrium/evolution analysis objects.

    Notes:
        - This class mirrors the old C++ analysis API shape.
        - Methods that previously depended on libMesh internals are intentionally
          left as VTK-grid placeholders and can be refined incrementally.
        - Only time-step logic is supported (no substeps).
    """

    def __init__(
        self,
        pytec3d: "pyTEC3DSolver",
        input_parameters: Optional[Dict[str, Any]] = None,
    ) -> None:
        """Initialize algorithm state and optional parameter overrides.

        Args:
            pytec3d: Owning pyTEC3D object.
            input_parameters: Optional input-parameter overrides.

        Returns:
            None
        """
        self.pytec3d: "pyTEC3DSolver" = pytec3d
        self.input_parameters: Dict[str, Any] = copy.deepcopy(
            DEFAULT_pyTEC3D_ALGORITHM_INPUT_PARAMETERS
        )
        if input_parameters:
            self.input_parameters.update(input_parameters)

        self.status: Dict[str, Any] = copy.deepcopy(DEFAULT_pyTEC3D_ALGORITHM_STATUS)

        self._first_absolute_iteration: bool = True
        self._previous_flux_vector: Optional[np.ndarray] = None
        self._old_control_current_density_base: Optional[np.ndarray] = None
        self._b_target_axis_cache: Optional[
            Tuple[Tuple[Any, ...], np.ndarray, np.ndarray]
        ] = None
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)
        self._halo_node_groups: Tuple[np.ndarray, ...] = ()
        self._halo_boundary_graph_cache: Optional[
            Tuple[np.ndarray, Dict[int, np.ndarray]]
        ] = None
        self._equilibrium_mesh: Optional[vtk.vtkUnstructuredGrid] = None
        LOGGER.debug(
            "pyTEC3DAlgorithm initialized%s.",
            " with parameter overrides" if input_parameters else "",
        )

    def reset(self) -> None:
        """Reset iterative state and cached data to default values.

        Returns:
            None
        """
        self.status = copy.deepcopy(DEFAULT_pyTEC3D_ALGORITHM_STATUS)
        self._first_absolute_iteration = True
        self._previous_flux_vector = None
        self._old_control_current_density_base = None
        self._b_target_axis_cache = None
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)
        self._halo_node_groups = ()
        self._halo_boundary_graph_cache = None
        LOGGER.debug("pyTEC3DAlgorithm state reset.")

    @property
    def vacuum_chamber_mesh(self) -> vtk.vtkUnstructuredGrid:
        """Return the extracted vacuum-chamber equilibrium mesh.

        Returns:
            vtk.vtkUnstructuredGrid: Validated equilibrium mesh.
        """
        if self._equilibrium_mesh is None:
            raise ValueError("equilibrium_mesh is not set.")
        return self._equilibrium_mesh

    @vacuum_chamber_mesh.setter
    def vacuum_chamber_mesh(self, value: vtk.vtkUnstructuredGrid) -> None:
        """Set the equilibrium mesh from a source mesh by extracting VACUUM_CHAMBER cells."""

        if self._equilibrium_mesh is not None:
            raise ValueError("equilibrium_mesh is already set.")

        source_mesh = vtk.vtkUnstructuredGrid()
        source_mesh.DeepCopy(value)

        Selector.select(source_mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        selection_cell_array = self._require_array(source_mesh.GetCellData(), "selection_cell")
        selected_cell_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_cell_array) == 1).astype(np.int64)
        if selected_cell_ids.size == 0:
            raise ValueError("No VACUUM_CHAMBER elements found while setting equilibrium_mesh.")

        selected_cell_ids_vtk = vtknp.numpy_to_vtkIdTypeArray(selected_cell_ids, deep=True)

        selection_node = vtk.vtkSelectionNode()
        selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
        selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
        selection_node.SetSelectionList(selected_cell_ids_vtk)

        selection = vtk.vtkSelection()
        selection.AddNode(selection_node)

        extract = vtk.vtkExtractSelection()
        extract.SetInputData(0, source_mesh)
        extract.SetInputData(1, selection)
        extract.Update()

        extracted = vtk.vtkUnstructuredGrid.SafeDownCast(extract.GetOutput())
        if extracted is None or extracted.GetNumberOfCells() == 0:
            raise ValueError("Failed to extract VACUUM_CHAMBER subgrid for equilibrium_mesh.")

        Selector.select_all(extracted, type="elem")
        Selector.nsle(extracted, mode="s")

        self._equilibrium_mesh = extracted
        model_mesh = self.pytec3d.model.mesh
        self._toroidal_slices_nodes = np.zeros((0, 0), dtype=np.int64)
        self._toroidal_slices_elems = np.zeros((0, 0), dtype=np.int64)
        self._halo_boundary_graph_cache = None

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
        extrusion_angle: float = 360.0,
        n_divisions: Optional[int] = None,
        coordinate_system: Optional["CylindricalCoordinateSystem"] = None,
    ) -> np.ndarray:
        """
        Build a toroidal projection map from a 2D input mesh to a target mesh.

        The first column stores source IDs. The remaining columns store the
        closest target IDs corresponding to the expected toroidal stations.

        For a 3D target:
            - node mapping -> n_divisions + 1 stations
            - element mapping -> n_divisions stations

        For a 2D target:
            - one target match is returned.

        Matching is performed using a vtkStaticPointLocator in Cartesian space.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of
                revolution.
            target_mesh: Target unstructured grid.
            location: "node" to map point IDs or "elem" to map cell IDs.
            extrusion_angle: Total toroidal extrusion angle in degrees.
            n_divisions: Number of toroidal element divisions.
            coordinate_system: Cylindrical coordinate system used for the
                projection.

        Returns:
            np.ndarray:
                Integer map with shape (num_input, 1 + num_matches).
        """
        import vtk

        # ------------------------------------------------------------------
        # Validation
        # ------------------------------------------------------------------

        if location not in ("node", "elem"):
            raise ValueError("location must be 'node' or 'elem'.")

        if extrusion_angle <= 0.0:
            raise ValueError("extrusion_angle must be positive.")

        if (
            input_mesh.GetNumberOfCells() > 0
            and input_mesh.GetMaxSpatialDimension() > 2
        ):
            raise ValueError("Input mesh must be 2D.")

        if coordinate_system is None:
            coordinate_system = YAxisCylindricalCoordinateSystem()

        if coordinate_system.type != "cylindrical":
            raise ValueError(
                "coordinate_system must be cylindrical for toroidal projection."
            )

        # ------------------------------------------------------------------
        # Determine target dimensionality and toroidal stations
        # ------------------------------------------------------------------

        target_is_3d = target_mesh.GetMaxSpatialDimension() > 2

        if target_is_3d:
            if n_divisions is None:
                raise ValueError(
                    "n_divisions must be provided for 3D target meshes."
                )

            n_divisions = int(n_divisions)

            if n_divisions <= 0:
                raise ValueError(
                    "n_divisions must be a positive integer."
                )

            if location == "node":
                num_matches = n_divisions + 1
                expected_theta_deg = np.linspace(
                    0.0,
                    float(extrusion_angle),
                    num_matches,
                    dtype=np.float64,
                )
            else:
                num_matches = n_divisions
                expected_theta_deg = (
                    np.arange(
                        n_divisions,
                        dtype=np.float64,
                    )
                    + 0.5
                ) * (
                    float(extrusion_angle) / n_divisions
                )

        else:
            num_matches = 1
            expected_theta_deg = np.zeros(
                1,
                dtype=np.float64,
            )

        # ------------------------------------------------------------------
        # Convert meshes to cylindrical coordinates
        # ------------------------------------------------------------------

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
            return np.empty(
                (0, 1 + num_matches),
                dtype=np.int64,
            )

        if target_local.size == 0:
            raise ValueError(
                "Target mesh has no points/cells to map."
            )

        input_local = np.asarray(
            input_local,
            dtype=np.float64,
        )
        target_local = np.asarray(
            target_local,
            dtype=np.float64,
        )

        input_ids = np.asarray(input_ids)
        target_ids = np.asarray(target_ids)

        input_radius = input_local[:, 0]
        input_axial = input_local[:, 2]

        target_radius = target_local[:, 0]
        target_theta_deg = target_local[:, 1]
        target_axial = target_local[:, 2]

        num_input = input_radius.size
        target_count = target_radius.size

        # ------------------------------------------------------------------
        # Convert target cylindrical coordinates to Cartesian coordinates.
        #
        # Y-axis cylindrical system:
        #
        #   x = R cos(theta)
        #   y = Z
        #   z = R sin(theta)
        #
        # The locator therefore works in the actual 3D spatial geometry.
        # ------------------------------------------------------------------

        target_theta_rad = np.deg2rad(target_theta_deg)

        target_xyz = np.empty(
            (target_count, 3),
            dtype=np.float64,
        )

        target_xyz[:, 0] = (
            target_radius * np.cos(target_theta_rad)
        )
        target_xyz[:, 1] = target_axial
        target_xyz[:, 2] = (
            target_radius * np.sin(target_theta_rad)
        )

        # ------------------------------------------------------------------
        # Build static locator
        # ------------------------------------------------------------------

        vtk_target_points = vtknp.numpy_to_vtk(
            target_xyz,
            deep=True,
        )

        locator_points = vtk.vtkPoints()
        locator_points.SetData(vtk_target_points)

        locator_dataset = vtk.vtkPolyData()
        locator_dataset.SetPoints(locator_points)

        locator = vtk.vtkStaticPointLocator()
        locator.SetDataSet(locator_dataset)
        locator.AutomaticOn()
        locator.BuildLocator()

        # ------------------------------------------------------------------
        # Generate all query points in one NumPy operation.
        #
        # Shape:
        #
        #     (num_input, num_matches, 3)
        #
        # Flatten afterwards so the locator queries can be performed in a
        # simple contiguous loop.
        # ------------------------------------------------------------------

        theta_rad = np.deg2rad(expected_theta_deg)

        query_x = (
            input_radius[:, None]
            * np.cos(theta_rad)[None, :]
        )

        query_y = input_axial[:, None]

        query_z = (
            input_radius[:, None]
            * np.sin(theta_rad)[None, :]
        )

        num_queries = num_input * num_matches

        query_points = np.empty(
            (num_queries, 3),
            dtype=np.float64,
        )

        query_points[:, 0] = query_x.ravel()
        query_points[:, 1] = np.broadcast_to(
            query_y,
            (num_input, num_matches),
        ).ravel()
        query_points[:, 2] = query_z.ravel()

        # ------------------------------------------------------------------
        # Query closest target point.
        #
        # vtkStaticPointLocator performs the actual spatial search in C++.
        # Only the scalar lookup crosses the Python/VTK boundary.
        # ------------------------------------------------------------------

        target_indices = np.empty(
            num_queries,
            dtype=np.int64,
        )

        for query_idx in range(num_queries):
            target_indices[query_idx] = locator.FindClosestPoint(
                query_points[query_idx]
            )

        # ------------------------------------------------------------------
        # Convert target locator indices to actual target IDs.
        # ------------------------------------------------------------------

        matched_ids = target_ids[target_indices]

        # ------------------------------------------------------------------
        # Assemble result.
        # ------------------------------------------------------------------

        mapping = np.empty(
            (num_input, 1 + num_matches),
            dtype=np.int64,
        )

        mapping[:, 0] = input_ids.astype(
            np.int64,
            copy=False,
        )

        mapping[:, 1:] = matched_ids.reshape(
            num_input,
            num_matches,
        ).astype(
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

        equilibrium_mesh = self.vacuum_chamber_mesh

        target_dimension = self.pytec3d.model.mesh.GetMaxSpatialDimension()
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

        Selector.select(self.pytec3d.model.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
        Selector.nsle(self.pytec3d.model.mesh, mode="s")
        Selector.select(equilibrium_mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")

        self._toroidal_slices_nodes = self.build_toroidal_slices(
            equilibrium_mesh,
            self.pytec3d.model.mesh,
            location="node",
            extrusion_angle=extrusion_angle,
            n_divisions=n_toroidal_divisions,
        )
        self._toroidal_slices_elems = self.build_toroidal_slices(
            equilibrium_mesh,
            self.pytec3d.model.mesh,
            location="elem",
            extrusion_angle=extrusion_angle,
            n_divisions=n_toroidal_divisions,
        )

    def _diagnostics_csv_path(self) -> str:
        """Return diagnostics CSV output path under the __pytec3d__ folder.

        Returns:
            str: Diagnostics CSV file path.
        """
        diagnostics_dir = os.path.join(os.getcwd(), "__pytec3d__")
        os.makedirs(diagnostics_dir, exist_ok=True)
        return os.path.join(diagnostics_dir, "pytec3d_iteration_diagnostics.csv")

    def _clear_diagnostics_file(self) -> None:
        """Remove diagnostics CSV on rank 0 before a new run.

        Returns:
            None
        """
        if self.pytec3d._world_rank() != 0:
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
            "time": float(self.pytec3d.model.time_step["time"]),
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

    def _compute_psi_derivatives(self) -> None:
        """Compute and cache PSI gradient and Hessian on the equilibrium mesh."""
        mesh = self.vacuum_chamber_mesh

        grad_filter = vtk.vtkGradientFilter()
        grad_filter.SetInputData(mesh)
        grad_filter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "PSI")
        grad_filter.SetResultArrayName("__pytec3d_psi_gradient__")
        if hasattr(grad_filter, "SetComputeGradient"):
            grad_filter.SetComputeGradient(True)
        grad_filter.Update()

        grad_output = grad_filter.GetOutput()
        gradient_array = grad_output.GetPointData().GetArray("__pytec3d_psi_gradient__")
        if gradient_array is None:
            raise RuntimeError("Failed to compute PSI gradient on equilibrium mesh.")

        hessian_filter = vtk.vtkGradientFilter()
        hessian_filter.SetInputConnection(grad_filter.GetOutputPort())
        hessian_filter.SetInputScalars(
            vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,
            "__pytec3d_psi_gradient__",
        )
        hessian_filter.SetResultArrayName("__pytec3d_psi_hessian__")
        if hasattr(hessian_filter, "SetComputeGradient"):
            hessian_filter.SetComputeGradient(True)
        hessian_filter.Update()

        hessian_output = hessian_filter.GetOutput()
        hessian_array = hessian_output.GetPointData().GetArray("__pytec3d_psi_hessian__")
        if hessian_array is None:
            raise RuntimeError("Failed to compute PSI Hessian on equilibrium mesh.")

        mesh_point_data = mesh.GetPointData()
        existing_grad = mesh_point_data.GetArray("__pytec3d_psi_gradient__")
        if existing_grad is not None:
            mesh_point_data.RemoveArray("__pytec3d_psi_gradient__")
        existing_hess = mesh_point_data.GetArray("__pytec3d_psi_hessian__")
        if existing_hess is not None:
            mesh_point_data.RemoveArray("__pytec3d_psi_hessian__")

        mesh_point_data.AddArray(gradient_array)
        mesh_point_data.AddArray(hessian_array)

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
        mesh = self.vacuum_chamber_mesh
        Selector.select_all(mesh)
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

        offsets = vtknp.vtk_to_numpy(offsets_array)
        connectivity = vtknp.vtk_to_numpy(connectivity_array)

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

    def _set_psi_magnetic_axis(self) -> None:
        """Compute and store magnetic-axis location and PSI from all mesh nodes.

        Returns:
            None
        """
        mesh = self.vacuum_chamber_mesh

        psi_array = pyTEC3DAlgorithm._require_array(mesh.GetPointData(), "PSI")

        psi = vtknp.vtk_to_numpy(psi_array)
        axis_index = int(np.argmax(psi))
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
        """Reset PLASMA/AMPS and Js fields over the full mesh.

        Returns:
            None
        """
        mesh = self.vacuum_chamber_mesh
        cell_data = mesh.GetCellData()
        point_data = mesh.GetPointData()

        plasma_array = pyTEC3DAlgorithm._require_array(cell_data, "PLASMA")
        amps_array = pyTEC3DAlgorithm._require_array(point_data, "AMPS")
        js_array = pyTEC3DAlgorithm._require_array(cell_data, "Js")

        plasma_np = vtknp.vtk_to_numpy(plasma_array)
        amps_np = vtknp.vtk_to_numpy(amps_array)
        js_np = vtknp.vtk_to_numpy(js_array)

        plasma_np[...] = 1
        amps_np[...] = 0.0
        js_np[...] = 0.0

        plasma_array.Modified()
        amps_array.Modified()
        js_array.Modified()

    def _identify_plasma_region(self) -> None:
        """Identify plasma elements from PSI thresholds and optional halo settings.

        Returns:
            None
        """
        mesh = self.vacuum_chamber_mesh

        psi_2_min = float(self.status["psi_2_min"])
        psi_min = float(self.status["psi_min"])
        is_limiter = self.status["separatrix_configuration"] == "limiter"
        halo_delta1 = float(self.input_parameters.get("halo_delta1", 0.0))
        t0_halo = float(self.input_parameters.get("t0_halo", 0.0))
        halo_requested = (
            is_limiter
            and halo_delta1 > 0.0
            and self.pytec3d.model.time_step["time"] >= t0_halo
        )

        point_data = mesh.GetPointData()
        cell_data = mesh.GetCellData()
        self._halo_node_groups = ()

        if halo_requested and psi_2_min >= psi_min:
            raise ValueError(
                "psi_2_min must be lower than psi_min when identifying plasma region with halo."
            )
        has_halo = halo_requested and psi_2_min < psi_min

        selected_elem_count = 0

        vac_cell_mask = np.ones(mesh.GetNumberOfCells(), dtype=bool)

        if np.any(vac_cell_mask):
            psi_x = self.status.get("last_psi_x_point")
            axis_node_id = int(self.status.get("magnetic_axis_node_id", -1))

            if psi_x is not None and axis_node_id >= 0:
                psi_x = float(psi_x)
                psi_point_array = pyTEC3DAlgorithm._require_array(point_data, "PSI")
                psi_point = vtknp.vtk_to_numpy(psi_point_array)

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
                offsets = vtknp.vtk_to_numpy(cells.GetOffsetsArray())
                connectivity = vtknp.vtk_to_numpy(cells.GetConnectivityArray())

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
                    selected_elem_count = int(np.count_nonzero(remove_mask))
                    Selector.apply_boolean_selection_mask(
                        remove_mask,
                        cell_data,
                        "elem",
                        "s",
                    )
                    Selector.set_value_on_selection(mesh, "PLASMA", 0, type="elem")

        if has_halo:
            lower_strict = np.nextafter(psi_2_min, np.inf)
            upper_strict = np.nextafter(psi_min, -np.inf)
            Selector.select(
                mesh,
                "PSI",
                (lower_strict, upper_strict),
                type="node",
                mode="s",
            )
            # Keep only PSI-band nodes attached to currently active plasma elements.
            Selector.select(mesh, "PLASMA", 1, type="elem", mode="s")
            Selector.nsle(mesh, mode="r")

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
        mesh = self.vacuum_chamber_mesh

        psi_array = pyTEC3DAlgorithm._require_array(mesh.GetPointData(), "PSI")
        psi = vtknp.vtk_to_numpy(psi_array)

        # 1. Initial selection of external boundary nodes attached to PLASMA.
        Selector.select_all(mesh, type="elem")
        Selector.nsel_ext(mesh, mode="s")

        selection_array = pyTEC3DAlgorithm._require_array(mesh.GetPointData(), "selection_point")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)
        if selected_ids.size == 0:
            raise ValueError(
                "No limiter candidates found on the external vacuum-chamber "
                "boundary attached to PLASMA elements."
            )

        # 2. Use precomputed PSI point gradient to evaluate field rotation.
        grad_array = pyTEC3DAlgorithm._require_array(mesh.GetPointData(), "__pytec3d_psi_gradient__")
        grad = vtknp.vtk_to_numpy(grad_array)

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

    def _check_x_point(self) -> None:
        """Locate the X-point as a PSI saddle-point candidate.

        Returns:
            None
        """
        mesh = self.vacuum_chamber_mesh
        point_data = mesh.GetPointData()
        psi_array = pyTEC3DAlgorithm._require_array(point_data, "PSI")
        hessian_array = pyTEC3DAlgorithm._require_array(point_data, "__pytec3d_psi_hessian__")
        gradient_array = pyTEC3DAlgorithm._require_array(point_data, "__pytec3d_psi_gradient__")

        psi = vtknp.vtk_to_numpy(psi_array)
        hessian = vtknp.vtk_to_numpy(hessian_array)
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
        saddle_ids = np.flatnonzero(saddle_mask)

        if saddle_ids.size == 0:
            raise ValueError("No PSI saddle-point candidates found for the X point.")

        gradient = vtknp.vtk_to_numpy(gradient_array)
        if gradient.ndim != 2 or gradient.shape[1] < 2:
            raise ValueError("Unexpected PSI gradient shape while enforcing X-point stationarity criterion.")
        grad_mag = np.sqrt(np.square(gradient[:, 0]) + np.square(gradient[:, 1]))
        finite_grad_ids = saddle_ids[np.isfinite(grad_mag[saddle_ids])]
        if finite_grad_ids.size == 0:
            raise ValueError("No finite gradient magnitudes found among PSI saddle-point candidates.")

        x_point_node_id = int(finite_grad_ids[int(np.argmin(grad_mag[finite_grad_ids]))])
        psi_x = float(psi[x_point_node_id])

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
        """Interpolate poloidal magnetic field (BR, BZ) at target axis using vtkGradientFilter on PSI.

        Returns:
            None
        """
        mesh = self.vacuum_chamber_mesh
        grad_array = pyTEC3DAlgorithm._require_array(mesh.GetPointData(), "__pytec3d_psi_gradient__")
        grad_np = vtknp.vtk_to_numpy(grad_array)

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

        # 2. Interpolate dPSI/dR and dPSI/dZ at the target point using cell weights
        grad_target = weights_np @ grad_np[point_ids_np, :2]
        dpsi_dr = float(grad_target[0])
        dpsi_dz = float(grad_target[1])

        r_target = target_point[0]
        if r_target <= 0.0:
            raise ValueError(f"Invalid non-positive radial coordinate R={r_target} for target axis.")

        # 3. Convert gradient components to magnetic field: BR = -(1/R) * dPSI/dZ, BZ = (1/R) * dPSI/dR
        br = -dpsi_dz / r_target
        bz = dpsi_dr / r_target

        target = (float(br), float(bz))
        self.status["b_target_axis"] = target
        LOGGER.debug("B at target axis via vtkGradientFilter: BR=%.6g, BZ=%.6g.", target[0], target[1])

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
        mesh = self.pytec3d.model.mesh
        control_coils = self._normalized_control_coils()
        if not control_coils:
            LOGGER.debug("No control coils configured; skipping control-current update.")
            return

        cell_data = mesh.GetCellData()
        js_array = pyTEC3DAlgorithm._require_array(cell_data, "Js")
        js_np = vtknp.vtk_to_numpy(js_array)
        radial_current = float(self.status["radial_correction_current"])
        vertical_current = float(self.status["vertical_correction_current"])

        selected_mask = np.zeros(js_np.shape[0], dtype=bool)
        correction_current_all = np.zeros(js_np.shape[0], dtype=float)
        for coil in control_coils:
            component = coil["component"]
            component_array = pyTEC3DAlgorithm._require_array(cell_data, component)
            component_mask = vtknp.vtk_to_numpy(component_array) == 1
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

        real_constant_array = pyTEC3DAlgorithm._require_array(cell_data, "CONSTANT")
        real_constant_ids = vtknp.vtk_to_numpy(real_constant_array)
        constants_df = self.pytec3d.model.constants
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
            and self.pytec3d.model.time_step["time"] >= t0_halo
            and self.status["separatrix_configuration"] == "limiter"
        ):
            halo_scale = halo_delta1 * (
                1.0 + (self.pytec3d.model.time_step["time"] - t0_halo) * halo_delta2
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
        """Compute Grad-Shafranov plasma current profile and update Js using PSI gradients.

        Returns:
            None
        """
        mesh = self.vacuum_chamber_mesh

        point_data = mesh.GetPointData()
        cell_data = mesh.GetCellData()

        js_array = pyTEC3DAlgorithm._require_array(cell_data, "Js")
        js_np = vtknp.vtk_to_numpy(js_array)

        psi_point_array = pyTEC3DAlgorithm._require_array(point_data, "PSI")
        psi_point = vtknp.vtk_to_numpy(psi_point_array)

        psie_array = pyTEC3DAlgorithm._require_array(cell_data, "PSIe")
        psi_elem = vtknp.vtk_to_numpy(psie_array)

        # -------------------------------------------------------------------------
        # 1. Interpolate precomputed PSI point gradients to cell centroids
        # -------------------------------------------------------------------------
        pyTEC3DAlgorithm._require_array(point_data, "__pytec3d_psi_gradient__")

        p2c_filter = vtk.vtkPointDataToCellData()
        p2c_filter.SetInputData(mesh)
        p2c_filter.SetProcessAllArrays(False)
        p2c_filter.AddPointDataArray("__pytec3d_psi_gradient__")
        p2c_filter.Update()

        grad_cell_data = p2c_filter.GetOutput().GetCellData().GetArray("__pytec3d_psi_gradient__")
        grad_cell_np = vtknp.vtk_to_numpy(grad_cell_data)

        dpsi_dr = grad_cell_np[:, 0]
        dpsi_dz = grad_cell_np[:, 1]

        centroids_array = pyTEC3DAlgorithm._require_array(cell_data, "cell_centroids")
        centroids = vtknp.vtk_to_numpy(centroids_array)
        elem_centroid_rad = centroids[:, 0]

        # Grad-Shafranov magnetic field: B_R = -(1/R) * dPSI/dZ, B_Z = (1/R) * dPSI/dR
        # Note: If your PSI definition includes a 2*pi factor, adjust denominator to (2.0 * np.pi * elem_centroid_rad)
        bx = -dpsi_dz / elem_centroid_rad
        by = dpsi_dr / elem_centroid_rad

        # -------------------------------------------------------------------------
        # 2. Extract domain setup and parameters
        # -------------------------------------------------------------------------
        volume_array = pyTEC3DAlgorithm._require_array(cell_data, "cell_volume")
        cell_volume = vtknp.vtk_to_numpy(volume_array)

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

        # -------------------------------------------------------------------------
        # 3. Select active PLASMA elements
        # -------------------------------------------------------------------------
        Selector.select(mesh, "PLASMA", 1, type="elem", mode="s")
        selection_array = cell_data.GetArray("selection_cell")
        selected_ids = np.flatnonzero(vtknp.vtk_to_numpy(selection_array) == 1)

        # -------------------------------------------------------------------------
        # 4. Compute toroidal and poloidal current densities
        # -------------------------------------------------------------------------
        psi_segn = (psi_elem[selected_ids] - psi_magnetic_axis) / psi_den
        psi_segn = np.nan_to_num(psi_segn, nan=0.0, posinf=1.0, neginf=0.0)
        psi_segn = np.clip(psi_segn, 0.0, 1.0)
        one_minus = np.maximum(1.0 - np.power(psi_segn, alpha_mp), 0.0)
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

        # Update output current density array
        js_np[selected_ids, 0] = js_radial
        js_np[selected_ids, 1] = js_vertical
        js_np[selected_ids, 2] = js_toroidal
        js_array.Modified()

        amps_array = point_data.GetArray("AMPS")
        halo_currents_enabled = (
            self.pytec3d.model.mesh.GetMaxSpatialDimension() > 2
            and self.status["separatrix_configuration"] == "limiter"
            and float(self.input_parameters.get("halo_delta1", 0.0)) > 0.0
            and self.pytec3d.model.time_step["time"]
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
        """Check Picard convergence using max absolute PSI change on all nodes.

        Returns:
            bool: True when convergence criterion is satisfied.
        """
        mesh = self.vacuum_chamber_mesh

        point_data = mesh.GetPointData()
        psi_array = point_data.GetArray("PSI")

        current_flux_vector = vtknp.vtk_to_numpy(psi_array).copy()

        tolerance = float(self.input_parameters["tolerance_value"])
        self.status["flux_tolerance"] = tolerance

        if self._first_absolute_iteration or self._previous_flux_vector is None:
            self._previous_flux_vector = current_flux_vector.copy()
            self.status["flux_norm"] = np.nan
            self.status["converged"] = False
            return False

        if self._previous_flux_vector.size != current_flux_vector.size:
            raise ValueError("Size mismatch between current and previous flux vectors.")

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

    def toroidal_average(
        self,
        array_name: str,
        location: str,
        component: Optional[int] = None,
    ) -> None:
        """Average a mapped array component from model mesh into equilibrium mesh."""
        if location not in ("node", "elem"):
            raise ValueError("location must be 'node' or 'elem'.")

        equilibrium_mesh = self.vacuum_chamber_mesh
        model_mesh = self.pytec3d.model.mesh
        self._set_toroidal_slices()

        mapping = self._toroidal_slices_nodes if location == "node" else self._toroidal_slices_elems
        src_ids = np.asarray(mapping[:, 0], dtype=np.int64)
        dst_ids = np.asarray(mapping[:, 1:], dtype=np.int64)

        eq_data = equilibrium_mesh.GetPointData() if location == "node" else equilibrium_mesh.GetCellData()
        model_data = model_mesh.GetPointData() if location == "node" else model_mesh.GetCellData()
        eq_array = pyTEC3DAlgorithm._require_array(eq_data, array_name)
        model_array = pyTEC3DAlgorithm._require_array(model_data, array_name)

        eq_np = vtknp.vtk_to_numpy(eq_array)
        model_np = vtknp.vtk_to_numpy(model_array)

        if eq_np.ndim == 1:
            avg_values = np.mean(model_np[dst_ids], axis=1)
            eq_np[src_ids] = avg_values
            eq_array.Modified()
            return

        n_components = eq_np.shape[1]
        if component is None:
            for comp in range(n_components):
                self.toroidal_average(array_name, location, component=comp)
            return
        if component < 0 or component >= n_components:
            raise ValueError(f"Invalid component index {component} for array '{array_name}'.")

        avg_values = np.mean(model_np[dst_ids, component], axis=1)
        eq_np[src_ids, component] = avg_values
        eq_array.Modified()

    def toroidal_projection(
        self,
        array_name: str,
        location: str,
        component: Optional[int] = None,
    ) -> None:
        """Project an equilibrium-mesh array component to the mapped model mesh entities."""
        if location not in ("node", "elem"):
            raise ValueError("location must be 'node' or 'elem'.")

        equilibrium_mesh = self.vacuum_chamber_mesh
        model_mesh = self.pytec3d.model.mesh
        self._set_toroidal_slices()

        mapping = self._toroidal_slices_nodes if location == "node" else self._toroidal_slices_elems
        src_ids = np.asarray(mapping[:, 0], dtype=np.int64)
        dst_ids = np.asarray(mapping[:, 1:], dtype=np.int64)

        eq_data = equilibrium_mesh.GetPointData() if location == "node" else equilibrium_mesh.GetCellData()
        model_data = model_mesh.GetPointData() if location == "node" else model_mesh.GetCellData()
        eq_array = pyTEC3DAlgorithm._require_array(eq_data, array_name)
        model_array = pyTEC3DAlgorithm._require_array(model_data, array_name)

        eq_np = vtknp.vtk_to_numpy(eq_array)
        model_np = vtknp.vtk_to_numpy(model_array)

        flat_dst = dst_ids.reshape(-1)
        valid = flat_dst >= 0
        flat_dst = flat_dst[valid]
        if flat_dst.size == 0:
            return

        if eq_np.ndim == 1:
            projected = np.broadcast_to(eq_np[src_ids, np.newaxis], dst_ids.shape).reshape(-1)[valid]
            model_np[flat_dst] = projected
            model_array.Modified()
            return

        n_components = eq_np.shape[1]
        if component is None:
            for comp in range(n_components):
                self.toroidal_projection(array_name, location, component=comp)
            return
        if component < 0 or component >= n_components:
            raise ValueError(f"Invalid component index {component} for array '{array_name}'.")

        projected = np.broadcast_to(eq_np[src_ids, component][:, np.newaxis], dst_ids.shape).reshape(-1)[valid]
        model_np[flat_dst, component] = projected
        model_array.Modified()

    def _swap_projected_js_components_on_model_slices(self) -> None:
        """Swap projected model Js components on mapped destination element ids.

        Applies [0, 2, 1] component reordering only for 3D static meshes.
        """
        model_mesh = self.pytec3d.model.mesh
        static_mesh_is_3d = model_mesh.GetMaxSpatialDimension() > 2
        if not static_mesh_is_3d:
            return

        self._set_toroidal_slices()

        mapping = self._toroidal_slices_elems
        dst_ids = np.asarray(mapping[:, 1:], dtype=np.int64)

        model_data = model_mesh.GetCellData()
        model_array = pyTEC3DAlgorithm._require_array(model_data, "Js")

        model_np = vtknp.vtk_to_numpy(model_array)
        if model_np.ndim != 2 or model_np.shape[1] < 3:
            raise ValueError("Expected vector Js array with at least 3 components.")

        flat_dst = dst_ids.reshape(-1)
        valid = flat_dst >= 0
        if not np.any(valid):
            return

        dst = flat_dst[valid]
        model_np[dst, :] = model_np[dst][:, [0, 2, 1]]
        Selector.select(model_mesh, "ID", np.unique(dst), type="elem")
        ModelEdit.map_vector_field_to_global_components(model_mesh, "Js", YAxisCylindricalCoordinateSystem(), "elem", True)
        model_array.Modified()

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
        if self.pytec3d._world_rank() != 0:
            return

        row = self._build_diagnostics_row(converged)
        row_df = pd.DataFrame([row])
        csv_path = self._diagnostics_csv_path()
        write_header = (not os.path.exists(csv_path)) or os.path.getsize(csv_path) == 0
        row_df.to_csv(csv_path, mode="a", header=write_header, index=False)

    def calculate_poloidal_flux(self) -> None:
        """Update Ae/PSI/PSIe fields from current vector-potential solution arrays.

        Returns:
            None
        """
        mesh = self.pytec3d.model.mesh
        if mesh.GetPointData().GetArray("A") is None:
            raise ValueError("Required PointData array 'A' is missing before poloidal flux calculation.")

        ModelEdit.element_solution(mesh, "A")
        ModelEdit.poloidal_flux(
            mesh=mesh,
            coordinate_system=YAxisCylindricalCoordinateSystem(),
        )

    def _run_picard_loop(
        self,
        *,
        solve_label: str,
        equilibrium_flag: bool,
    ) -> bool:
        """Shared Picard-loop runner for equilibrium/evolution solves.

        Args:
            solve_label: Human-readable solve label for logs/errors.
            include_control_currents: Whether to run control-current update.

        Returns:
            bool: True when convergence is reached.
        """
        max_picard_iterations = max(
            1, int(self.input_parameters["max_picard_iterations"])
        )
        self.status["converged"] = False

        if self._first_absolute_iteration:
            self._set_toroidal_slices()
            self._cache_halo_boundary_graph()

        for iteration in range(max_picard_iterations):
            self.status["picard_iteration"] = iteration + 1
            LOGGER.info(
                "%s Picard iteration %d/%d.",
                solve_label,
                self.status["picard_iteration"],
                max_picard_iterations,
            )

            if not self._first_absolute_iteration:
                self._compute_psi_derivatives()
                self.find_separatrix()
                self.compute_jplasma_gradshafranov()
                if equilibrium_flag:
                    self.control_currents()
            else:
                self.check_convergence()

            self.toroidal_projection("Js", "elem")
            self._swap_projected_js_components_on_model_slices()
            self.toroidal_projection("AMPS", "node")
            self.pytec3d.solve()
            self.calculate_poloidal_flux()
            self.toroidal_average("PSI", location="node")
            self.toroidal_average("PSIe", "elem")

            converged = self.check_convergence()
            if not equilibrium_flag and self._first_absolute_iteration:
                converged = True
            self._append_iteration_diagnostics(converged)
            self._first_absolute_iteration = False
            if converged:
                self.status["converged"] = True
                LOGGER.info(
                    "%s solve converged at iteration %d (flux_norm=%s).",
                    solve_label,
                    self.status["picard_iteration"],
                    self.status["flux_norm"],
                )
                return True

        self.status["converged"] = False
        LOGGER.error(
            "%s solve did not converge after %d iterations (last flux_norm=%s).",
            solve_label,
            max_picard_iterations,
            self.status["flux_norm"],
        )
        raise RuntimeError(
            f"{solve_label} solve did not converge after "
            f"{max_picard_iterations} Picard iterations "
            f"(last flux metric: {self.status['flux_norm']})."
        )

    def equilibrium_solve(self) -> bool:
        """Run the equilibrium Picard loop until convergence or iteration limit.

        Returns:
            bool: True when convergence is reached.
        """
        LOGGER.info("Starting equilibrium solve.")
        if self.pytec3d.solver is None:
            raise RuntimeError("No solver session available. Call load_model_into_solver() before equilibrium_solve().")

        if self._first_absolute_iteration:
            self._clear_diagnostics_file()
        self.vacuum_chamber_mesh = self.pytec3d.model.mesh
        return self._run_picard_loop(
            solve_label="Equilibrium",
            equilibrium_flag=True,
        )

    def evolution_solve(self) -> bool:
        """Run the evolution Picard loop, including toroidal transfer operations.

        Returns:
            bool: True when evolution solve completes successfully.
        """
        LOGGER.info("Starting evolution solve.")
        if self.pytec3d.solver is None:
            raise RuntimeError("No solver session available. Call load_model_into_solver() before evolution_solve().")

        if self._first_absolute_iteration:
            self._clear_diagnostics_file()
        return self._run_picard_loop(
            solve_label="Evolution",
            equilibrium_flag=False,
        )

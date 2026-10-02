from typing import Any, Optional, Union, TYPE_CHECKING

import numpy as np
import vtk
import vtk.util.numpy_support as vtknp

from .logging_utils import get_logger
from .validation_utils import ensure_vtk_unstructured_grid


LOGGER = get_logger("pytec3d.selector")

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    from pytec3d.coordinate_system import CoordinateSystem


class Selector:
    """Selection and value-edit operations for VTK point/cell data arrays."""
    MODE_LABELS = {0: "s", 1: "a", 2: "r", 3: "u"}
    TYPE_LABELS = {0: "node", 1: "elem"}

    @staticmethod
    def mode_label(value: int) -> str:
        """Map numeric mode codes to selection mode characters.

        Args:
            value: Integer mode code used by ParaView filter properties.

        Returns:
            str: One of "s", "a", "r", or "u".
        """
        return Selector.MODE_LABELS.get(int(value), "s")

    @staticmethod
    def type_label(value: int) -> str:
        """Map numeric selection type to node/element labels.

        Args:
            value: Integer type code where 0=node and 1=element.

        Returns:
            str: "node" or "elem".
        """
        return Selector.TYPE_LABELS.get(int(value), "node")

    @staticmethod
    def ensure_selection_arrays(mesh: "vtk.vtkUnstructuredGrid") -> None:
        """
        Ensures that selection arrays ('selection_cell' and 'selection_point') exist in the mesh.
        If not present, initializes them to zero for all cells/points.

        Args:
            mesh: The mesh to operate on.

        Returns:
            None

        Side Effects:
            Adds selection arrays to mesh CellData and PointData if missing.
        """
        ensure_vtk_unstructured_grid(mesh)
        cell_data = mesh.GetCellData()
        if cell_data.GetArray("selection_cell") is None:
            arr = np.zeros(mesh.GetNumberOfCells(), dtype=np.uint8)
            vtk_arr = vtknp.numpy_to_vtk(
                arr, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR
            )
            vtk_arr.SetName("selection_cell")
            cell_data.AddArray(vtk_arr)
        point_data = mesh.GetPointData()
        if point_data.GetArray("selection_point") is None:
            arr = np.zeros(mesh.GetNumberOfPoints(), dtype=np.uint8)
            vtk_arr = vtknp.numpy_to_vtk(
                arr, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR
            )
            vtk_arr.SetName("selection_point")
            point_data.AddArray(vtk_arr)

    @staticmethod
    def select(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        value: Any,
        component: Union[int, str] = 0,
        type: str = "node",
        mode: str = "s",
    ) -> None:
        """
        Selects cells or points in the mesh where the specified array's component matches the given value(s).

        Supports:
            - Single value (equality)
            - Range (tuple of two values, inclusive)
            - List/set/1D NumPy array of values (membership)
            - Magnitude (for multi-component arrays, use component='magnitude')

        Selection is stored as a mask array (0/1) in the VTK grid.

        Args:
            mesh: The mesh to operate on.
            array_name: Name of the array to use for selection.
            value: Value(s) to match for selection
                (scalar, tuple range, or list/set/1D NumPy array membership).
            component: Component index or 'magnitude' for vector arrays (default 0).
            type: 'elem' or 'node' to specify selection type (default 'node').
            mode: Selection mode:
                's': select from all (replace selection)
                'a': add new selection to current selection
                'r': reselect from current selection (intersection)
                'u': remove masked ids from current selection

        Returns:
            None

        Side Effects:
            Updates selection arrays in mesh data.

        Raises:
            ValueError: If the specified array is missing in the chosen data type.
        """
        Selector.ensure_selection_arrays(mesh)
        data = mesh.GetCellData() if type == "elem" else mesh.GetPointData()
        if array_name == "ID":
            count = (
                mesh.GetNumberOfCells() if type == "elem" else mesh.GetNumberOfPoints()
            )
            arr = np.arange(count, dtype=np.int64)
        else:
            vtk_array = data.GetArray(array_name)
            if vtk_array is None:
                raise ValueError(f"Array '{array_name}' not found in {type}.")
            arr_np = vtknp.vtk_to_numpy(vtk_array)
            if component == "magnitude":
                if arr_np.ndim == 1:
                    arr = np.abs(arr_np)
                else:
                    arr = np.linalg.norm(arr_np, axis=1)
            else:
                arr = arr_np if arr_np.ndim == 1 else arr_np[:, component]
        Selector.apply_selection_mask(arr, value, data, type, mode)

    @staticmethod
    def select_pos(
        mesh: "vtk.vtkUnstructuredGrid",
        axis: int,
        value: Any,
        mode: str = "s",
        coordinate_system: Optional["CoordinateSystem"] = None,
        type: str = "node",
    ) -> None:
        """
        Selects points or cell centroids in the mesh based on their position along a given axis.

        Args:
            mesh: The mesh to operate on.
            axis: Axis index (0 for x, 1 for y, 2 for z).
            value: Value(s) to match (scalar, tuple for range, or set).
            mode: Selection mode ('s', 'a', 'r', 'u').
            coordinate_system: Optional coordinate system to transform coordinates before selection.
            type: 'elem' or 'node' to specify selection type (default 'node').

        Returns:
            None

        Side Effects:
            Updates selection arrays in mesh data.
        """
        Selector.ensure_selection_arrays(mesh)
        if type == "elem":
            data = mesh.GetCellData()
            centroids_array = data.GetArray("cell_centroids")
            coords = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
        else:
            coords = vtknp.vtk_to_numpy(mesh.GetPoints().GetData())
            data = mesh.GetPointData()
        if coordinate_system is not None:
            coords = coordinate_system.to_local(coords)
        arr = coords[:, axis]
        Selector.apply_selection_mask(arr, value, data, type, mode)

    @staticmethod
    def apply_selection_mask(
        arr: np.ndarray,
        value: Any,
        data: "vtk.vtkDataSetAttributes",
        type: str,
        mode: str,
    ) -> None:
        """
        Builds the selection mask based on value type and applies it to the selection array in the VTK grid.

        Args:
            arr: Array of values to compare.
            value: Value(s) to match (scalar, tuple range, or membership container).
            data: VTK data object (CellData or PointData).
            type: 'elem' or 'node'.
            mode: Selection mode ('s', 'a', 'r', 'u').

        Returns:
            None

        Side Effects:
            Modifies the selection array in the VTK grid.

        Raises:
            ValueError: If an invalid mode is provided or range tuple is not length 2.
        """
        if isinstance(value, tuple):
            if len(value) != 2:
                raise ValueError(
                    "Tuple for range selection must have exactly two elements."
                )
            mask = (arr >= value[0]) & (arr <= value[1])
        elif isinstance(value, (list, set)):
            mask = np.isin(arr, list(value))
        elif isinstance(value, np.ndarray):
            if value.ndim != 1:
                raise ValueError("NumPy array membership selection must be 1D.")
            mask = np.isin(arr, value)
        else:
            mask = arr == value
        mask = mask.astype(np.uint8)
        sel_array = data.GetArray(
            "selection_cell" if type == "elem" else "selection_point"
        )
        sel_np = vtknp.vtk_to_numpy(sel_array)
        if mode == "s":
            sel_np[:] = mask
        elif mode == "a":
            sel_np[:] = np.logical_or(sel_np, mask).astype(np.uint8)
        elif mode == "r":
            sel_np[:] = np.logical_and(sel_np, mask).astype(np.uint8)
        elif mode == "u":
            sel_np[:] = np.logical_and(sel_np, np.logical_not(mask)).astype(np.uint8)
        else:
            raise ValueError(f"Unknown selection mode '{mode}'.")
        sel_array.Modified()

    @staticmethod
    def apply_boolean_selection_mask(
        mask: np.ndarray,
        data: "vtk.vtkDataSetAttributes",
        type: str,
        mode: str,
    ) -> None:
        """Apply a precomputed boolean mask to node/element selection arrays.

        Args:
            mask: Boolean mask to apply.
            data: Target VTK dataset attributes.
            type: "elem" or "node" selection type.
            mode: Selection mode ('s', 'a', 'r', 'u').

        Returns:
            None
        """
        mask = np.asarray(mask, dtype=bool)
        sel_array = data.GetArray(
            "selection_cell" if type == "elem" else "selection_point"
        )
        sel_np = vtknp.vtk_to_numpy(sel_array)
        if mask.shape[0] != sel_np.shape[0]:
            raise ValueError("Selection mask has unexpected length.")
        if mode == "s":
            sel_np[:] = mask
        elif mode == "a":
            current = sel_np.astype(bool)
            np.logical_or(current, mask, out=mask)
            sel_np[:] = mask
        elif mode == "r":
            current = sel_np.astype(bool)
            np.logical_and(current, mask, out=mask)
            sel_np[:] = mask
        elif mode == "u":
            current = sel_np.astype(bool)
            np.logical_and(current, np.logical_not(mask), out=mask)
            sel_np[:] = mask
        else:
            raise ValueError(f"Unknown selection mode '{mode}'.")
        sel_array.Modified()

    @staticmethod
    def nsle(mesh: "vtk.vtkUnstructuredGrid", mode: str = "s") -> None:
        """
        Select nodes attached to the currently selected cells.

        Args:
            mesh: The mesh to operate on.
            mode: Selection mode ('s', 'a', 'r', 'u'). Defaults to 's'.

        Returns:
            None

        Raises:
            ValueError: If required selection arrays are missing.
        """
        Selector.ensure_selection_arrays(mesh)

        cell_sel_array = mesh.GetCellData().GetArray("selection_cell")
        if cell_sel_array is None:
            raise ValueError("Cell selection array 'selection_cell' not found.")
        cell_sel = vtknp.vtk_to_numpy(cell_sel_array).astype(bool)
        if not cell_sel.any():
            return

        cells = mesh.GetCells()
        offsets_array = cells.GetOffsetsArray()
        connectivity_array = cells.GetConnectivityArray()
        if offsets_array is None or connectivity_array is None:
            raise ValueError(
                "VTK connectivity arrays are required for nsle() "
                "(GetOffsetsArray/GetConnectivityArray)."
            )

        offsets = np.asarray(vtknp.vtk_to_numpy(offsets_array), dtype=np.int64)
        connectivity = np.asarray(vtknp.vtk_to_numpy(connectivity_array), dtype=np.int64)
        if offsets.size < 2 or connectivity.size == 0:
            return

        cell_lengths = np.diff(offsets)
        if np.any(cell_lengths < 0):
            raise ValueError("Invalid VTK cell connectivity offsets.")

        selected_conn = np.repeat(cell_sel, cell_lengths)
        point_mask = np.zeros(mesh.GetNumberOfPoints(), dtype=bool)
        point_mask[connectivity[selected_conn]] = True
        if not np.any(point_mask):
            return
        Selector.apply_boolean_selection_mask(
            point_mask,
            mesh.GetPointData(),
            "node",
            mode,
        )

    @staticmethod
    def esln(
        mesh: "vtk.vtkUnstructuredGrid", all_nodes: bool = True, mode: str = "s"
    ) -> None:
        """
        Select cells based on currently selected nodes.

        Args:
            mesh: The mesh to operate on.
            all_nodes: If True, select cells where all nodes are selected; otherwise any node.
            mode: Selection mode ('s', 'a', 'r', 'u'). Defaults to 's'.

        Returns:
            None

        Raises:
            ValueError: If required selection arrays are missing.
        """
        Selector.ensure_selection_arrays(mesh)

        point_sel_array = mesh.GetPointData().GetArray("selection_point")
        if point_sel_array is None:
            raise ValueError("Point selection array 'selection_point' not found.")
        point_sel = vtknp.vtk_to_numpy(point_sel_array).astype(bool)
        if not point_sel.any():
            return

        cells = mesh.GetCells()
        offsets_array = cells.GetOffsetsArray()
        connectivity_array = cells.GetConnectivityArray()
        if offsets_array is None or connectivity_array is None:
            raise ValueError(
                "VTK connectivity arrays are required for esln() "
                "(GetOffsetsArray/GetConnectivityArray)."
            )

        offsets = np.asarray(vtknp.vtk_to_numpy(offsets_array), dtype=np.int64)
        connectivity = np.asarray(vtknp.vtk_to_numpy(connectivity_array), dtype=np.int64)

        if offsets.size < 2:
            return

        n_cells = mesh.GetNumberOfCells()
        cell_lengths = np.diff(offsets)
        if np.any(cell_lengths < 0):
            raise ValueError("Invalid VTK cell connectivity offsets.")

        node_selected = point_sel[connectivity].astype(np.int8, copy=False)
        selected_counts = np.add.reduceat(node_selected, offsets[:-1])

        if all_nodes:
            cell_mask = (cell_lengths > 0) & (selected_counts == cell_lengths)
        else:
            cell_mask = selected_counts > 0

        if not np.any(cell_mask):
            return
        Selector.apply_boolean_selection_mask(
            cell_mask,
            mesh.GetCellData(),
            "elem",
            mode,
        )

    @staticmethod
    def select_all(mesh: "vtk.vtkUnstructuredGrid", type: Optional[str] = None) -> None:
        """
        Selects all cells or points in the mesh.

        Args:
            mesh: The mesh to operate on.
            type: 'elem' or 'node'. If None, selects both.

        Returns:
            None

        Side Effects:
            Sets all selection array values to 1 (selected).
        """
        Selector.ensure_selection_arrays(mesh)
        if type is None or type not in ["elem", "node"]:
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 1
            mesh.GetCellData().GetArray("selection_cell").Modified()
            point_sel = vtknp.vtk_to_numpy(
                mesh.GetPointData().GetArray("selection_point")
            )
            point_sel[:] = 1
            mesh.GetPointData().GetArray("selection_point").Modified()
        elif type == "elem":
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 1
            mesh.GetCellData().GetArray("selection_cell").Modified()
        elif type == "node":
            point_sel = vtknp.vtk_to_numpy(
                mesh.GetPointData().GetArray("selection_point")
            )
            point_sel[:] = 1
            mesh.GetPointData().GetArray("selection_point").Modified()

    @staticmethod
    def unselect_all(
        mesh: "vtk.vtkUnstructuredGrid", type: Optional[str] = None
    ) -> None:
        """
        Unselects all cells or points in the mesh.

        Args:
            mesh: The mesh to operate on.
            type: 'elem' or 'node'. If None, unselects both.

        Returns:
            None

        Side Effects:
            Sets all selection array values to 0 (unselected).
        """
        Selector.ensure_selection_arrays(mesh)
        if type is None or type not in ["elem", "node"]:
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 0
            mesh.GetCellData().GetArray("selection_cell").Modified()
            point_sel = vtknp.vtk_to_numpy(
                mesh.GetPointData().GetArray("selection_point")
            )
            point_sel[:] = 0
            mesh.GetPointData().GetArray("selection_point").Modified()
        elif type == "elem":
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 0
            mesh.GetCellData().GetArray("selection_cell").Modified()
        elif type == "node":
            point_sel = vtknp.vtk_to_numpy(
                mesh.GetPointData().GetArray("selection_point")
            )
            point_sel[:] = 0
            mesh.GetPointData().GetArray("selection_point").Modified()

    @staticmethod
    def _extract_selected_external_edges_2d(
        mesh: "vtk.vtkUnstructuredGrid",
        selected_cell_ids: np.ndarray,
        temp_id_name: str = "__pytec3d_original_point_ids__",
    ) -> "vtk.vtkPolyData":
        """Extract external boundary edges of selected 2D cells as polydata.

        Args:
            mesh: Source 2D mesh.
            selected_cell_ids: Selected cell ids.
            temp_id_name: Temporary point-id array name.

        Returns:
            vtk.vtkPolyData: Boundary edges as polydata.
        """
        ensure_vtk_unstructured_grid(mesh)
        mesh_point_data = mesh.GetPointData()
        created_temp_ids = False
        if mesh_point_data.GetArray(temp_id_name) is None:
            tmp_ids = np.arange(mesh.GetNumberOfPoints(), dtype=np.int64)
            vtk_tmp_ids = vtknp.numpy_to_vtkIdTypeArray(tmp_ids, deep=True)
            vtk_tmp_ids.SetName(temp_id_name)
            mesh_point_data.AddArray(vtk_tmp_ids)
            mesh_point_data.Modified()
            created_temp_ids = True

        try:
            id_array = vtknp.numpy_to_vtkIdTypeArray(selected_cell_ids, deep=True)
            selection_node = vtk.vtkSelectionNode()
            selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            selection_node.SetSelectionList(id_array)
            selection = vtk.vtkSelection()
            selection.AddNode(selection_node)

            extract_selection = vtk.vtkExtractSelection()
            extract_selection.SetInputData(0, mesh)
            extract_selection.SetInputData(1, selection)
            extract_selection.Update()

            to_poly = vtk.vtkGeometryFilter()
            to_poly.SetInputData(extract_selection.GetOutput())
            to_poly.Update()

            feature_edges = vtk.vtkFeatureEdges()
            feature_edges.SetInputData(to_poly.GetOutput())
            feature_edges.BoundaryEdgesOn()
            feature_edges.FeatureEdgesOff()
            feature_edges.ManifoldEdgesOff()
            feature_edges.NonManifoldEdgesOff()
            feature_edges.Update()

            edges = vtk.vtkPolyData()
            edges.DeepCopy(feature_edges.GetOutput())
            return edges
        finally:
            if created_temp_ids:
                mesh_point_data.RemoveArray(temp_id_name)
                mesh_point_data.Modified()

    @staticmethod
    def nsel_ext(mesh: "vtk.vtkUnstructuredGrid", mode: str = "s") -> None:
        """
        Selects all external nodes (points) of the currently selected cells.

        This method finds the surface nodes of the selected cells and updates the
        'selection_point' array in PointData to select those nodes.

        Args:
            mesh: The mesh to operate on.

        Returns:
            None

        Side Effects:
            Updates 'selection_point' array in mesh PointData to select external nodes.

        Raises:
            ValueError: If required selection arrays are missing.
        """
        Selector.ensure_selection_arrays(mesh)
        cell_sel_array = mesh.GetCellData().GetArray("selection_cell")
        if cell_sel_array is None:
            raise ValueError(
                "Cell selection array 'selection_cell' not found. Run ensure_selection_arrays first."
            )
        cell_sel = vtknp.vtk_to_numpy(cell_sel_array)
        selected_cell_ids = np.flatnonzero(cell_sel == 1).astype(np.int64)
        if selected_cell_ids.size == 0:
            LOGGER.info("No cells are selected.")
            return

        temp_id_name = "__pytec3d_original_point_ids__"

        if mesh.GetMaxSpatialDimension() <= 2:
            edges = Selector._extract_selected_external_edges_2d(
                mesh,
                selected_cell_ids,
                temp_id_name=temp_id_name,
            )
            point_data = edges.GetPointData()
            vtk_ids_array = point_data.GetArray(temp_id_name)
            if vtk_ids_array is None:
                vtk_ids_array = point_data.GetArray("vtkOriginalPointIds")
            if vtk_ids_array is None:
                LOGGER.info("No external nodes found (%s missing).", temp_id_name)
                return
            point_ids = vtknp.vtk_to_numpy(vtk_ids_array)
            if point_ids.size == 0:
                LOGGER.info("No external nodes found (feature edges returned empty).")
                return
        else:
            id_array = vtknp.numpy_to_vtkIdTypeArray(selected_cell_ids, deep=True)
            selection_node = vtk.vtkSelectionNode()
            selection_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            selection_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            selection_node.SetSelectionList(id_array)
            selection = vtk.vtkSelection()
            selection.AddNode(selection_node)
            extract_selection = vtk.vtkExtractSelection()
            extract_selection.SetInputData(0, mesh)
            extract_selection.SetInputData(1, selection)
            mesh_point_data = mesh.GetPointData()
            created_temp_ids = False
            if mesh_point_data.GetArray(temp_id_name) is None:
                tmp_ids = np.arange(mesh.GetNumberOfPoints(), dtype=np.int64)
                vtk_tmp_ids = vtknp.numpy_to_vtkIdTypeArray(tmp_ids, deep=True)
                vtk_tmp_ids.SetName(temp_id_name)
                mesh_point_data.AddArray(vtk_tmp_ids)
                mesh_point_data.Modified()
                created_temp_ids = True
            try:
                extract_selection.Update()
                selected_mesh = extract_selection.GetOutput()
                extract_surface = vtk.vtkGeometryFilter()
                extract_surface.SetInputData(selected_mesh)
                extract_surface.Update()
                surface = extract_surface.GetOutput()
                point_data = surface.GetPointData()
                vtk_ids_array = point_data.GetArray(temp_id_name)
                if vtk_ids_array is None:
                    vtk_ids_array = point_data.GetArray("vtkOriginalPointIds")
                if vtk_ids_array is None:
                    LOGGER.info("No external nodes found (%s missing).", temp_id_name)
                    return
                point_ids = vtknp.vtk_to_numpy(vtk_ids_array)
            finally:
                if created_temp_ids:
                    mesh_point_data.RemoveArray(temp_id_name)
                    mesh_point_data.Modified()

        Selector.select(
            mesh,
            "ID",
            set(point_ids.tolist()),
            type="node",
            mode=mode,
        )

    @staticmethod
    def get_selected_external_node_adjacency(
        mesh: "vtk.vtkUnstructuredGrid",
    ) -> tuple[np.ndarray, dict[int, np.ndarray]]:
        """
        Return external node ids and edge adjacency for selected 2D cells.

        The boundary is the 1D polydata produced by the same GeometryFilter +
        FeatureEdges pipeline used by nsel_ext().

        Args:
            mesh: Source 2D mesh.

        Returns:
            tuple[np.ndarray, dict[int, np.ndarray]]: Boundary node ids and adjacency map.
        """
        Selector.ensure_selection_arrays(mesh)
        if mesh.GetMaxSpatialDimension() > 2:
            raise ValueError(
                "get_selected_external_node_adjacency() is intended for 2D meshes only."
            )

        cell_sel_array = mesh.GetCellData().GetArray("selection_cell")
        selected_cell_ids = np.flatnonzero(vtknp.vtk_to_numpy(cell_sel_array) == 1)
        if selected_cell_ids.size == 0:
            return np.empty((0,), dtype=np.int64), {}

        temp_id_name = "__pytec3d_original_point_ids__"
        edges = Selector._extract_selected_external_edges_2d(
            mesh,
            selected_cell_ids,
            temp_id_name=temp_id_name,
        )
        original_point_ids = np.empty((0,), dtype=np.int64)
        point_data = edges.GetPointData()
        vtk_ids_array = point_data.GetArray(temp_id_name)
        if vtk_ids_array is None:
            vtk_ids_array = point_data.GetArray("vtkOriginalPointIds")
        if vtk_ids_array is None:
            LOGGER.info("No external nodes found (%s missing).", temp_id_name)
            boundary_node_ids = np.empty((0,), dtype=np.int64)
            edge_connectivity = None
            edge_offsets = None
        else:
            original_point_ids = np.asarray(
                vtknp.vtk_to_numpy(vtk_ids_array),
                dtype=np.int64,
            )
            boundary_node_ids = np.unique(original_point_ids)
            lines = edges.GetLines()
            edge_connectivity = (
                lines.GetConnectivityArray() if lines is not None else None
            )
            edge_offsets = lines.GetOffsetsArray() if lines is not None else None

        if boundary_node_ids.size == 0 or edge_connectivity is None or edge_offsets is None:
            return boundary_node_ids, {}

        adjacency_sets: dict[int, set[int]] = {
            int(node_id): set() for node_id in boundary_node_ids
        }
        offsets = np.asarray(vtknp.vtk_to_numpy(edge_offsets), dtype=np.int64)
        connectivity = np.asarray(vtknp.vtk_to_numpy(edge_connectivity), dtype=np.int64)
        for start, stop in zip(offsets[:-1], offsets[1:]):
            path = original_point_ids[connectivity[int(start) : int(stop)]]
            for left, right in zip(path[:-1], path[1:]):
                if left == right:
                    continue
                left = int(left)
                right = int(right)
                adjacency_sets.setdefault(left, set()).add(right)
                adjacency_sets.setdefault(right, set()).add(left)

        adjacency = {
            node_id: np.fromiter(sorted(neighbors), dtype=np.int64, count=len(neighbors))
            for node_id, neighbors in adjacency_sets.items()
        }
        return boundary_node_ids, adjacency

    @staticmethod
    def set_value_on_selection(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        value: Any,
        component: int = 0,
        type: str = "node",
    ) -> None:
        """
        Sets the value of the specified array for all selected cells or points in the mesh.

        Args:
            mesh: The mesh to operate on.
            array_name: Name of the array to update.
            value: Value to set for selected elements.
            component: Component index for multi-component arrays (default 0).
            type: 'elem' or 'node' to specify selection type (default 'node').

        Returns:
            None

        Side Effects:
            Updates array values for selected elements in mesh data.

        Raises:
            ValueError: If the target array is missing.
        """
        Selector.ensure_selection_arrays(mesh)
        data = mesh.GetCellData() if type == "elem" else mesh.GetPointData()
        vtk_array = data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {type}.")
        arr_np = vtknp.vtk_to_numpy(vtk_array)
        sel_array = data.GetArray(
            "selection_cell" if type == "elem" else "selection_point"
        )
        sel_np = vtknp.vtk_to_numpy(sel_array)
        selected_ids = np.where(sel_np == 1)[0]
        if arr_np.ndim == 1:
            arr_np[selected_ids] = value
        else:
            arr_np[selected_ids, component] = value
        vtk_array.Modified()

    @staticmethod
    def get_selection_values(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        component: Union[int, str] = 0,
        type: str = "node",
    ) -> np.ndarray:
        """
        Retrieves array values for all currently selected cells or points.

        Args:
            mesh: The mesh to operate on.
            array_name: Name of the array to retrieve values from (or "ID").
            component: Component index or 'magnitude' for vector arrays (default 0).
            type: 'elem' or 'node' to specify selection type (default 'node').

        Returns:
            np.ndarray: 1D NumPy array containing array values for selected elements.

        Raises:
            ValueError: If the specified array is missing in cell or point data.
        """
        Selector.ensure_selection_arrays(mesh)
        data = mesh.GetCellData() if type == "elem" else mesh.GetPointData()

        sel_array_name = "selection_cell" if type == "elem" else "selection_point"
        sel_array = data.GetArray(sel_array_name)
        sel_np = vtknp.vtk_to_numpy(sel_array)
        selected_ids = np.flatnonzero(sel_np == 1)

        if array_name == "ID":
            count = (
                mesh.GetNumberOfCells() if type == "elem" else mesh.GetNumberOfPoints()
            )
            arr = np.arange(count, dtype=np.int64)
        else:
            vtk_array = data.GetArray(array_name)
            if vtk_array is None:
                raise ValueError(f"Array '{array_name}' not found in {type}.")
            arr_np = vtknp.vtk_to_numpy(vtk_array)
            if component == "magnitude":
                if arr_np.ndim == 1:
                    arr = np.abs(arr_np)
                else:
                    arr = np.linalg.norm(arr_np, axis=1)
            else:
                arr = arr_np if arr_np.ndim == 1 else arr_np[:, component]

        return arr[selected_ids]

    @staticmethod
    def inv(mesh: "vtk.vtkUnstructuredGrid", type: Optional[str] = None) -> None:
        """Inverts the current selection status for cells, points, or both.

        Args:
            mesh: The mesh to operate on.
            type: 'elem' or 'node'. If None, inverts both cell and point selections.

        Returns:
            None

        Side Effects:
            Toggles values (0 -> 1, 1 -> 0) in the selection array(s) and calls Modified().
        """
        Selector.ensure_selection_arrays(mesh)

        if type is None or type not in ["elem", "node"]:
            cell_arr = mesh.GetCellData().GetArray("selection_cell")
            cell_sel = vtknp.vtk_to_numpy(cell_arr)
            cell_sel[:] = 1 - cell_sel
            cell_arr.Modified()

            point_arr = mesh.GetPointData().GetArray("selection_point")
            point_sel = vtknp.vtk_to_numpy(point_arr)
            point_sel[:] = 1 - point_sel
            point_arr.Modified()
        elif type == "elem":
            cell_arr = mesh.GetCellData().GetArray("selection_cell")
            cell_sel = vtknp.vtk_to_numpy(cell_arr)
            cell_sel[:] = 1 - cell_sel
            cell_arr.Modified()
        elif type == "node":
            point_arr = mesh.GetPointData().GetArray("selection_point")
            point_sel = vtknp.vtk_to_numpy(point_arr)
            point_sel[:] = 1 - point_sel
            point_arr.Modified()

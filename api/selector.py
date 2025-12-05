from typing import Any, Optional, Union, TYPE_CHECKING

import numpy as np
import vtk
import vtk.util.numpy_support as vtknp

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import vtkmodules.vtkCommonDataModel as vtk_dm  # type: ignore[import-not-found]


class Selector:
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
        cell_data = mesh.GetCellData()
        if cell_data.GetArray("selection_cell") is None:
            arr = np.zeros(mesh.GetNumberOfCells(), dtype=np.uint8)
            vtk_arr = vtknp.numpy_to_vtk(arr, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
            vtk_arr.SetName("selection_cell")
            cell_data.AddArray(vtk_arr)
        point_data = mesh.GetPointData()
        if point_data.GetArray("selection_point") is None:
            arr = np.zeros(mesh.GetNumberOfPoints(), dtype=np.uint8)
            vtk_arr = vtknp.numpy_to_vtk(arr, deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
            vtk_arr.SetName("selection_point")
            point_data.AddArray(vtk_arr)

    @staticmethod
    def select(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        value: Any,
        component: Union[int, str] = 0,
        type: str = "point_data",
        mode: str = "s",
    ) -> None:
        """
        Selects cells or points in the mesh where the specified array's component matches the given value(s).

        Supports:
            - Single value (equality)
            - Range (tuple of two values, inclusive)
            - Set of values (membership)
            - Magnitude (for multi-component arrays, use component='magnitude')

        Selection is stored as a mask array (0/1) in the VTK grid.

        Args:
            mesh: The mesh to operate on.
            array_name: Name of the array to use for selection.
            value: Value(s) to match for selection (scalar, tuple, or set).
            component: Component index or 'magnitude' for vector arrays (default 0).
            type: 'cell_data' or 'point_data' to specify selection type (default 'point_data').
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
        data = mesh.GetCellData() if type == "cell_data" else mesh.GetPointData()
        if array_name == "ID":
            count = mesh.GetNumberOfCells() if type == "cell_data" else mesh.GetNumberOfPoints()
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
    ) -> None:
        """
        Selects points in the mesh based on their position along a given axis.

        Args:
            mesh: The mesh to operate on.
            axis: Axis index (0 for x, 1 for y, 2 for z).
            value: Value(s) to match (scalar, tuple for range, or set).
            mode: Selection mode ('s', 'a', 'r', 'u').

        Returns:
            None

        Side Effects:
            Updates selection arrays in mesh PointData.
        """
        Selector.ensure_selection_arrays(mesh)
        coords = vtknp.vtk_to_numpy(mesh.GetPoints().GetData())
        arr = coords[:, axis]
        data = mesh.GetPointData()
        Selector.apply_selection_mask(arr, value, data, "point_data", mode)

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
            value: Value(s) to match (scalar, tuple, or set).
            data: VTK data object (CellData or PointData).
            type: 'cell_data' or 'point_data'.
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
                raise ValueError("Tuple for range selection must have exactly two elements.")
            mask = (arr >= value[0]) & (arr <= value[1])
        elif isinstance(value, set):
            mask = np.isin(arr, list(value))
        else:
            mask = arr == value
        mask = mask.astype(np.uint8)
        sel_array = data.GetArray("selection_cell" if type == "cell_data" else "selection_point")
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

        selected_ids = np.nonzero(cell_sel)[0]
        if selected_ids.size == 0:
            return

        # Build a selection of the chosen cell ids
        id_array = vtknp.numpy_to_vtkIdTypeArray(selected_ids, deep=True)
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
        selected_mesh = extract_selection.GetOutput()

        original_ids_array = selected_mesh.GetPointData().GetArray("vtkOriginalPointIds")
        if original_ids_array is None:
            # Fallback: if ids are missing, nothing to select
            return
        original_point_ids = vtknp.vtk_to_numpy(original_ids_array)

        point_sel_array = mesh.GetPointData().GetArray("selection_point")
        if point_sel_array is None:
            raise ValueError("Point selection array 'selection_point' not found.")
        point_sel = vtknp.vtk_to_numpy(point_sel_array)
        point_sel[:] = 0
        point_sel[original_point_ids] = 1
        point_sel_array.Modified()

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

        selected_cells = []
        for cid in range(mesh.GetNumberOfCells()):
            cell = mesh.GetCell(cid)
            pts = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
            if not pts:
                continue
            node_mask = point_sel[np.asarray(pts, dtype=int)]
            if all_nodes:
                if node_mask.all():
                    selected_cells.append(cid)
            else:
                if node_mask.any():
                    selected_cells.append(cid)

        if not selected_cells:
            return

        Selector.select(
            mesh,
            "ID",
            set(selected_cells),
            type="cell_data",
            mode=mode,
        )

    @staticmethod
    def select_all(mesh: "vtk.vtkUnstructuredGrid", type: Optional[str] = None) -> None:
        """
        Selects all cells or points in the mesh.

        Args:
            mesh: The mesh to operate on.
            type: 'cell_data' or 'point_data'. If None, selects both.

        Returns:
            None

        Side Effects:
            Sets all selection array values to 1 (selected).
        """
        Selector.ensure_selection_arrays(mesh)
        if type is None or type not in ["cell_data", "point_data"]:
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 1
            mesh.GetCellData().GetArray("selection_cell").Modified()
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 1
            mesh.GetPointData().GetArray("selection_point").Modified()
        elif type == "cell_data":
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 1
            mesh.GetCellData().GetArray("selection_cell").Modified()
        elif type == "point_data":
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 1
            mesh.GetPointData().GetArray("selection_point").Modified()

    @staticmethod
    def unselect_all(mesh: "vtk.vtkUnstructuredGrid", type: Optional[str] = None) -> None:
        """
        Unselects all cells or points in the mesh.

        Args:
            mesh: The mesh to operate on.
            type: 'cell_data' or 'point_data'. If None, unselects both.

        Returns:
            None

        Side Effects:
            Sets all selection array values to 0 (unselected).
        """
        Selector.ensure_selection_arrays(mesh)
        if type is None or type not in ["cell_data", "point_data"]:
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 0
            mesh.GetCellData().GetArray("selection_cell").Modified()
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 0
            mesh.GetPointData().GetArray("selection_point").Modified()
        elif type == "cell_data":
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 0
            mesh.GetCellData().GetArray("selection_cell").Modified()
        elif type == "point_data":
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 0
            mesh.GetPointData().GetArray("selection_point").Modified()

    @staticmethod
    def nsel_ext(mesh: "vtk.vtkUnstructuredGrid") -> None:
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
            raise ValueError("Cell selection array 'selection_cell' not found. Run ensure_selection_arrays first.")
        cell_sel = vtknp.vtk_to_numpy(cell_sel_array)
        selected_cell_ids = np.where(cell_sel == 1)[0]
        if selected_cell_ids.size == 0:
            print("No cells are selected.")
            return
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
        selected_mesh = extract_selection.GetOutput()
        extract_surface = vtk.vtkGeometryFilter()
        extract_surface.SetInputData(selected_mesh)
        extract_surface.Update()
        surface = extract_surface.GetOutput()
        vtk_ids_array = surface.GetPointData().GetArray("vtkOriginalPointIds")
        if vtk_ids_array is None:
            print("No external nodes found (vtkOriginalPointIds missing).")
            return
        surface_point_ids = vtknp.vtk_to_numpy(vtk_ids_array)
        point_sel_array = mesh.GetPointData().GetArray("selection_point")
        if point_sel_array is None:
            raise ValueError("Point selection array 'selection_point' not found. Run ensure_selection_arrays first.")
        point_sel = vtknp.vtk_to_numpy(point_sel_array)
        point_sel[:] = 0
        point_sel[surface_point_ids] = 1
        point_sel_array.Modified()

    @staticmethod
    def set_value_on_selection(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        value: Any,
        component: int = 0,
        type: str = "point_data",
    ) -> None:
        """
        Sets the value of the specified array for all selected cells or points in the mesh.

        Args:
            mesh: The mesh to operate on.
            array_name: Name of the array to update.
            value: Value to set for selected elements.
            component: Component index for multi-component arrays (default 0).
            type: 'cell_data' or 'point_data' to specify selection type (default 'point_data').

        Returns:
            None

        Side Effects:
            Updates array values for selected elements in mesh data.

        Raises:
            ValueError: If the target array is missing.
        """
        Selector.ensure_selection_arrays(mesh)
        data = mesh.GetCellData() if type == "cell_data" else mesh.GetPointData()
        vtk_array = data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {type}.")
        arr_np = vtknp.vtk_to_numpy(vtk_array)
        sel_array = data.GetArray("selection_cell" if type == "cell_data" else "selection_point")
        sel_np = vtknp.vtk_to_numpy(sel_array)
        selected_ids = np.where(sel_np == 1)[0]
        if arr_np.ndim == 1:
            arr_np[selected_ids] = value
        else:
            arr_np[selected_ids, component] = value
        vtk_array.Modified()

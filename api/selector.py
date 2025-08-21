import numpy as np
import vtk
import vtk.util.numpy_support as vtknp

class Selector:
    @staticmethod
    def ensure_selection_arrays(mesh):
        """
        Ensures that selection arrays ('selection_cell' and 'selection_point') exist in the mesh.
        If not present, initializes them to zero for all cells/points.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.

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
    def select(mesh, array_name, value, component=0, cell_or_point="point_data", mode="s"):
        """
        Selects cells or points in the mesh where the specified array's component matches the given value(s).

        Supports:
            - Single value (equality)
            - Range (tuple of two values, inclusive)
            - Set of values (membership)
            - Magnitude (for multi-component arrays, use component='magnitude')

        Selection is stored as a mask array (0/1) in the VTK grid.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.
            array_name (str): Name of the array to use for selection.
            value: Value(s) to match for selection (scalar, tuple, or set).
            component (int or str): Component index or 'magnitude' for vector arrays (default 0).
            cell_or_point (str): 'cell_data' or 'point_data' to specify selection type (default 'point_data').
            mode (str): Selection mode:
                's': select from all (replace selection)
                'a': add new selection to current selection
                'r': reselect from current selection (intersection)
                'u': remove masked ids from current selection

        Side Effects:
            Updates selection arrays in mesh data.
        """
        Selector.ensure_selection_arrays(mesh)
        data = mesh.GetCellData() if cell_or_point == "cell_data" else mesh.GetPointData()
        vtk_array = data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {cell_or_point}.")
        arr_np = vtknp.vtk_to_numpy(vtk_array)
        if component == "magnitude":
            if arr_np.ndim == 1:
                arr = np.abs(arr_np)
            else:
                arr = np.linalg.norm(arr_np, axis=1)
        else:
            arr = arr_np if arr_np.ndim == 1 else arr_np[:, component]
        Selector.apply_selection_mask(arr, value, data, cell_or_point, mode)

    @staticmethod
    def select_pos(mesh, axis, value, mode="s"):
        """
        Selects points in the mesh based on their position along a given axis.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.
            axis (int): Axis index (0 for x, 1 for y, 2 for z).
            value: Value(s) to match (scalar, tuple for range, or set).
            mode (str): Selection mode ('s', 'a', 'r', 'u').

        Side Effects:
            Updates selection arrays in mesh PointData.
        """
        Selector.ensure_selection_arrays(mesh)
        coords = vtknp.vtk_to_numpy(mesh.GetPoints().GetData())
        arr = coords[:, axis]
        data = mesh.GetPointData()
        Selector.apply_selection_mask(arr, value, data, "point_data", mode)

    @staticmethod
    def apply_selection_mask(arr, value, data, cell_or_point, mode):
        """
        Builds the selection mask based on value type and applies it to the selection array in the VTK grid.

        Parameters:
            arr (np.ndarray): Array of values to compare.
            value: Value(s) to match (scalar, tuple, or set).
            data (vtk.vtkDataSetAttributes): VTK data object (CellData or PointData).
            cell_or_point (str): 'cell_data' or 'point_data'.
            mode (str): Selection mode ('s', 'a', 'r', 'u').

        Side Effects:
            Modifies the selection array in the VTK grid.
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
        sel_array = data.GetArray("selection_cell" if cell_or_point == "cell_data" else "selection_point")
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
    def select_all(mesh, cell_or_point=None):
        """
        Selects all cells or points in the mesh.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.
            cell_or_point (str, optional): 'cell_data' or 'point_data'. If None, selects both.

        Side Effects:
            Sets all selection array values to 1 (selected).
        """
        Selector.ensure_selection_arrays(mesh)
        if cell_or_point is None or cell_or_point not in ["cell_data", "point_data"]:
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 1
            mesh.GetCellData().GetArray("selection_cell").Modified()
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 1
            mesh.GetPointData().GetArray("selection_point").Modified()
        elif cell_or_point == "cell_data":
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 1
            mesh.GetCellData().GetArray("selection_cell").Modified()
        elif cell_or_point == "point_data":
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 1
            mesh.GetPointData().GetArray("selection_point").Modified()

    @staticmethod
    def unselect_all(mesh, cell_or_point=None):
        """
        Unselects all cells or points in the mesh.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.
            cell_or_point (str, optional): 'cell_data' or 'point_data'. If None, unselects both.

        Side Effects:
            Sets all selection array values to 0 (unselected).
        """
        Selector.ensure_selection_arrays(mesh)
        if cell_or_point is None or cell_or_point not in ["cell_data", "point_data"]:
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 0
            mesh.GetCellData().GetArray("selection_cell").Modified()
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 0
            mesh.GetPointData().GetArray("selection_point").Modified()
        elif cell_or_point == "cell_data":
            cell_sel = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("selection_cell"))
            cell_sel[:] = 0
            mesh.GetCellData().GetArray("selection_cell").Modified()
        elif cell_or_point == "point_data":
            point_sel = vtknp.vtk_to_numpy(mesh.GetPointData().GetArray("selection_point"))
            point_sel[:] = 0
            mesh.GetPointData().GetArray("selection_point").Modified()

    @staticmethod
    def nsel_ext(mesh):
        """
        Selects all external nodes (points) of the currently selected cells.

        This method finds the surface nodes of the selected cells and updates the
        'selection_point' array in PointData to select those nodes.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.

        Side Effects:
            Updates 'selection_point' array in mesh PointData to select external nodes.
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
    def set_value_on_selection(mesh, array_name, value, component=0, cell_or_point="point_data"):
        """
        Sets the value of the specified array for all selected cells or points in the mesh.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh to operate on.
            array_name (str): Name of the array to update.
            value: Value to set for selected elements.
            component (int): Component index for multi-component arrays (default 0).
            cell_or_point (str): 'cell_data' or 'point_data' to specify selection type (default 'point_data').

        Side Effects:
            Updates array values for selected elements in mesh data.
        """
        Selector.ensure_selection_arrays(mesh)
        data = mesh.GetCellData() if cell_or_point == "cell_data" else mesh.GetPointData()
        vtk_array = data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {cell_or_point}.")
        arr_np = vtknp.vtk_to_numpy(vtk_array)
        sel_array = data.GetArray("selection_cell" if cell_or_point == "cell_data" else "selection_point")
        sel_np = vtknp.vtk_to_numpy(sel_array)
        selected_ids = np.where(sel_np == 1)[0]
        if arr_np.ndim == 1:
            arr_np[selected_ids] = value
        else:
            arr_np[selected_ids, component] = value
        vtk_array.Modified()

    # def nsle(self):
    #     """
    #     Selects the nodes attached to currently selected cells.
    #     Updates the 'selection_point' array in PointData.
    #     """
    #     cell_sel_array = self.mesh.GetCellData().GetArray("selection_cell")
    #     if cell_sel_array is None:
    #         raise ValueError(
    #             "Cell selection array 'selection_cell' not found. Run _ensure_selection_arrays first."
    #         )
    #     cell_sel = vtknp.vtk_to_numpy(cell_sel_array)
    #     selected_cell_ids = np.where(cell_sel == 1)[0]
    #     if len(selected_cell_ids) == 0:
    #         print("No cells are selected.")
    #         return

    #     cell_array = vtknp.vtk_to_numpy(self.mesh.GetCells().GetData())
    #     cell_offsets = vtknp.vtk_to_numpy(self.mesh.GetCellLocationsArray())

    #     # Efficiently collect all node IDs attached to selected cells
    #     offsets = cell_offsets[selected_cell_ids]
    #     npts = cell_array[offsets]
    #     # Build a list of all node indices for selected cells
    #     all_node_ids = [
    #         cell_array[offset + 1 : offset + 1 + npt]
    #         for offset, npt in zip(offsets, npts)
    #     ]
    #     # Flatten and get unique node IDs
    #     attached_nodes = np.unique(np.concatenate(all_node_ids))

    #     point_sel_array = self.mesh.GetPointData().GetArray("selection_point")
    #     if point_sel_array is None:
    #         raise ValueError(
    #             "Point selection array 'selection_point' not found. Run _ensure_selection_arrays first."
    #         )
    #     point_sel = vtknp.vtk_to_numpy(point_sel_array)
    #     point_sel[:] = 0
    #     point_sel[attached_nodes] = 1
    #     point_sel_array.Modified()

    # def esln(self, all_nodes=False):
    #     """
    #     Selects cells based on node selection:
    #     - If all_nodes=False: selects cells with at least one selected node.
    #     - If all_nodes=True: selects cells where all nodes are selected.
    #     Updates the 'selection_cell' array in CellData.
    #     """
    #     point_sel_array = self.mesh.GetPointData().GetArray("selection_point")
    #     if point_sel_array is None:
    #         raise ValueError(
    #             "Point selection array 'selection_point' not found. Run _ensure_selection_arrays first."
    #         )
    #     point_sel = vtknp.vtk_to_numpy(point_sel_array)
    #     selected_points = np.where(point_sel == 1)[0]
    #     if selected_points.size == 0:
    #         print("No nodes are selected.")
    #         return

    #     # Get connectivity arrays as NumPy
    #     cell_array = vtknp.vtk_to_numpy(self.mesh.GetCells().GetData())
    #     cell_offsets = vtknp.vtk_to_numpy(self.mesh.GetCellLocationsArray())
    #     num_cells = self.mesh.GetNumberOfCells()
    #     cell_sel_array = self.mesh.GetCellData().GetArray("selection_cell")
    #     if cell_sel_array is None:
    #         raise ValueError(
    #             "Cell selection array 'selection_cell' not found. Run _ensure_selection_arrays first."
    #         )
    #     cell_sel = vtknp.vtk_to_numpy(cell_sel_array)
    #     cell_sel[:] = 0

    #     # Efficient cell classification
    #     mask_any = np.zeros(num_cells, dtype=bool)
    #     mask_all = np.zeros(num_cells, dtype=bool)
    #     selected_set = set(selected_points)

    #     for cid in range(num_cells):
    #         start = cell_offsets[cid]
    #         npts = cell_array[start]
    #         pts = cell_array[start + 1 : start + 1 + npts]
    #         sel_mask = np.isin(pts, selected_points)
    #         mask_any[cid] = sel_mask.any()
    #         mask_all[cid] = sel_mask.all()

    #     if all_nodes:
    #         cell_sel[:] = mask_all.astype(np.uint8)
    #     else:
    #         cell_sel[:] = mask_any.astype(np.uint8)
    #     cell_sel_array.Modified()
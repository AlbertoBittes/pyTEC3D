from typing import Any, Dict, Optional, TYPE_CHECKING

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp

from constants import VTK_ARRAYS_TO_INIT

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import vtkmodules.vtkCommonDataModel as vtk_dm  # type: ignore[import-not-found]

class VTKGridHelper:
    @staticmethod
    def add_array_to_grid(
        mesh: "vtk.vtkUnstructuredGrid",
        name: str,
        location: str = "point",
        components: int = 1,
        fill: Any = 0,
        dtype: Optional[int] = vtk.VTK_DOUBLE,
        values: Optional[np.ndarray] = None,
    ) -> "vtk.vtkDataArray":
        """
        Adds a VTK array to the mesh at the specified location.

        Args:
            mesh: The VTK unstructured grid to add the array to.
            name: Name of the array.
            location: "cell" or "point" (default "point").
            components: Number of components (default 1).
            fill: Default value to fill the array if values are not provided (default 0).
            dtype: VTK data type (default vtk.VTK_DOUBLE). If None, dtype validation is skipped.
            values: Explicit values to store in the array. Must match expected shape.

        Returns:
            vtkDataArray: The created or existing VTK array.

        Raises:
            ValueError: If location is invalid, shapes/components mismatch, or dtype differs.
        """
        if location not in ["cell", "point"]:
            raise ValueError(f"Unknown location '{location}'. Use 'cell' or 'point'.")

        data_obj = mesh.GetCellData() if location == "cell" else mesh.GetPointData()
        num_entries = mesh.GetNumberOfCells() if location == "cell" else mesh.GetNumberOfPoints()
        num_components = max(1, components)
        expected_shape = (num_entries,) if num_components == 1 else (num_entries, num_components)

        values_np = None
        if values is not None:
            values_np = np.asarray(values)
            if values_np.shape == (num_entries, 1) and num_components == 1:
                values_np = values_np.reshape(num_entries)
            if values_np.shape != expected_shape:
                raise ValueError(
                    f"Values for array '{name}' must have shape {expected_shape}, got {values_np.shape}."
                )

        vtk_array = data_obj.GetArray(name)
        if vtk_array is None:
            arr_to_store = values_np if values_np is not None else np.full(expected_shape, fill)
            vtk_array = vtknp.numpy_to_vtk(arr_to_store, deep=True, array_type=dtype)
            vtk_array.SetName(name)
            vtk_array.SetNumberOfComponents(num_components)
            data_obj.AddArray(vtk_array)
        else:
            if vtk_array.GetNumberOfComponents() != num_components:
                raise ValueError(
                    f"Existing array '{name}' has {vtk_array.GetNumberOfComponents()} components, expected {num_components}."
                )
            if dtype is not None and vtk_array.GetDataType() != dtype:
                raise ValueError(
                    f"Existing array '{name}' has VTK type {vtk_array.GetDataType()}, expected {dtype}."
                )

        if values_np is not None:
            vtk_np = vtknp.vtk_to_numpy(vtk_array)
            vtk_np[:] = values_np
            vtk_array.Modified()

        return vtk_array

    @staticmethod
    def initialize_vtk_arrays(
        mesh: "vtk.vtkUnstructuredGrid",
        arrays_to_init: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        """
        Initializes required VTK arrays in cell/point data as specified in arrays_to_init.

        Args:
            mesh: The VTK unstructured grid to initialize arrays for.
            arrays_to_init: Dictionary specifying arrays to initialize; defaults to VTK_ARRAYS_TO_INIT.

        Returns:
            None
        """
        if arrays_to_init is None:
            arrays_to_init = VTK_ARRAYS_TO_INIT

        for name, info in arrays_to_init.items():
            VTKGridHelper.add_array_to_grid(
                mesh,
                name=name,
                location=info["location"],
                components=info["components"],
                fill=info.get("fill", 0),
                dtype=info.get("dtype", vtk.VTK_DOUBLE),
            )

    @staticmethod
    def initialize_vtk_array(
        mesh: "vtk.vtkUnstructuredGrid",
        name: str,
        location: str = "point",
        components: int = 1,
        fill: Any = 0,
        dtype: int = vtk.VTK_DOUBLE,
    ) -> None:
        """
        Initializes a single VTK array in cell or point data.

        Args:
            mesh: The VTK unstructured grid to initialize the array for.
            name: Name of the array.
            location: "cell" or "point" (default "point").
            components: Number of components (default 1).
            fill: Default value to fill the array (default 0).
            dtype: VTK data type (default vtk.VTK_DOUBLE).
        """
        arrays_to_init = {
            name: {
                "location": location,
                "components": components,
                "fill": fill,
                "dtype": dtype,
            }
        }
        VTKGridHelper.initialize_vtk_arrays(mesh, arrays_to_init)

    @staticmethod
    def reset_vtk_arrays_to_default(
        mesh: "vtk.vtkUnstructuredGrid",
        arrays_to_init: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> None:
        """
        Resets arrays specified in arrays_to_init to their default 'fill' values if present.

        Args:
            mesh: The VTK unstructured grid whose arrays will be reset.
            arrays_to_init: Dictionary specifying arrays to reset; defaults to VTK_ARRAYS_TO_INIT.
        """
        if arrays_to_init is None:
            from constants import VTK_ARRAYS_TO_INIT
            arrays_to_init = VTK_ARRAYS_TO_INIT

        cell_data = mesh.GetCellData()
        point_data = mesh.GetPointData()

        for name, info in arrays_to_init.items():
            data_obj = cell_data if info["location"] == "cell" else point_data
            vtk_array = data_obj.GetArray(name)
            if vtk_array is not None and "fill" in info:
                num = (
                    mesh.GetNumberOfCells()
                    if info["location"] == "cell"
                    else mesh.GetNumberOfPoints()
                )
                shape = (num, info["components"]) if info["components"] > 1 else num
                arr = np.full(shape, info["fill"])
                # Overwrite the VTK array with the default values
                vtk_np = vtknp.vtk_to_numpy(vtk_array)
                vtk_np[:] = arr
                vtk_array.Modified()

    @staticmethod
    def safe_remove_array(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        location: Optional[str] = None,
    ) -> None:
        """
        Safely removes an array from CellData or PointData of the mesh if it exists.

        Args:
            mesh (vtk.vtkUnstructuredGrid): The mesh from which to remove the array.
            array_name: The name of the array to remove.
            location: "cell", "point", or None. If None, tries to guess location.

        Side Effects:
            Removes the array from the detected or specified data location if present.
        """
        # Try to guess location if not provided
        if location is None:
            cell_data = mesh.GetCellData()
            point_data = mesh.GetPointData()
            if cell_data.HasArray(array_name):
                cell_data.RemoveArray(array_name)
                return
            elif point_data.HasArray(array_name):
                point_data.RemoveArray(array_name)
                return
            else:
                print(f"Array '{array_name}' not found in CellData or PointData.")
                return
        elif location == "cell":
            data = mesh.GetCellData()
        elif location == "point":
            data = mesh.GetPointData()
        else:
            raise ValueError(f"Unknown location '{location}'. Use 'cell', 'point', or None.")
        if data.HasArray(array_name):
            data.RemoveArray(array_name)

    @staticmethod
    def add_boundary_array(mesh: "vtk.vtkUnstructuredGrid") -> str:
        """Add a new boundary_id_N array to point data with default zeros."""
        point_data = mesh.GetPointData()
        existing = []
        for i in range(point_data.GetNumberOfArrays()):
            name = point_data.GetArrayName(i)
            if name and name.startswith("boundary_id_"):
                try:
                    existing.append(int(name[len("boundary_id_"):]))
                except ValueError:
                    continue

        next_index = 0
        if existing:
            next_index = max(existing) + 1

        array_name = f"boundary_id_{next_index}"
        VTKGridHelper.add_array_to_grid(
            mesh,
            name=array_name,
            location="point",
            components=1,
            fill=0,
            dtype=vtk.VTK_INT,
        )
        return array_name


    @staticmethod
    def validate_dataframe_with_vtk_data(
        vtk_data: "vtk.vtkDataSetAttributes",
        df: pd.DataFrame,
        array_name: str,
    ) -> None:
        """
        Validates that the first column of a DataFrame has unique values,
        all values have a correspondence in the specified VTK data array,
        and there are no undefined values in the VTK data array. Useful
        to check 1-to-1 correspondence between IDs.

        Args:
            vtk_data: The VTK cell or point data to validate against.
            df (pd.DataFrame): The DataFrame to validate.
            array_name: The name of the VTK data array to check.

        Raises:
            ValueError: If the DataFrame is empty, the first column is not unique, the VTK data array is missing,
                or the checks fail.
        """
        if df.empty:
            print("The DataFrame is empty.")
            return
        first_column = df.iloc[:, 0]
        if not first_column.is_unique:
            raise ValueError(
                "The first column of the DataFrame must have unique values."
            )
        vtk_array = vtk_data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"The '{array_name}' array is missing in the VTK data.")
        vtk_array_np = vtknp.vtk_to_numpy(vtk_array)
        if not set(first_column).issubset(set(vtk_array_np)):
            raise ValueError(
                f"Some values in the first column of the DataFrame do not have a correspondence in the '{array_name}' array."
            )
        undefined_values = set(vtk_array_np) - set(first_column)
        if undefined_values:
            raise ValueError(
                f"The VTK data contains undefined values: {undefined_values}."
            )
        print(f"Validation passed for '{array_name}' array and DataFrame.")

    @staticmethod
    def update_vtk_grid(data: "vtk.vtkDataSetAttributes", df) -> None:
        """
        Updates the VTK grid's cell/point data arrays using the DataFrame resulting from `replace_ansys_with_vtk_ids`.
        The function assigns values to the appropriate components at the corresponding VTK IDs.

        Args:
            data: VTK DataSetAttributes (CellData or PointData)
            df (pd.DataFrame): The DataFrame containing columns 'vtk_id', 'label', 'component', and 'value'.

        Returns:
            None

        Raises:
            ValueError: If any expected label array is missing in the VTK data.
        """

        loads_bc_labels = df["label"].unique()
        for load_bc_label in loads_bc_labels:
            load_bc_df = df[df["label"] == load_bc_label]
            vtk_array = data.GetArray(load_bc_label)
            if vtk_array is None:
                raise ValueError(f"Array '{load_bc_label}' not found in the VTK file.")
            vtk_array_np = vtknp.vtk_to_numpy(vtk_array)
            if vtk_array_np.ndim == 1:
                vtk_array_np[:] = np.nan
            else:
                vtk_array_np[:] = np.full(vtk_array_np.shape, np.nan)
            if vtk_array_np.ndim == 1:
                for component in load_bc_df["component"].unique():
                    component_mask = load_bc_df["component"] == component
                    vtk_ids = load_bc_df.loc[component_mask, "vtk_id"].to_numpy()
                    values = load_bc_df.loc[component_mask, "value"].to_numpy()
                    vtk_array_np[vtk_ids] = values
            else:
                for component in load_bc_df["component"].unique():
                    component_mask = load_bc_df["component"] == component
                    vtk_ids = load_bc_df.loc[component_mask, "vtk_id"].to_numpy()
                    values = load_bc_df.loc[component_mask, "value"].to_numpy()
                    vtk_array_np[vtk_ids, component] = values
            vtk_array.Modified()

    @staticmethod
    def compute_cell_centroids(mesh: "vtk.vtkUnstructuredGrid") -> np.ndarray:
        """
        Computes the geometric centroids for all cells in the given VTK unstructured grid.

        Args:
            mesh: The mesh whose cell centroids are to be computed.

        Returns:
            np.ndarray: Array of shape (num_cells, 3) containing the centroid coordinates (x, y, z) for each cell.

        Notes:
            Uses vtkCellCenters to compute centroids. The output array is ordered by cell index.
        """
        cellCenters = vtk.vtkCellCenters()
        cellCenters.SetInputData(mesh)
        cellCenters.Update()
        vtk_points = cellCenters.GetOutput().GetPoints().GetData()
        np_centers = vtknp.vtk_to_numpy(vtk_points)
        return np_centers

    @staticmethod
    def store_cell_centroids(
        mesh: "vtk.vtkUnstructuredGrid", array_name: str = "cell_centroids"
    ) -> None:
        """
        Computes cell centroids and stores them in a VTK array in CellData.

        Args:
            mesh: The mesh whose cell centroids are to be computed and stored.
            array_name: Name of the array to store centroids in CellData (default "cell_centroids").

        Side Effects:
            Adds/overwrites the centroid array on mesh.GetCellData().
        """
        np_centers = VTKGridHelper.compute_cell_centroids(mesh)
        VTKGridHelper.add_array_to_grid(
            mesh,
            name=array_name,
            location="cell",
            components=3,
            dtype=vtk.VTK_DOUBLE,
            values=np_centers,
        )

    @staticmethod
    def compute_cell_volumes(mesh: "vtk.vtkUnstructuredGrid") -> np.ndarray:
        """
        Computes the volume (for 3D cells) or area (for 2D cells) for each cell in the mesh.

        Args:
            mesh: The mesh whose cell volumes/areas are to be computed.

        Returns:
            np.ndarray: Array of shape (num_cells,) containing the volume (or area) for each cell.

        Notes:
            Uses ComputeVolume or ComputeArea if available on the cell object. If not available, sets value to 0.0.
            The output array is ordered by cell index.
        """
        num_cells = mesh.GetNumberOfCells()
        volumes = np.zeros(num_cells)
        for i in range(num_cells):
            cell = mesh.GetCell(i)
            if hasattr(cell, "ComputeVolume"):
                try:
                    volumes[i] = cell.ComputeVolume()
                except Exception:
                    volumes[i] = 0.0
            elif hasattr(cell, "ComputeArea"):
                try:
                    volumes[i] = cell.ComputeArea()
                except Exception:
                    volumes[i] = 0.0
            else:
                volumes[i] = 0.0
        return volumes

    @staticmethod
    def store_cell_volumes(
        mesh: "vtk.vtkUnstructuredGrid", array_name: str = "cell_volume"
    ) -> None:
        """
        Computes cell volumes/areas and stores them in a VTK array in CellData.

        Args:
            mesh: The mesh whose cell volumes/areas are to be computed and stored.
            array_name: Name of the array to store volumes/areas in CellData (default "cell_volume").

        Side Effects:
            Adds/overwrites the volume/area array on mesh.GetCellData().
        """
        volumes = VTKGridHelper.compute_cell_volumes(mesh)
        VTKGridHelper.add_array_to_grid(
            mesh,
            name=array_name,
            location="cell",
            components=1,
            dtype=vtk.VTK_DOUBLE,
            values=volumes,
        )

    @staticmethod
    def add_component(
        mesh: "vtk.vtkUnstructuredGrid", component_name: str, type: str = "cell"
    ) -> "vtk.vtkDataArray":
        """
        Adds a component array marking the current selection.

        Args:
            mesh: The mesh to operate on.
            component_name: Name of the component array to create.
            type: "cell" to use CellData selection, "point" for PointData selection.

        Returns:
            vtkDataArray: The created component array.
        """
        if type not in ["cell", "point"]:
            raise ValueError("type must be 'cell' or 'point'.")

        # Avoid circular import at module load time
        from selector import Selector

        Selector.ensure_selection_arrays(mesh)
        data_obj = mesh.GetCellData() if type == "cell" else mesh.GetPointData()
        selection_array = data_obj.GetArray("selection_cell" if type == "cell" else "selection_point")
        if selection_array is None:
            raise ValueError("Selection array not found; ensure selection arrays are initialized.")

        selection_mask = vtknp.vtk_to_numpy(selection_array)
        new_array = VTKGridHelper.add_array_to_grid(
            mesh,
            name=component_name,
            location=type,
            components=1,
            fill=0,
            dtype=vtk.VTK_UNSIGNED_CHAR,
        )
        target_np = vtknp.vtk_to_numpy(new_array)
        target_np[:] = selection_mask
        new_array.Modified()
        return new_array

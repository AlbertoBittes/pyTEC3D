from collections import defaultdict
import os
from os import PathLike
from typing import Any, Dict, List, Optional, TYPE_CHECKING, Tuple, Union
import warnings

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp

from constants import VTK_ARRAYS_TO_INIT

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import vtkmodules.vtkCommonDataModel as vtk_dm  # type: ignore[import-not-found]

class VTKGridHelper:
    @staticmethod
    def read_mesh(
        vtk_file_path: Union[str, PathLike[str]],
        reset_arrays: bool = True,
    ) -> "vtk.vtkUnstructuredGrid":
        """
        Reads a VTK unstructured grid from a file in XML format (.vtu or .pvtu).

        Args:
            vtk_file_path: Path to the VTK file in XML format.
            reset_arrays: When True, reset initialized arrays to defaults after loading.

        Returns:
            vtk.vtkUnstructuredGrid: The loaded mesh.
        """
        path = os.fspath(vtk_file_path)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Mesh file '{path}' not found.")

        _, ext = os.path.splitext(path)
        ext = ext.lower()
        if ext == ".pvtu":
            reader = vtk.vtkXMLPUnstructuredGridReader()
        else:
            reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(path)
        reader.Update()
        mesh = reader.GetOutput()

        if reset_arrays:
            # Avoid circular import at module load time
            from selector import Selector
            from tree_gauge_helper import TreeGaugeHelper

            VTKGridHelper.initialize_vtk_arrays(mesh)
            Selector.ensure_selection_arrays(mesh)
            VTKGridHelper.store_cell_centroids(mesh)
            TreeGaugeHelper.add_midside_nodes_mask_array(mesh)
            VTKGridHelper.reset_vtk_arrays_to_default(mesh)

        print(f"VTK unstructured grid loaded from '{path}'.")
        print(f"Number of points: {mesh.GetNumberOfPoints()}")
        print(f"Number of cells: {mesh.GetNumberOfCells()}")
        print(f"Number of point data arrays: {mesh.GetPointData().GetNumberOfArrays()}")
        print(f"Number of cell data arrays: {mesh.GetCellData().GetNumberOfArrays()}")
        print(f"Number of field data arrays: {mesh.GetFieldData().GetNumberOfArrays()}")

        return mesh

    @staticmethod
    def add_array_to_grid(
        mesh: "vtk.vtkUnstructuredGrid",
        name: str,
        location: str = "point",
        components: int = 1,
        fill: Any = 0,
        dtype: Optional[int] = vtk.VTK_DOUBLE,
        values: Optional[np.ndarray] = None,
        replace: bool = False,
    ) -> Union[np.ndarray, "vtk.vtkDataArray"]:
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
            vtkDataArray or np.ndarray: The created VTK array, or an existing array's NumPy view.

        Raises:
            ValueError: If location is invalid, shapes/components mismatch, or dtype differs.
        """
        if location not in ["cell", "point"]:
            raise ValueError(f"Unknown location '{location}'. Use 'cell' or 'point'.")

        data_obj = mesh.GetCellData() if location == "cell" else mesh.GetPointData()
        num_entries = mesh.GetNumberOfCells() if location == "cell" else mesh.GetNumberOfPoints()
        num_components = max(1, components)
        expected_shape = (num_entries,) if num_components == 1 else (num_entries, num_components)

        vtk_array = data_obj.GetArray(name)
        if vtk_array is not None and not replace:
            return vtknp.vtk_to_numpy(vtk_array)

        if replace and vtk_array is not None:
            data_obj.RemoveArray(name)
            vtk_array = None

        values_np = None
        if not replace and values is not None:
            values_np = np.asarray(values)
            if values_np.shape == (num_entries, 1) and num_components == 1:
                values_np = values_np.reshape(num_entries)
            if values_np.shape != expected_shape:
                raise ValueError(
                    f"Values for array '{name}' must have shape {expected_shape}, got {values_np.shape}."
                )

        arr_to_store = values_np
        if replace:
            arr_to_store = np.full(expected_shape, fill)
        elif arr_to_store is None:
            arr_to_store = np.full(expected_shape, fill)

        vtk_array = vtknp.numpy_to_vtk(arr_to_store, deep=True, array_type=dtype)
        vtk_array.SetName(name)
        vtk_array.SetNumberOfComponents(num_components)
        data_obj.AddArray(vtk_array)
        vtk_array.Modified()

        return vtk_array

    @staticmethod
    def save_mesh(
        mesh: "vtk.vtkUnstructuredGrid",
        output_file: Union[str, PathLike[str]],
        binary: bool = True,
    ) -> None:
        """
        Save an unstructured grid to a .vtu file.

        Args:
            mesh: The VTK unstructured grid to save.
            output_file: Destination file path.
            binary: When True, write in binary mode; otherwise ASCII.

        Returns:
            None
        """
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(os.fspath(output_file))
        writer.SetInputData(mesh)
        if binary:
            writer.SetDataModeToBinary()
        else:
            writer.SetDataModeToAscii()
        writer.Write()
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
    def build_toroidal_projection_map(
        input_mesh: "vtk.vtkUnstructuredGrid",
        target_mesh: "vtk.vtkUnstructuredGrid",
        location: str = "point",
        tolerance: float = 1e-6,
        origin: Optional[np.ndarray] = None,
        axis: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Build a toroidal projection map from a 2D input mesh to a 3D target mesh.

        Each input point/cell centroid is associated with all target points/cell centroids that share
        the same radial distance from the rotation axis and the same axial coordinate within
        the given tolerance.
        The resulting map is returned as a 2D NumPy array where the first column holds input IDs and
        the remaining columns hold the corresponding target IDs (padded with NaNs when needed).
        If the input mesh has a selection array for the chosen location, only selected IDs are used.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            target_mesh: Target unstructured grid (typically 3D).
            location: "point" to map point IDs or "cell" to map cell IDs.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            origin: A point on the rotation axis. Defaults to (0, 0, 0).
            axis: Rotation axis direction. Defaults to (0, 1, 0).

        Returns:
            np.ndarray: 2D array with shape (num_input, 1 + max_matches).

        Raises:
            ValueError: If inputs are invalid or the input mesh is not 2D.
        """
        from selector import Selector

        if location not in ["point", "cell"]:
            raise ValueError("location must be 'point' or 'cell'.")
        if tolerance <= 0:
            raise ValueError("tolerance must be positive.")

        Selector.ensure_selection_arrays(input_mesh)
        axis_vec = np.array(
            axis if axis is not None else [0.0, 1.0, 0.0], dtype=float
        ).reshape(-1)
        if axis_vec.size != 3:
            raise ValueError("axis must be a 3-element vector.")
        axis_norm = np.linalg.norm(axis_vec)
        if axis_norm <= 0.0:
            raise ValueError("axis must be a non-zero vector.")
        axis_unit = axis_vec / axis_norm
        origin_vec = np.array(
            origin if origin is not None else [0.0, 0.0, 0.0], dtype=float
        ).reshape(-1)
        if origin_vec.size != 3:
            raise ValueError("origin must be a 3-element vector.")

        input_points = input_mesh.GetPoints()
        if input_points is None:
            return np.empty((0, 1), dtype=float)
        input_points_np = vtknp.vtk_to_numpy(input_points.GetData())

        if input_mesh.GetNumberOfCells() > 0 and input_mesh.GetMaxSpatialDimension() > 2:
            raise ValueError("Input mesh must be 2D.")

        def get_positions(mesh: "vtk.vtkUnstructuredGrid") -> np.ndarray:
            if location == "point":
                points = mesh.GetPoints()
                if points is None:
                    return np.empty((0, 3))
                return vtknp.vtk_to_numpy(points.GetData())
            return VTKGridHelper.compute_cell_centroids(mesh)

        def get_input_positions_and_ids(
            mesh: "vtk.vtkUnstructuredGrid",
        ) -> Tuple[np.ndarray, np.ndarray]:
            positions = get_positions(mesh)
            if location == "point":
                count = mesh.GetNumberOfPoints()
                data = mesh.GetPointData()
                selection_name = "selection_point"
            else:
                count = mesh.GetNumberOfCells()
                data = mesh.GetCellData()
                selection_name = "selection_cell"
            ids = np.arange(count, dtype=int)
            selection_array = data.GetArray(selection_name) if data is not None else None
            if selection_array is None:
                return positions, ids
            selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
            if selection_mask.size != ids.size:
                raise ValueError(
                    f"Selection array '{selection_name}' has unexpected length."
                )
            if np.any(selection_mask):
                return positions[selection_mask], ids[selection_mask]
            return positions[:0], ids[:0]

        input_pos, input_ids = get_input_positions_and_ids(input_mesh)
        target_pos = get_positions(target_mesh)

        if input_pos.size == 0:
            return np.empty((0, 1), dtype=float)

        input_vec = input_pos - origin_vec
        target_vec = target_pos - origin_vec
        input_axial = np.dot(input_vec, axis_unit)
        target_axial = np.dot(target_vec, axis_unit)
        input_radial_vec = input_vec - np.outer(input_axial, axis_unit)
        target_radial_vec = target_vec - np.outer(target_axial, axis_unit)
        input_radius = np.linalg.norm(input_radial_vec, axis=1)
        target_radius = np.linalg.norm(target_radial_vec, axis=1)

        scale = 1.0 / tolerance
        def quantize(r: float, axial: float) -> tuple:
            return (int(round(r * scale)), int(round(axial * scale)))

        target_index: Dict[tuple, List[int]] = defaultdict(list)
        for idx, (r, axial) in enumerate(zip(target_radius, target_axial)):
            target_index[quantize(r, axial)].append(idx)

        num_input = input_radius.shape[0]
        matches_per_input: List[List[int]] = []
        for input_id in range(num_input):
            matches = target_index.get(
                quantize(input_radius[input_id], input_axial[input_id]), []
            )
            matches_per_input.append(sorted(matches))

        max_matches = max((len(m) for m in matches_per_input), default=0)
        mapping = np.full((num_input, 1 + max_matches), np.nan)
        mapping[:, 0] = input_ids.astype(int)

        counts = []
        for row_idx, matches in enumerate(matches_per_input):
            counts.append(len(matches))
            if matches:
                mapping[row_idx, 1 : 1 + len(matches)] = matches

        if max_matches > 0 and any(count != max_matches for count in counts):
            warnings.warn(
                "Toroidal projection map has varying numbers of target IDs per source ID; "
                "rows were padded with NaNs. This should not happen for an extruded mesh.",
                RuntimeWarning,
            )

        return mapping

    @staticmethod
    def toroidal_projection(
        input_mesh: "vtk.vtkUnstructuredGrid",
        input_array_name: str,
        input_component: Optional[int],
        target_mesh: "vtk.vtkUnstructuredGrid",
        target_array_name: str,
        target_component: Optional[int] = None,
        projection_map: Optional[np.ndarray] = None,
        location: str = "point",
        tolerance: float = 1e-6,
        origin: Optional[np.ndarray] = None,
        axis: Optional[np.ndarray] = None,
    ) -> None:
        """
        Project scalar or vector fields from a 2D input mesh onto a 3D target mesh.

        Modes:
            - input_component is None: treat the input as a 3-component vector field and rotate
              it toroidally to each matched target position.
            - input_component and target_component are set: treat as a scalar value and copy it
              into the target component for each matched position.
            - input_component is set and target_component is None: treat the selected component
              as a toroidal (theta) component and write it as a vector aligned with the local
              toroidal unit direction around the given axis.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            input_array_name: Name of the input array on the input mesh.
            input_component: Component index to read from the input array, or None for full vector.
            target_mesh: Target unstructured grid (typically 3D).
            target_array_name: Name of the target array to write on the target mesh.
            projection_map: Optional precomputed toroidal projection map.
            location: "point" to map point data or "cell" to map cell data.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            origin: A point on the rotation axis. Defaults to (0, 0, 0).
            axis: Rotation axis direction. Defaults to (0, 1, 0).
            target_component: Component index to write for scalar copy mode.

        Returns:
            None
        """
        if location not in ["point", "cell"]:
            raise ValueError("location must be 'point' or 'cell'.")
        if input_component is None and target_component is not None:
            raise ValueError("target_component cannot be set when input_component is None.")

        input_data = input_mesh.GetPointData() if location == "point" else input_mesh.GetCellData()
        input_array = input_data.GetArray(input_array_name) if input_data is not None else None
        if input_array is None:
            raise ValueError(f"Array '{input_array_name}' not found in input {location} data.")

        input_np = vtknp.vtk_to_numpy(input_array)
        input_components = input_array.GetNumberOfComponents()
        is_vector_mode = input_component is None
        is_scalar_copy_mode = input_component is not None and target_component is not None
        is_theta_mode = input_component is not None and target_component is None

        if is_vector_mode:
            if input_components != 3:
                raise ValueError("input_array_name must refer to a 3-component array.")
            input_vectors = input_np
        else:
            if input_components == 1:
                if input_component != 0:
                    raise ValueError("input_component must be 0 for scalar arrays.")
                input_values = input_np
            else:
                if input_component < 0 or input_component >= input_np.shape[1]:
                    raise ValueError("input_component out of range for input array.")
                input_values = input_np[:, input_component]

        if projection_map is None:
            projection_map = VTKGridHelper.build_toroidal_projection_map(
                input_mesh,
                target_mesh,
                location=location,
                tolerance=tolerance,
                origin=origin,
                axis=axis,
            )

        if projection_map.size == 0:
            return

        target_data = target_mesh.GetPointData() if location == "point" else target_mesh.GetCellData()
        if target_data is None:
            raise ValueError(f"Target mesh has no {location} data.")
        target_array = target_data.GetArray(target_array_name)
        if target_array is None:
            target_components = 3 if (is_vector_mode or is_theta_mode) else max(1, target_component + 1)
            target_array = VTKGridHelper.add_array_to_grid(
                target_mesh,
                name=target_array_name,
                location=location,
                components=target_components,
                fill=0.0,
                dtype=vtk.VTK_DOUBLE,
            )
        if is_vector_mode or is_theta_mode:
            if target_array.GetNumberOfComponents() != 3:
                raise ValueError("target_array_name must refer to a 3-component array.")
        else:
            if target_component < 0:
                raise ValueError("target_component must be non-negative.")
            if target_array.GetNumberOfComponents() == 1 and target_component != 0:
                raise ValueError("target_component must be 0 for scalar target arrays.")
            if target_component >= target_array.GetNumberOfComponents():
                raise ValueError("target_component out of range for target array.")

        target_np = vtknp.vtk_to_numpy(target_array)

        axis_unit = None
        origin_vec = None
        target_radial_unit = None
        target_toroidal_unit = None
        input_radial_unit = None
        input_toroidal_unit = None
        if is_vector_mode or is_theta_mode:
            axis_vec = np.array(
                axis if axis is not None else [0.0, 1.0, 0.0], dtype=float
            ).reshape(-1)
            if axis_vec.size != 3:
                raise ValueError("axis must be a 3-element vector.")
            axis_norm = np.linalg.norm(axis_vec)
            if axis_norm <= 0.0:
                raise ValueError("axis must be a non-zero vector.")
            axis_unit = axis_vec / axis_norm
            origin_vec = np.array(
                origin if origin is not None else [0.0, 0.0, 0.0], dtype=float
            ).reshape(-1)
            if origin_vec.size != 3:
                raise ValueError("origin must be a 3-element vector.")

            if location == "point":
                points = target_mesh.GetPoints()
                target_pos = np.empty((0, 3)) if points is None else vtknp.vtk_to_numpy(points.GetData())
                input_points = input_mesh.GetPoints()
                input_pos = np.empty((0, 3)) if input_points is None else vtknp.vtk_to_numpy(input_points.GetData())
            else:
                target_pos = VTKGridHelper.compute_cell_centroids(target_mesh)
                input_pos = VTKGridHelper.compute_cell_centroids(input_mesh)

            if target_pos.size == 0:
                return

            target_vec = target_pos - origin_vec
            target_axial = np.dot(target_vec, axis_unit)
            target_radial_vec = target_vec - np.outer(target_axial, axis_unit)
            target_radial_norm = np.linalg.norm(target_radial_vec, axis=1)
            target_radial_unit = np.zeros_like(target_radial_vec)
            target_nonzero = target_radial_norm > 0.0
            target_radial_unit[target_nonzero] = target_radial_vec[target_nonzero] / target_radial_norm[target_nonzero, None]
            target_toroidal_unit = np.cross(axis_unit, target_radial_unit)

            if is_vector_mode:
                if input_pos.size == 0:
                    return
                input_vec = input_pos - origin_vec
                input_axial = np.dot(input_vec, axis_unit)
                input_radial_vec = input_vec - np.outer(input_axial, axis_unit)
                input_radial_norm = np.linalg.norm(input_radial_vec, axis=1)
                input_radial_unit = np.zeros_like(input_radial_vec)
                input_nonzero = input_radial_norm > 0.0
                input_radial_unit[input_nonzero] = input_radial_vec[input_nonzero] / input_radial_norm[input_nonzero, None]
                input_toroidal_unit = np.cross(axis_unit, input_radial_unit)

        for row in projection_map:
            input_id = row[0]
            if np.isnan(input_id):
                continue
            input_id_int = int(input_id)
            if is_vector_mode:
                if input_id_int < 0 or input_id_int >= input_vectors.shape[0]:
                    continue
                input_vec = input_vectors[input_id_int]
                vr = np.dot(input_vec, input_radial_unit[input_id_int])
                vt = np.dot(input_vec, input_toroidal_unit[input_id_int])
                va = np.dot(input_vec, axis_unit)
            else:
                if input_id_int < 0 or input_id_int >= input_values.shape[0]:
                    continue
                value = input_values[input_id_int]
            for target_id in row[1:]:
                if np.isnan(target_id):
                    continue
                target_id_int = int(target_id)
                if target_id_int < 0 or target_id_int >= target_np.shape[0]:
                    continue
                if is_vector_mode:
                    target_np[target_id_int, :] = (
                        vr * target_radial_unit[target_id_int]
                        + va * axis_unit
                        + vt * target_toroidal_unit[target_id_int]
                    )
                elif is_scalar_copy_mode:
                    if target_np.ndim == 1:
                        target_np[target_id_int] = value
                    else:
                        target_np[target_id_int, target_component] = value
                else:
                    target_np[target_id_int, :] = value * target_toroidal_unit[target_id_int]

        target_array.Modified()

    @staticmethod
    def toroidal_projection_inverse(
        input_mesh: "vtk.vtkUnstructuredGrid",
        input_array_name: str,
        input_component: Optional[int],
        target_mesh: "vtk.vtkUnstructuredGrid",
        target_array_name: str,
        target_component: Optional[int] = None,
        projection_map: Optional[np.ndarray] = None,
        location: str = "point",
        tolerance: float = 1e-6,
        origin: Optional[np.ndarray] = None,
        axis: Optional[np.ndarray] = None,
    ) -> None:
        """
        Average toroidal fields from a 3D input mesh onto a 2D target mesh.

        Modes:
            - input_component is None and target_component is None: treat both arrays as 3-component
              vectors. Decompose input into radial/toroidal/axial components, average across
              toroidal positions, and reconstruct on the target mesh.
            - input_component and target_component are set: treat as scalar values and average
              across toroidal positions into the target component.

        Args:
            input_mesh: 3D unstructured grid providing the vector field.
            input_array_name: Name of the input array on the input mesh.
            input_component: Component index to read from the input array, or None for full vector.
            target_mesh: 2D unstructured grid used to generate a surface of revolution.
            target_array_name: Name of the target array on the target mesh.
            target_component: Component index to write for scalar targets (default None/0).
            projection_map: Optional precomputed toroidal projection map.
            location: "point" to map point data or "cell" to map cell data.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            origin: A point on the rotation axis. Defaults to (0, 0, 0).
            axis: Rotation axis direction. Defaults to (0, 1, 0).

        Returns:
            None
        """
        if location not in ["point", "cell"]:
            raise ValueError("location must be 'point' or 'cell'.")
        if (input_component is None) != (target_component is None):
            raise ValueError("input_component and target_component must both be None or both set.")

        input_data = input_mesh.GetPointData() if location == "point" else input_mesh.GetCellData()
        input_array = input_data.GetArray(input_array_name) if input_data is not None else None
        if input_array is None:
            raise ValueError(f"Array '{input_array_name}' not found in input {location} data.")
        input_np = vtknp.vtk_to_numpy(input_array)
        input_components = input_array.GetNumberOfComponents()
        if input_components not in [1, 3]:
            raise ValueError("input_array_name must refer to a scalar or 3-component array.")

        if projection_map is None:
            projection_map = VTKGridHelper.build_toroidal_projection_map(
                target_mesh,
                input_mesh,
                location=location,
                tolerance=tolerance,
                origin=origin,
                axis=axis,
            )

        if projection_map.size == 0:
            return

        target_data = target_mesh.GetPointData() if location == "point" else target_mesh.GetCellData()
        if target_data is None:
            raise ValueError(f"Target mesh has no {location} data.")
        target_array = target_data.GetArray(target_array_name)
        if target_array is None:
            target_components = 1 if input_components == 1 else 3
            target_array = VTKGridHelper.add_array_to_grid(
                target_mesh,
                name=target_array_name,
                location=location,
                components=target_components,
                fill=np.nan,
                dtype=vtk.VTK_DOUBLE,
            )
        target_components = target_array.GetNumberOfComponents()
        if input_component is None:
            if input_components != 3 or target_components != 3:
                raise ValueError("Both input and target arrays must be 3-component for vector mode.")
        else:
            if target_components == 1 and target_component != 0:
                raise ValueError("target_component must be 0 for scalar target arrays.")
            if target_component < 0 or target_component >= target_components:
                raise ValueError("target_component out of range for target array.")
            if input_components == 1 and input_component != 0:
                raise ValueError("input_component must be 0 for scalar input arrays.")
            if input_components > 1 and input_component >= input_components:
                raise ValueError("input_component out of range for input array.")

        target_np = vtknp.vtk_to_numpy(target_array)

        if input_component is None:
            axis_vec = np.array(
                axis if axis is not None else [0.0, 1.0, 0.0], dtype=float
            ).reshape(-1)
            if axis_vec.size != 3:
                raise ValueError("axis must be a 3-element vector.")
            axis_norm = np.linalg.norm(axis_vec)
            if axis_norm <= 0.0:
                raise ValueError("axis must be a non-zero vector.")
            axis_unit = axis_vec / axis_norm
            origin_vec = np.array(
                origin if origin is not None else [0.0, 0.0, 0.0], dtype=float
            ).reshape(-1)
            if origin_vec.size != 3:
                raise ValueError("origin must be a 3-element vector.")

            if location == "point":
                points = input_mesh.GetPoints()
                input_pos = np.empty((0, 3)) if points is None else vtknp.vtk_to_numpy(points.GetData())
                target_points = target_mesh.GetPoints()
                target_pos = np.empty((0, 3)) if target_points is None else vtknp.vtk_to_numpy(target_points.GetData())
            else:
                input_pos = VTKGridHelper.compute_cell_centroids(input_mesh)
                target_pos = VTKGridHelper.compute_cell_centroids(target_mesh)

            if input_pos.size == 0:
                return

            input_vec = input_pos - origin_vec
            input_axial = np.dot(input_vec, axis_unit)
            input_radial_vec = input_vec - np.outer(input_axial, axis_unit)
            radial_norm = np.linalg.norm(input_radial_vec, axis=1)
            input_radial_unit = np.zeros_like(input_radial_vec)
            nonzero = radial_norm > 0.0
            input_radial_unit[nonzero] = input_radial_vec[nonzero] / radial_norm[nonzero, None]
            input_toroidal_unit = np.cross(axis_unit, input_radial_unit)

            target_vec = target_pos - origin_vec
            target_axial = np.dot(target_vec, axis_unit)
            target_radial_vec = target_vec - np.outer(target_axial, axis_unit)
            target_radial_norm = np.linalg.norm(target_radial_vec, axis=1)
            target_radial_unit = np.zeros_like(target_radial_vec)
            target_nonzero = target_radial_norm > 0.0
            target_radial_unit[target_nonzero] = target_radial_vec[target_nonzero] / target_radial_norm[target_nonzero, None]
            target_toroidal_unit = np.cross(axis_unit, target_radial_unit)

        for row in projection_map:
            target_id = row[0]
            if np.isnan(target_id):
                continue
            target_id_int = int(target_id)
            target_ids = row[1:]
            target_ids = target_ids[~np.isnan(target_ids)]
            if target_ids.size == 0:
                if target_np.ndim == 1:
                    target_np[target_id_int] = np.nan
                else:
                    target_np[target_id_int, target_component] = np.nan
                continue
            input_ids = target_ids.astype(int)
            if input_component is None:
                input_vecs = input_np[input_ids]
                vr = np.mean(np.einsum("ij,ij->i", input_vecs, input_radial_unit[input_ids]))
                vt = np.mean(np.einsum("ij,ij->i", input_vecs, input_toroidal_unit[input_ids]))
                va = np.mean(np.dot(input_vecs, axis_unit))
                target_np[target_id_int, :] = (
                    vr * target_radial_unit[target_id_int]
                    + vt * target_toroidal_unit[target_id_int]
                    + va * axis_unit
                )
            else:
                if input_components == 1:
                    value = float(np.mean(input_np[input_ids]))
                else:
                    if input_component < 0 or input_component >= input_components:
                        raise ValueError("input_component out of range for input array.")
                    value = float(np.mean(input_np[input_ids, input_component]))
                if target_np.ndim == 1:
                    target_np[target_id_int] = value
                else:
                    target_np[target_id_int, target_component] = value

        target_array.Modified()

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

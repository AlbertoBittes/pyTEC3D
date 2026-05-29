import os
from os import PathLike
from typing import Any, Dict, List, Optional, TYPE_CHECKING, Tuple, Union

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

from .constants import VTK_ARRAYS_TO_INIT
from .logging_utils import get_logger
from .validation_utils import ensure_pathlike, ensure_vtk_unstructured_grid


LOGGER = get_logger("apec_m.model_edit")

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import vtkmodules.vtkCommonDataModel as vtk_dm  # type: ignore[import-not-found]
    from .coordinate_system import CoordinateSystem


class ModelEdit:
    @staticmethod
    def _get_positions(
        mesh: "vtk.vtkUnstructuredGrid",
        location: str,
    ) -> np.ndarray:
        """Return node coordinates or cached cell centroids for a mesh location.

        Args:
            mesh: Input unstructured grid.
            location: "node" for points or "elem" for cell centroids.

        Returns:
            np.ndarray: Array of positions with shape (N, 3).
        """
        ensure_vtk_unstructured_grid(mesh)
        if location == "node":
            points = mesh.GetPoints()
            if points is None:
                return np.empty((0, 3))
            return vtknp.vtk_to_numpy(points.GetData())
        cell_data = mesh.GetCellData()
        centroids_array = cell_data.GetArray("cell_centroids") if cell_data is not None else None
        if centroids_array is None:
            raise ValueError("CellData array 'cell_centroids' not found.")
        return vtknp.vtk_to_numpy(centroids_array)

    @staticmethod
    def to_local_coordinates(
        mesh: "vtk.vtkUnstructuredGrid",
        coordinate_system: "CoordinateSystem",
        location: str = "node",
        selected_only: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Map mesh points or cell centroids into the provided coordinate system.

        Args:
            mesh: The VTK unstructured grid to sample.
            coordinate_system: Coordinate system used for local coordinates.
            location: "node" or "elem".
            selected_only: When True, use selection arrays to filter IDs.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (local_coordinates, ids)
        """
        from .selector import Selector

        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")

        positions = ModelEdit._get_positions(mesh, location)
        if positions.size == 0:
            return np.empty((0, 3)), np.empty((0,), dtype=int)

        if selected_only:
            Selector.ensure_selection_arrays(mesh)
            if location == "node":
                count = mesh.GetNumberOfPoints()
                data = mesh.GetPointData()
                selection_name = "selection_point"
            else:
                count = mesh.GetNumberOfCells()
                data = mesh.GetCellData()
                selection_name = "selection_cell"
            ids = np.arange(count, dtype=int)
            selection_array = (
                data.GetArray(selection_name) if data is not None else None
            )
            if selection_array is None:
                return positions, ids
            selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
            if selection_mask.size != ids.size:
                raise ValueError(
                    f"Selection array '{selection_name}' has unexpected length."
                )
            if np.any(selection_mask):
                positions = positions[selection_mask]
                ids = ids[selection_mask]
            else:
                return positions[:0], ids[:0]
        else:
            count = (
                mesh.GetNumberOfPoints()
                if location == "node"
                else mesh.GetNumberOfCells()
            )
            ids = np.arange(count, dtype=int)

        local_positions = coordinate_system.to_local(positions)
        return local_positions, ids

    @staticmethod
    def to_global_coordinates(
        mesh: "vtk.vtkUnstructuredGrid",
        local_positions: np.ndarray,
        coordinate_system: "CoordinateSystem",
        location: str = "node",
        selected_only: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Map local coordinates back into the global coordinate system.

        Args:
            mesh: The VTK unstructured grid used for IDs/selection.
            local_positions: Local coordinates to map into global coordinates.
            coordinate_system: Coordinate system used for the mapping.
            location: "node" or "elem".
            selected_only: When True, use selection arrays to filter IDs.

        Returns:
            Tuple[np.ndarray, np.ndarray]: (global_coordinates, ids)
        """
        from .selector import Selector

        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")

        count = (
            mesh.GetNumberOfPoints()
            if location == "node"
            else mesh.GetNumberOfCells()
        )
        ids = np.arange(count, dtype=int)

        if selected_only:
            Selector.ensure_selection_arrays(mesh)
            if location == "node":
                data = mesh.GetPointData()
                selection_name = "selection_point"
            else:
                data = mesh.GetCellData()
                selection_name = "selection_cell"
            selection_array = data.GetArray(selection_name) if data is not None else None
            if selection_array is None:
                selection_mask = np.ones((ids.size,), dtype=bool)
            else:
                selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
                if selection_mask.size != ids.size:
                    raise ValueError(
                        f"Selection array '{selection_name}' has unexpected length."
                    )
            if selection_array is not None and np.any(selection_mask):
                ids = ids[selection_mask]
            elif selection_array is not None:
                return np.empty((0, 3), dtype=float), ids[:0]

        local_positions = np.asarray(local_positions, dtype=float)
        if local_positions.shape[0] != ids.size:
            raise ValueError(
                "local_positions must match the number of selected ids."
            )

        global_positions = coordinate_system.to_global(local_positions)
        return global_positions, ids

    @staticmethod
    def map_vector_field_to_local_components(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        coordinate_system: "CoordinateSystem",
        location: str = "node",
        selected_only: bool = False,
    ) -> np.ndarray:
        """
        Map a 3-component vector field into local coordinate components.

        For cartesian, this returns local x/y/z. For cylindrical, r/theta/z.
        For spherical, r/theta/phi.

        Args:
            mesh: Mesh containing the vector array.
            array_name: Name of the vector array.
            coordinate_system: Coordinate system used for the mapping.
            location: "node" or "elem".
            selected_only: When True, map only selected ids.

        Returns:
            np.ndarray: Updated NumPy view of the target array.
        """
        data = mesh.GetPointData() if location == "node" else mesh.GetCellData()
        vtk_array = data.GetArray(array_name) if data is not None else None
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {location} data.")

        mapped = ModelEdit._map_vector_field_to_local_components(
            mesh,
            array_name,
            coordinate_system,
            location=location,
            selected_only=selected_only,
        )
        vecs = vtknp.vtk_to_numpy(vtk_array)
        if selected_only:
            from .selector import Selector

            Selector.ensure_selection_arrays(mesh)
            selection_name = "selection_point" if location == "node" else "selection_cell"
            selection_array = data.GetArray(selection_name) if data is not None else None
            if selection_array is None:
                return vecs
            selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
            if np.any(selection_mask):
                vecs[selection_mask] = mapped
        else:
            vecs[:] = mapped
        vtk_array.Modified()
        return vecs

    @staticmethod
    def _map_vector_field_to_local_components(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        coordinate_system: "CoordinateSystem",
        location: str = "node",
        selected_only: bool = False,
    ) -> np.ndarray:
        """Compute local vector components without writing them back to VTK.

        Args:
            mesh: Mesh containing the vector array.
            array_name: Name of the vector array.
            coordinate_system: Coordinate system used for the mapping.
            location: "node" or "elem".
            selected_only: When True, return values only for selected ids.

        Returns:
            np.ndarray: Mapped vector components.
        """
        from .selector import Selector

        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")

        data = mesh.GetPointData() if location == "node" else mesh.GetCellData()
        vtk_array = data.GetArray(array_name) if data is not None else None
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {location} data.")
        if vtk_array.GetNumberOfComponents() != 3:
            raise ValueError("Array must have exactly 3 components.")

        vecs = np.asarray(vtknp.vtk_to_numpy(vtk_array), dtype=float)
        positions = ModelEdit._get_positions(mesh, location)
        if positions.size == 0:
            return np.empty((0, 3), dtype=float) if selected_only else vecs.copy()

        if selected_only:
            Selector.ensure_selection_arrays(mesh)
            if location == "node":
                data = mesh.GetPointData()
                selection_name = "selection_point"
            else:
                data = mesh.GetCellData()
                selection_name = "selection_cell"
            selection_array = data.GetArray(selection_name) if data is not None else None
            if selection_array is None:
                return np.empty((0, 3), dtype=float)
            selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
            if selection_mask.size != vecs.shape[0]:
                raise ValueError(
                    f"Selection array '{selection_name}' has unexpected length."
                )
            selected_indices = np.flatnonzero(selection_mask)
            if selected_indices.size == 0:
                return np.empty((0, 3), dtype=float)
            selected_vecs = vecs[selected_indices].copy()
            selected_positions = positions[selected_indices]
            valid_mask = np.any(~np.isnan(selected_vecs), axis=1)
            if not np.any(valid_mask):
                return selected_vecs
            selected_vecs[valid_mask] = np.nan_to_num(
                selected_vecs[valid_mask],
                nan=0.0,
            )
            local_components = coordinate_system.vector_components_to_local(
                selected_vecs[valid_mask],
                selected_positions[valid_mask],
            )
            selected_vecs[valid_mask] = local_components
            return selected_vecs

        valid_mask = np.any(~np.isnan(vecs), axis=1)
        mapped = vecs.copy()
        if not np.any(valid_mask):
            return mapped

        mapped[valid_mask] = np.nan_to_num(mapped[valid_mask], nan=0.0)
        local_components = coordinate_system.vector_components_to_local(
            mapped[valid_mask], positions[valid_mask]
        )
        mapped[valid_mask] = local_components
        return mapped

    @staticmethod
    def map_vector_field_to_global_components(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        coordinate_system: "CoordinateSystem",
        location: str = "node",
        selected_only: bool = False,
    ) -> np.ndarray:
        """
        Map a 3-component vector field from local components to global components.

        Args:
            mesh: Mesh containing the array.
            array_name: Name of the array to rotate.
            coordinate_system: Coordinate system used for rotation.
            location: "node" or "elem".
            selected_only: When True, apply to selected IDs only.

        Returns:
            np.ndarray: The rotated array (NumPy view).
        """
        data = mesh.GetPointData() if location == "node" else mesh.GetCellData()
        vtk_array = data.GetArray(array_name) if data is not None else None
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {location} data.")

        mapped = ModelEdit._map_vector_field_to_global_components(
            mesh,
            array_name,
            coordinate_system,
            location=location,
            selected_only=selected_only,
        )
        vecs = vtknp.vtk_to_numpy(vtk_array)
        if selected_only:
            from .selector import Selector

            Selector.ensure_selection_arrays(mesh)
            selection_name = "selection_point" if location == "node" else "selection_cell"
            selection_array = data.GetArray(selection_name) if data is not None else None
            if selection_array is None:
                return vecs
            selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
            if np.any(selection_mask):
                vecs[selection_mask] = mapped
        else:
            vecs[:] = mapped
        vtk_array.Modified()
        return vecs

    @staticmethod
    def _map_vector_field_to_global_components(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        coordinate_system: "CoordinateSystem",
        location: str = "node",
        selected_only: bool = False,
    ) -> np.ndarray:
        """Compute global vector components without writing them back to VTK.

        Args:
            mesh: Mesh containing the vector array.
            array_name: Name of the vector array.
            coordinate_system: Coordinate system used for the mapping.
            location: "node" or "elem".
            selected_only: When True, return values only for selected ids.

        Returns:
            np.ndarray: Mapped global vector components.
        """
        from .selector import Selector

        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")

        data = mesh.GetPointData() if location == "node" else mesh.GetCellData()
        vtk_array = data.GetArray(array_name) if data is not None else None
        if vtk_array is None:
            raise ValueError(f"Array '{array_name}' not found in {location} data.")
        if vtk_array.GetNumberOfComponents() != 3:
            raise ValueError("Array must have exactly 3 components.")

        vecs = np.asarray(vtknp.vtk_to_numpy(vtk_array), dtype=float)
        positions = ModelEdit._get_positions(mesh, location)
        if positions.size == 0:
            return np.empty((0, 3), dtype=float) if selected_only else vecs.copy()

        if selected_only:
            Selector.ensure_selection_arrays(mesh)
            if location == "node":
                data = mesh.GetPointData()
                selection_name = "selection_point"
            else:
                data = mesh.GetCellData()
                selection_name = "selection_cell"
            selection_array = data.GetArray(selection_name) if data is not None else None
            if selection_array is None:
                return np.empty((0, 3), dtype=float)
            selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
            if selection_mask.size != vecs.shape[0]:
                raise ValueError(
                    f"Selection array '{selection_name}' has unexpected length."
                )
            selected_indices = np.flatnonzero(selection_mask)
            if selected_indices.size == 0:
                return np.empty((0, 3), dtype=float)
            selected_vecs = vecs[selected_indices].copy()
            selected_positions = positions[selected_indices]
            valid_mask = np.any(~np.isnan(selected_vecs), axis=1)
            if not np.any(valid_mask):
                return selected_vecs
            selected_vecs[valid_mask] = np.nan_to_num(
                selected_vecs[valid_mask],
                nan=0.0,
            )
            global_components = coordinate_system.vector_components_to_global(
                selected_vecs[valid_mask],
                selected_positions[valid_mask],
            )
            selected_vecs[valid_mask] = global_components
            return selected_vecs

        valid_mask = np.any(~np.isnan(vecs), axis=1)
        mapped = vecs.copy()
        if not np.any(valid_mask):
            return mapped

        mapped[valid_mask] = np.nan_to_num(mapped[valid_mask], nan=0.0)
        global_components = coordinate_system.vector_components_to_global(
            mapped[valid_mask], positions[valid_mask]
        )
        mapped[valid_mask] = global_components
        return mapped

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
        path = ensure_pathlike(vtk_file_path, "vtk_file_path")
        if not os.path.exists(path):
            raise FileNotFoundError(f"Mesh file '{path}' not found.")

        _, ext = os.path.splitext(path)
        ext = ext.lower()
        if ext == ".pvtu":
            reader = vtk.vtkXMLPUnstructuredGridReader()
        elif ext == ".vtu":
            reader = vtk.vtkXMLUnstructuredGridReader()
        else:
            raise ValueError(
                f"Unsupported mesh extension '{ext}'. Expected '.vtu' or '.pvtu'."
            )
        reader.SetFileName(path)
        reader.Update()
        mesh = ensure_vtk_unstructured_grid(reader.GetOutput(), arg_name="mesh")

        if mesh.GetNumberOfPoints() == 0:
            raise ValueError(f"Loaded mesh '{path}' has zero points.")

        if mesh.GetNumberOfCells() == 0:
            raise ValueError(f"Loaded mesh '{path}' has zero cells.")

        ModelEdit.store_cell_centroids(mesh)
        ModelEdit.store_cell_volumes(mesh)

        if reset_arrays:
            # Avoid circular import at module load time
            from .selector import Selector
            from .tree_gauge import TreeGauge

            ModelEdit.initialize_vtk_arrays(mesh)
            Selector.ensure_selection_arrays(mesh)
            TreeGauge.add_midside_nodes_mask_array(mesh)
            ModelEdit.reset_vtk_arrays_to_default(mesh)

        LOGGER.info("VTK unstructured grid loaded from '%s'.", path)

        return mesh

    @staticmethod
    def add_array(
        mesh: "vtk.vtkUnstructuredGrid",
        name: str,
        location: str = "node",
        components: int = 1,
        fill: Any = 0,
        dtype: Optional[int] = vtk.VTK_DOUBLE,
        values: Optional[np.ndarray] = None,
        replace: bool = False,
    ) -> "vtk.vtkDataArray":
        """
        Adds a VTK array to the mesh at the specified location.

        Args:
            mesh: The VTK unstructured grid to add the array to.
            name: Name of the array.
            location: "elem" or "node" (default "node").
            components: Number of components (default 1).
            fill: Default value to fill the array if values are not provided (default 0).
            dtype: VTK data type (default vtk.VTK_DOUBLE). If None, dtype validation is skipped.
            values: Explicit values to store in the array. Must match expected shape.
            replace: When True, remove and recreate the array if it already exists.

        Returns:
            vtk.vtkDataArray: The created or pre-existing VTK array.

        Raises:
            ValueError: If location is invalid, shapes/components mismatch, or dtype differs.
        """
        ensure_vtk_unstructured_grid(mesh)
        if location not in ["elem", "node"]:
            raise ValueError(f"Unknown location '{location}'. Use 'elem' or 'node'.")

        data_obj = mesh.GetCellData() if location == "elem" else mesh.GetPointData()
        num_entries = (
            mesh.GetNumberOfCells() if location == "elem" else mesh.GetNumberOfPoints()
        )
        num_components = max(1, components)
        expected_shape = (
            (num_entries,) if num_components == 1 else (num_entries, num_components)
        )

        vtk_array = data_obj.GetArray(name)
        if vtk_array is not None and not replace:
            return vtk_array

        if replace and vtk_array is not None:
            data_obj.RemoveArray(name)
            vtk_array = None

        values_np = None
        if values is not None:
            values_np = np.asarray(values)
            if values_np.shape == (num_entries, 1) and num_components == 1:
                values_np = values_np.reshape(num_entries)
            if values_np.shape != expected_shape:
                raise ValueError(
                    f"Values for array '{name}' must have shape {expected_shape}, got {values_np.shape}."
                )

        if values_np is not None:
            arr_to_store = values_np
        else:
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
        ensure_vtk_unstructured_grid(mesh)

        output_path = ensure_pathlike(output_file, "output_file")

        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(output_path)
        writer.SetInputData(mesh)
        if binary:
            writer.SetDataModeToBinary()
        else:
            writer.SetDataModeToAscii()

        write_status = writer.Write()
        if write_status != 1:
            raise IOError(f"Failed to write VTK mesh to '{output_path}'.")

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
        ensure_vtk_unstructured_grid(mesh)
        if arrays_to_init is None:
            arrays_to_init = VTK_ARRAYS_TO_INIT

        for name, info in arrays_to_init.items():
            ModelEdit.add_array(
                mesh,
                name=name,
                location=info["location"],
                components=info["components"],
                fill=info.get("fill", 0),
                dtype=info.get("dtype", vtk.VTK_DOUBLE),
            )

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

        Returns:
            None
        """
        ensure_vtk_unstructured_grid(mesh)
        if arrays_to_init is None:
            arrays_to_init = VTK_ARRAYS_TO_INIT

        cell_data = mesh.GetCellData()
        point_data = mesh.GetPointData()

        for name, info in arrays_to_init.items():
            data_obj = cell_data if info["location"] == "elem" else point_data
            vtk_array = data_obj.GetArray(name)
            if vtk_array is not None:
                expected_components = max(1, info["components"])
                expected_dtype = info.get("dtype", vtk.VTK_DOUBLE)
                if (
                    vtk_array.GetNumberOfComponents() != expected_components
                    or (
                        expected_dtype is not None
                        and vtk_array.GetDataType() != expected_dtype
                    )
                ):
                    ModelEdit.remove_array(
                        mesh, name, location=info["location"]
                    )
                    ModelEdit.add_array(
                        mesh,
                        name=name,
                        location=info["location"],
                        components=expected_components,
                        fill=info.get("fill", 0),
                        dtype=expected_dtype,
                        replace=True,
                    )
                    continue
            if vtk_array is not None and "fill" in info:
                num = (
                    mesh.GetNumberOfCells()
                    if info["location"] == "elem"
                    else mesh.GetNumberOfPoints()
                )
                shape = (num, info["components"]) if info["components"] > 1 else num
                arr = np.full(shape, info["fill"])
                # Overwrite the VTK array with the default values
                vtk_np = vtknp.vtk_to_numpy(vtk_array)
                vtk_np[:] = arr
                vtk_array.Modified()

    @staticmethod
    def remove_array(
        mesh: "vtk.vtkUnstructuredGrid",
        array_name: str,
        location: Optional[str] = None,
    ) -> None:
        """
        Safely removes an array from CellData or PointData of the mesh if it exists.

        Args:
            mesh: The mesh from which to remove the array.
            array_name: The name of the array to remove.
            location: "elem", "node", or None. If None, guesses from the data that contains it.

        Side Effects:
            Removes the array from the detected or specified data location if present.

        Returns:
            None
        """
        ensure_vtk_unstructured_grid(mesh)
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
                LOGGER.warning("Array '%s' not found in CellData or PointData.", array_name)
                return
        elif location == "elem":
            data = mesh.GetCellData()
        elif location == "node":
            data = mesh.GetPointData()
        else:
            raise ValueError(
                f"Unknown location '{location}'. Use 'elem', 'node', or None."
            )
        if data.HasArray(array_name):
            data.RemoveArray(array_name)

    @staticmethod
    def element_solution(
        mesh: "vtk.vtkUnstructuredGrid",
        point_array_name: str,
    ) -> "vtk.vtkDataArray":
        """
        Project a point-data array to cells and store it as `<name>e` in CellData.

        Args:
            mesh: Input mesh.
            point_array_name: Name of the point-data array to transfer.

        Returns:
            vtkDataArray: The created/updated cell-data array.

        Raises:
            ValueError: If the source point-data array does not exist.
        """
        ensure_vtk_unstructured_grid(mesh)
        point_data = mesh.GetPointData()
        source_array = point_data.GetArray(point_array_name)
        if source_array is None:
            raise ValueError(
                f"PointData array '{point_array_name}' not found on the mesh."
            )

        p2c = vtk.vtkPointDataToCellData()
        p2c.SetInputData(mesh)
        p2c.SetProcessAllArrays(False)
        p2c.AddPointDataArray(point_array_name)
        p2c.PassPointDataOff()
        p2c.Update()

        output = p2c.GetOutput()
        projected_array = output.GetCellData().GetArray(point_array_name)
        if projected_array is None:
            raise ValueError(
                f"Could not project PointData array '{point_array_name}' to CellData."
            )

        target_name = f"{point_array_name}e"
        projected_np = vtknp.vtk_to_numpy(projected_array)
        n_components = projected_array.GetNumberOfComponents()
        dtype = projected_array.GetDataType()
        cell_data = mesh.GetCellData()
        target_array = cell_data.GetArray(target_name)

        if target_array is None:
            ModelEdit.add_array(
                mesh,
                name=target_name,
                location="elem",
                components=n_components,
                dtype=dtype,
                values=projected_np,
                replace=False,
            )
            target_array = cell_data.GetArray(target_name)
            target_array.Modified()
            return target_array

        # If shape/components match, update in-place.
        if target_array.GetNumberOfComponents() == n_components:
            target_np = vtknp.vtk_to_numpy(target_array)
            if target_np.shape == projected_np.shape:
                target_np[:] = projected_np
                target_array.Modified()
                return target_array

        # Recreate target array with projected values when layout changed.
        cell_data.RemoveArray(target_name)
        ModelEdit.add_array(
            mesh,
            name=target_name,
            location="elem",
            components=n_components,
            dtype=dtype,
            values=projected_np,
            replace=False,
        )

        target_array = cell_data.GetArray(target_name)
        target_array.Modified()
        return target_array

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
        ensure_vtk_unstructured_grid(mesh)
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

        Returns:
            None
        """
        np_centers = ModelEdit.compute_cell_centroids(mesh)
        ModelEdit.add_array(
            mesh,
            name=array_name,
            location="elem",
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
            Uses vtkCellSizeFilter so all supported cell types are handled
            consistently (including non-triangular 2D cells).
        """
        ensure_vtk_unstructured_grid(mesh)
        n_cells = mesh.GetNumberOfCells()
        if n_cells == 0:
            return np.zeros((0,), dtype=float)

        size_filter = vtk.vtkCellSizeFilter()
        size_filter.SetInputData(mesh)
        size_filter.SetComputeVertexCount(False)
        size_filter.SetComputeLength(True)
        size_filter.SetComputeArea(True)
        size_filter.SetComputeVolume(True)
        size_filter.SetComputeSum(False)

        length_name = "__cell_length_tmp__"
        area_name = "__cell_area_tmp__"
        volume_name = "__cell_volume_tmp__"
        size_filter.SetLengthArrayName(length_name)
        size_filter.SetAreaArrayName(area_name)
        size_filter.SetVolumeArrayName(volume_name)
        size_filter.Update()

        out_mesh = size_filter.GetOutput()
        out_cell_data = out_mesh.GetCellData()

        length_array = out_cell_data.GetArray(length_name)
        area_array = out_cell_data.GetArray(area_name)
        volume_array = out_cell_data.GetArray(volume_name)

        lengths = (
            np.asarray(vtknp.vtk_to_numpy(length_array), dtype=float).reshape(-1)
            if length_array is not None
            else np.zeros((n_cells,), dtype=float)
        )
        areas = (
            np.asarray(vtknp.vtk_to_numpy(area_array), dtype=float).reshape(-1)
            if area_array is not None
            else np.zeros((n_cells,), dtype=float)
        )
        volumes = (
            np.asarray(vtknp.vtk_to_numpy(volume_array), dtype=float).reshape(-1)
            if volume_array is not None
            else np.zeros((n_cells,), dtype=float)
        )

        if (
            lengths.shape[0] != n_cells
            or areas.shape[0] != n_cells
            or volumes.shape[0] != n_cells
        ):
            raise ValueError("vtkCellSizeFilter produced arrays with unexpected length.")

        # Select area for 2D cells and volume for 3D cells.
        cell_dims = np.fromiter(
            (mesh.GetCell(i).GetCellDimension() for i in range(n_cells)),
            dtype=np.int8,
            count=n_cells,
        )

        values = np.zeros((n_cells,), dtype=float)
        values[cell_dims >= 3] = volumes[cell_dims >= 3]
        values[cell_dims == 2] = areas[cell_dims == 2]
        return values

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

        Returns:
            None
        """
        volumes = ModelEdit.compute_cell_volumes(mesh)
        ModelEdit.add_array(
            mesh,
            name=array_name,
            location="elem",
            components=1,
            dtype=vtk.VTK_DOUBLE,
            values=volumes,
            replace=True,
        )

    @staticmethod
    def add_component(
        mesh: "vtk.vtkUnstructuredGrid", component_name: str, type: str = "elem"
    ) -> "vtk.vtkDataArray":
        """
        Adds a component array marking the current selection.

        Args:
            mesh: The mesh to operate on.
            component_name: Name of the component array to create.
            type: "elem" to use CellData selection, "node" for PointData selection.

        Returns:
            vtkDataArray: The created component array.
        """
        ensure_vtk_unstructured_grid(mesh)
        if type not in ["elem", "node"]:
            raise ValueError("type must be 'elem' or 'node'.")

        # Avoid circular import at module load time
        from .selector import Selector

        Selector.ensure_selection_arrays(mesh)
        data_obj = mesh.GetCellData() if type == "elem" else mesh.GetPointData()
        selection_array = data_obj.GetArray(
            "selection_cell" if type == "elem" else "selection_point"
        )
        if selection_array is None:
            raise ValueError(
                "Selection array not found; ensure selection arrays are initialized."
            )

        selection_mask = vtknp.vtk_to_numpy(selection_array)
        new_array = ModelEdit.add_array(
            mesh,
            name=component_name,
            location=type,
            components=1,
            fill=0,
            dtype=vtk.VTK_UNSIGNED_CHAR,
        )
        new_array_vtk = data_obj.GetArray(component_name)
        if new_array_vtk is None:
            raise RuntimeError(
                f"add_array did not create array '{component_name}' as expected."
            )
        target_np = vtknp.vtk_to_numpy(new_array_vtk)
        target_np[:] = selection_mask
        new_array_vtk.Modified()
        return new_array_vtk

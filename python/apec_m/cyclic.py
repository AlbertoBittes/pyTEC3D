from typing import Optional, Sequence, TYPE_CHECKING

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp


from .selector import Selector
from .model_edit import ModelEdit
from .constants import CONSTANTS_COLUMNS_TYPES

from .coordinate_system import CoordinateSystem, CylindricalCoordinateSystem


class Cyclic:
    @staticmethod
    def build_toroidal_slices(
        input_mesh: "vtk.vtkUnstructuredGrid",
        target_mesh: "vtk.vtkUnstructuredGrid",
        location: str = "node",
        tolerance: float = 1e-6,
        n_slices: int = 1,
        coordinate_system: Optional["CoordinateSystem"] = None,
    ) -> np.ndarray:
        """
        Build a toroidal projection map from a 2D input mesh to a 3D target mesh.

        Each input point/cell centroid is associated with target points/cell centroids that share
        the same radial distance from the rotation axis and the same axial coordinate within
        the given tolerance. For points on the rotation axis (r=0), the map assumes a single
        match. The tolerance is increased until each input ID has at least n_slices corresponding
        target IDs, then matches are sorted by poloidal distance and trimmed to n_slices.
        The resulting map is returned as a 2D NumPy array where the first column holds input IDs and
        the remaining columns hold the corresponding target IDs.
        If selection arrays exist on the input/target mesh for the chosen location, only selected
        IDs are used.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            target_mesh: Target unstructured grid (typically 3D).
            location: "node" to map point IDs or "elem" to map cell IDs.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            n_slices: Minimum number of corresponding target IDs per input ID.
            coordinate_system: Cylindrical coordinate system used for local coordinates.

        Returns:
            np.ndarray: Integer array with shape (num_input, 1 + n_slices), where
                column 0 stores 2D entity IDs and columns 1..n_slices store the
                corresponding toroidal slice IDs in the 3D mesh.

        Raises:
            ValueError: If inputs are invalid or the input mesh is not 2D.
        """
        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")
        if tolerance <= 0:
            raise ValueError("tolerance must be positive.")
        if not isinstance(n_slices, int) or n_slices <= 0:
            raise ValueError("n_slices must be a positive integer.")

        if (
            input_mesh.GetNumberOfCells() > 0
            and input_mesh.GetMaxSpatialDimension() > 2
        ):
            raise ValueError("Input mesh must be 2D.")

        if coordinate_system is None:
            from .coordinate_system import CylindricalCoordinateSystem

            coordinate_system = CylindricalCoordinateSystem()
        if coordinate_system.type != "cylindrical":
            raise ValueError(
                "coordinate_system must be cylindrical for toroidal projection."
            )

        from .selector import Selector

        Selector.ensure_selection_arrays(input_mesh)
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
            return np.empty((0, 1 + n_slices), dtype=np.int64)

        input_radius = input_local[:, 0]
        input_axial = input_local[:, 2]
        target_radius = target_local[:, 0]
        target_axial = target_local[:, 2]

        target_count = target_radius.shape[0]
        if target_count == 0:
            raise ValueError("Target mesh has no points/cells to map.")
        if n_slices > target_count:
            raise ValueError(
                "n_slices cannot exceed the number of available target points/cells."
            )

        num_input = input_radius.shape[0]
        tolerance_current = tolerance
        counts = np.empty(num_input, dtype=int)

        axis_mask = input_radius <= tolerance
        required_counts = np.where(axis_mask, 1, n_slices)

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

        mapping = np.full((num_input, 1 + n_slices), -1, dtype=np.int64)
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
                    (n_slices,), candidate_ids[0], dtype=candidate_ids.dtype
                )
            else:
                poloidal_dist = np.hypot(dr[candidate_idx], dz[candidate_idx])
                if candidate_idx.size > n_slices:
                    pick = np.argpartition(poloidal_dist, n_slices - 1)[:n_slices]
                    pick = pick[np.argsort(poloidal_dist[pick])]
                    trimmed = candidate_ids[pick]
                else:
                    order = np.argsort(poloidal_dist)
                    trimmed = candidate_ids[order]
            mapping[row_idx, 1 : 1 + trimmed.size] = trimmed[:n_slices].astype(
                np.int64, copy=False
            )

        return mapping

    @staticmethod
    def cpcyc(
        mesh: "vtk.vtkUnstructuredGrid",
        coordinate_system: "CoordinateSystem",
        theta_rotation: float,
        array_name: str = "CP_0",
        tolerance: float = 1e-6,
        variables: Optional[Sequence[str]] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Assign incremental pair IDs for selected points based on rotated positions.

        Selected points are rotated by theta_rotation in the provided cylindrical coordinate
        system, mapped to their nearest mesh points,
        and both points in each pair receive the same incrementing value in the array. Existing
        non-zero values are preserved and skipped when assigning new IDs.

        Args:
            mesh: VTK unstructured grid to operate on.
            coordinate_system: Cylindrical coordinate system used for the rotation.
            theta_rotation: Rotation (degrees) added to the local theta coordinate.
            array_name: Point data array name to store pair IDs (default "CP_0").
            tolerance: Absolute tolerance for theta == 0 check (default 1e-6).
            variables: Optional sequence of DOF variable names to couple.

        Returns:
            Optional[pd.DataFrame]: DataFrame with columns ["CP_ID", "VARIABLE"] when
                variables are provided, otherwise None.
        """
        if coordinate_system.type != "cylindrical":
            raise ValueError("coordinate_system must be cylindrical.")
        if tolerance <= 0:
            raise ValueError("tolerance must be positive.")

        variables_list: Optional[list[str]]
        if variables is None:
            variables_list = None
        elif isinstance(variables, str):
            variables_list = [variables]
        else:
            variables_list = [str(var) for var in variables]

        Selector.ensure_selection_arrays(mesh)

        point_data = mesh.GetPointData()
        selection_array = point_data.GetArray("selection_point")
        selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
        if not np.any(selection_mask):
            if variables_list is None:
                return None
            return pd.DataFrame(columns=["CP_ID", "VARIABLE"]).astype(
                {"CP_ID": int, "VARIABLE": str}
            )

        points = mesh.GetPoints()
        coords = vtknp.vtk_to_numpy(points.GetData())
        selected_ids = np.flatnonzero(selection_mask)
        selected_coords = coords[selected_ids]
        local_coords = coordinate_system.to_local(selected_coords)
        if not np.all(np.abs(local_coords[:, 1]) <= tolerance):
            raise ValueError("Selected points must have theta == 0 within tolerance.")
        local_coords[:, 1] += theta_rotation
        rotated_coords = coordinate_system.to_global(local_coords)

        locator = vtk.vtkPointLocator()
        locator.SetDataSet(mesh)
        locator.BuildLocator()

        cp_array = point_data.GetArray(array_name)
        cp_np = vtknp.vtk_to_numpy(cp_array)

        existing_values = set(np.unique(cp_np[cp_np != 0]).astype(int).tolist())
        next_value = 1
        while next_value in existing_values:
            next_value += 1

        new_cp_ids: list[int] = []
        assigned = set()
        for src_id, target_coord in zip(selected_ids, rotated_coords):
            if src_id in assigned or cp_np[src_id] != 0:
                continue
            tgt_id = locator.FindClosestPoint(target_coord)
            if tgt_id < 0 or tgt_id == src_id:
                continue
            if tgt_id in assigned or cp_np[tgt_id] != 0:
                continue

            cp_id = next_value
            cp_np[src_id] = cp_id
            cp_np[tgt_id] = cp_id
            assigned.add(src_id)
            assigned.add(tgt_id)
            existing_values.add(cp_id)
            new_cp_ids.append(cp_id)
            next_value += 1
            while next_value in existing_values:
                next_value += 1

        cp_array.Modified()

        if variables_list is None:
            return None

        rows = [(cp_id, var) for cp_id in new_cp_ids for var in variables_list]
        return pd.DataFrame(rows, columns=["CP_ID", "VARIABLE"]).astype(
            {"CP_ID": int, "VARIABLE": str}
        )

    @staticmethod
    def coil_winding(
        mesh: "vtk.vtkUnstructuredGrid",
        coordinate_system: "CoordinateSystem",
        winding_parameters: Optional[Sequence[float]] = None,
    ) -> pd.DataFrame:
        """
        Assign per-cell real constants for selected coil elements in a cylindrical frame.

        For each selected cell, a new increasing real-constant ID is written to the mesh
        cell-data array, and a matching row is returned in a constants DataFrame with
        directional cosines TX/TY/TZ aligned with the local theta axis.

        Args:
            mesh: Mesh containing selected cells (selection_cell == 1).
            coordinate_system: Cylindrical coordinate system used for theta direction.
            winding_parameters: Optional tuple (SC, NC, VC, R, SYM).

        Returns:
            pd.DataFrame: Rows ready to assign into apec.constants.
        """
        if coordinate_system.type != "cylindrical":
            raise ValueError("coordinate_system must be cylindrical.")

        if winding_parameters is not None and len(winding_parameters) != 5:
            raise ValueError(
                "winding_values must be a tuple of 5 values: (SC, NC, VC, R, SYM)."
            )

        Selector.ensure_selection_arrays(mesh)

        cell_data = mesh.GetCellData()
        selection_array = cell_data.GetArray("selection_cell")
        selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
        selected_ids = np.flatnonzero(selection_mask)

        dtype_map = dict(CONSTANTS_COLUMNS_TYPES)
        columns = list(dtype_map.keys())
        if selected_ids.size == 0:
            return pd.DataFrame(columns=columns).astype(dtype_map)

        real_constant_np = vtknp.vtk_to_numpy(cell_data.GetArray("CONSTANT"))
        current_max = int(np.max(real_constant_np)) if real_constant_np.size else 0
        new_constants = np.arange(
            current_max + 1, current_max + 1 + selected_ids.size, dtype=int
        )
        real_constant_np[selected_ids] = new_constants
        cell_data.GetArray("CONSTANT").Modified()

        centroids = ModelEdit.compute_cell_centroids(mesh)
        selected_centroids = centroids[selected_ids]
        local_components = np.zeros((selected_ids.size, 3), dtype=float)
        local_components[:, 1] = 1.0  # theta axis
        coil_directions = coordinate_system.vector_components_to_global(
            local_components, selected_centroids
        )

        constants_df = pd.DataFrame(
            0.0, index=np.arange(selected_ids.size), columns=columns
        )
        constants_df["CONSTANT"] = new_constants
        constants_df["TX"] = coil_directions[:, 0]
        constants_df["TY"] = coil_directions[:, 1]
        constants_df["TZ"] = coil_directions[:, 2]

        if winding_parameters is not None:
            sc, nc, vc, r_value, sym = winding_parameters
            constants_df["SC"] = sc
            constants_df["NC"] = nc
            constants_df["VC"] = vc
            constants_df["R"] = r_value
            constants_df["SYM"] = sym

        return constants_df.astype(dtype_map)

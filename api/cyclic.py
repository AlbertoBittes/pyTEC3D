from typing import Optional, Sequence, TYPE_CHECKING

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp


from selector import Selector
from model_edit import ModelEdit
from constants import REAL_CONSTANTS_COLUMNS_TYPES

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    from coordinate_system import CoordinateSystem


class Cyclic:
    @staticmethod
    def _to_slice_layout(projection_map: np.ndarray) -> np.ndarray:
        if projection_map.ndim != 2:
            raise ValueError("projection_map must be a 2D array.")
        if projection_map.shape[0] == 0 or projection_map.shape[1] == 0:
            return projection_map
        if projection_map.shape[0] > projection_map.shape[1]:
            return projection_map.T
        return projection_map

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
        Build toroidal `slice_x` arrays from a 2D input mesh to a 3D target mesh.

        The output stores one input ID per column (full input size). Row `slice_0` (row 0)
        holds 2D IDs for selected inputs at theta=0, while non-selected inputs remain `-1`.
        Rows `slice_1..slice_{n-1}` hold matching IDs from the target mesh, ordered by
        increasing wrapped toroidal angle (0..360) when multiple candidates are available.
        The returned array is initialized with `-1`, and only selected input IDs are populated.
        Integer arrays `slice_0..slice_{n-1}` are also written on the target mesh at the
        requested location, with values equal to the reference 2D ID.
        Matching is performed on cylindrical radius and axial coordinate within tolerance.
        Tolerance is expanded until each non-axis 2D ID has at least `n_slices` matches.
        If selection arrays exist on either mesh for the chosen location, only selected IDs
        are used.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            target_mesh: Target unstructured grid (typically 3D).
            location: "node" to map point IDs or "elem" to map cell IDs.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            n_slices: Total number of toroidal slices to build, including `slice_0`.
            coordinate_system: Cylindrical coordinate system used for local coordinates.

        Returns:
            np.ndarray: Integer array with shape (n_slices, num_input_total_2d).

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
            from coordinate_system import CylindricalCoordinateSystem

            coordinate_system = CylindricalCoordinateSystem()
        if coordinate_system.type != "cylindrical":
            raise ValueError(
                "coordinate_system must be cylindrical for toroidal projection."
            )

        input_count = (
            input_mesh.GetNumberOfPoints()
            if location == "node"
            else input_mesh.GetNumberOfCells()
        )
        slices = np.full((n_slices, input_count), -1, dtype=np.int64)

        Selector.ensure_selection_arrays(input_mesh)
        Selector.ensure_selection_arrays(target_mesh)
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
            target_count_total = (
                target_mesh.GetNumberOfPoints()
                if location == "node"
                else target_mesh.GetNumberOfCells()
            )
            slice_arrays = np.full((n_slices, target_count_total), -1, dtype=np.int64)
            for slice_idx in range(n_slices):
                slice_name = f"slice_{slice_idx}"
                data = (
                    target_mesh.GetPointData()
                    if location == "node"
                    else target_mesh.GetCellData()
                )
                data.RemoveArray(slice_name)
                vtk_arr = vtknp.numpy_to_vtk(
                    slice_arrays[slice_idx], deep=True, array_type=vtk.VTK_LONG_LONG
                )
                vtk_arr.SetName(slice_name)
                data.AddArray(vtk_arr)
            return slices

        input_radius = input_local[:, 0]
        input_axial = input_local[:, 2]
        target_radius = target_local[:, 0]
        target_theta = np.mod(target_local[:, 1], 360.0)
        target_axial = target_local[:, 2]

        target_count = target_radius.shape[0]
        if target_count == 0:
            raise ValueError("Target mesh has no points/cells to map.")
        if n_slices > target_count and np.any(input_radius > tolerance):
            raise ValueError(
                "n_slices cannot exceed the number of available target points/cells."
            )

        target_count_total = (
            target_mesh.GetNumberOfPoints()
            if location == "node"
            else target_mesh.GetNumberOfCells()
        )
        slice_arrays = np.full((n_slices, target_count_total), -1, dtype=np.int64)

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

        slices[0, input_ids.astype(int)] = input_ids.astype(np.int64)

        for row_idx in range(num_input):
            dr = np.abs(target_radius - input_radius[row_idx])
            dz = np.abs(target_axial - input_axial[row_idx])
            mask = (dr <= tolerance_current) & (dz <= tolerance_current)
            candidate_idx = np.where(mask)[0]
            if candidate_idx.size == 0:
                continue
            candidate_ids = target_ids[candidate_idx]
            input_id = int(input_ids[row_idx])
            if axis_mask[row_idx]:
                axis_target_id = int(candidate_ids[0])
                if n_slices > 1:
                    slices[1:, input_id] = axis_target_id
                for slice_idx in range(n_slices):
                    if 0 <= axis_target_id < target_count_total:
                        slice_arrays[slice_idx, axis_target_id] = input_id
                continue

            candidate_theta = target_theta[candidate_idx]
            poloidal_dist = np.hypot(dr[candidate_idx], dz[candidate_idx])
            ordered_candidate_idx = np.lexsort((poloidal_dist, candidate_theta))
            ordered_candidates = candidate_ids[ordered_candidate_idx][:n_slices].astype(
                np.int64
            )
            slice0_target_id = int(ordered_candidates[0])
            if 0 <= slice0_target_id < target_count_total:
                slice_arrays[0, slice0_target_id] = input_id
            remaining = ordered_candidates[1:]
            if remaining.size > 0:
                slices[1 : 1 + remaining.size, input_id] = remaining
                for offset, target_id in enumerate(remaining, start=1):
                    target_idx = int(target_id)
                    if 0 <= target_idx < target_count_total:
                        slice_arrays[offset, target_idx] = input_id

        data = (
            target_mesh.GetPointData()
            if location == "node"
            else target_mesh.GetCellData()
        )
        for slice_idx in range(n_slices):
            slice_name = f"slice_{slice_idx}"
            data.RemoveArray(slice_name)
            vtk_arr = vtknp.numpy_to_vtk(
                slice_arrays[slice_idx], deep=True, array_type=vtk.VTK_LONG_LONG
            )
            vtk_arr.SetName(slice_name)
            data.AddArray(vtk_arr)

        return slices

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
            Optional[pd.DataFrame]: DataFrame with columns ["cp_id", "variable"] when
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
            return pd.DataFrame(columns=["cp_id", "variable"]).astype(
                {"cp_id": int, "variable": str}
            )

        points = mesh.GetPoints()
        coords = vtknp.vtk_to_numpy(points.GetData())
        selected_ids = np.flatnonzero(selection_mask)
        selected_coords = coords[selected_ids]
        local_coords = coordinate_system.to_local(selected_coords)
        if not np.all(np.abs(local_coords[:, 1]) <= tolerance):
            raise ValueError(
                "Selected points must have theta == 0 within tolerance."
            )
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
        return pd.DataFrame(rows, columns=["cp_id", "variable"]).astype(
            {"cp_id": int, "variable": str}
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
        cell-data array, and a matching row is returned in a real_constants DataFrame with
        directional cosines TX/TY/TZ aligned with the local theta axis.

        Args:
            mesh: Mesh containing selected cells (selection_cell == 1).
            coordinate_system: Cylindrical coordinate system used for theta direction.
            winding_parameters: Optional tuple (SC, NC, VC, R, SYM).

        Returns:
            pd.DataFrame: Rows ready to assign into apec.real_constants.
        """
        if coordinate_system.type != "cylindrical":
            raise ValueError("coordinate_system must be cylindrical.")

        if winding_parameters is not None and len(winding_parameters) != 5:
            raise ValueError("winding_values must be a tuple of 5 values: (SC, NC, VC, R, SYM).")

        Selector.ensure_selection_arrays(mesh)

        cell_data = mesh.GetCellData()
        selection_array = cell_data.GetArray("selection_cell")
        selection_mask = vtknp.vtk_to_numpy(selection_array).astype(bool)
        selected_ids = np.flatnonzero(selection_mask)

        dtype_map = dict(REAL_CONSTANTS_COLUMNS_TYPES)
        columns = list(dtype_map.keys())
        if selected_ids.size == 0:
            return pd.DataFrame(columns=columns).astype(dtype_map)

        real_constant_np = vtknp.vtk_to_numpy(cell_data.GetArray("ansys_real_constant"))
        current_max = int(np.max(real_constant_np)) if real_constant_np.size else 0
        new_real_constants = np.arange(
            current_max + 1, current_max + 1 + selected_ids.size, dtype=int
        )
        real_constant_np[selected_ids] = new_real_constants
        cell_data.GetArray("ansys_real_constant").Modified()

        centroids = ModelEdit.compute_cell_centroids(mesh)
        selected_centroids = centroids[selected_ids]
        local_components = np.zeros((selected_ids.size, 3), dtype=float)
        local_components[:, 1] = 1.0  # theta axis
        coil_directions = coordinate_system.vector_components_to_global(
            local_components, selected_centroids
        )

        real_constants_df = pd.DataFrame(0.0, index=np.arange(selected_ids.size), columns=columns)
        real_constants_df["ansys_real_constant"] = new_real_constants
        real_constants_df["TX"] = coil_directions[:, 0]
        real_constants_df["TY"] = coil_directions[:, 1]
        real_constants_df["TZ"] = coil_directions[:, 2]

        if winding_parameters is not None:
            sc, nc, vc, r_value, sym = winding_parameters
            real_constants_df["SC"] = sc
            real_constants_df["NC"] = nc
            real_constants_df["VC"] = vc
            real_constants_df["R"] = r_value
            real_constants_df["SYM"] = sym

        return real_constants_df.astype(dtype_map)

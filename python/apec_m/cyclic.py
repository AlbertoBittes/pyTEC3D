from typing import Optional, Sequence, TYPE_CHECKING

import numpy as np
import pandas as pd
import vtk
import vtk.util.numpy_support as vtknp


from .selector import Selector
from .model_edit import ModelEdit
from .constants import CONSTANTS_COLUMNS_TYPES
from .validation_utils import ensure_vtk_unstructured_grid

from .coordinate_system import CoordinateSystem, CylindricalCoordinateSystem


class Cyclic:
    @staticmethod
    def poloidal_flux(
        mesh: "vtk.vtkUnstructuredGrid",
        coordinate_system: "CylindricalCoordinateSystem",
    ):
        """
        Compute poloidal flux from available vector-potential arrays.

        - If PointData `A` exists, write `PSI`.
        - If CellData `Ae` exists, write `PSIe`.

        Args:
            mesh: Target mesh.
            coordinate_system: Cylindrical coordinate system.

        Returns:
            dict: Contains `"node"` and/or `"elem"` keys for computed flux arrays.
        """
        ensure_vtk_unstructured_grid(mesh)
        if coordinate_system.type != "cylindrical":
            raise ValueError("coordinate_system must be cylindrical.")
        spatial_dim = mesh.GetMaxSpatialDimension()

        def _compute_and_store(
            data: "vtk.vtkDataSetAttributes",
            positions: np.ndarray,
            source_name: str,
            target_name: str,
            target_location: str,
        ) -> np.ndarray:
            source_vtk = data.GetArray(source_name)
            if source_vtk.GetNumberOfComponents() != 3:
                raise ValueError(
                    f"Array '{source_name}' must have 3 components."
                )

            source_np = np.asarray(vtknp.vtk_to_numpy(source_vtk), dtype=float)
            if source_np.ndim != 2 or source_np.shape[1] != 3:
                raise ValueError(
                    f"Array '{source_name}' has invalid shape {source_np.shape}; expected (N, 3)."
                )

            if spatial_dim != 3:
                radial = np.asarray(positions[:, 0], dtype=float)
                flux_values = radial * source_np[:, 2]
            else:
                local_positions = coordinate_system.to_local(positions)
                local_vectors = coordinate_system.vector_components_to_local(
                    source_np, positions
                )
                radial = local_positions[:, 0]
                a_theta = local_vectors[:, 1]
                flux_values = radial * a_theta

            target_vtk = data.GetArray(target_name)
            if target_vtk is None or target_vtk.GetNumberOfComponents() != 1:
                ModelEdit.add_array(
                    mesh,
                    name=target_name,
                    location=target_location,
                    components=1,
                    dtype=vtk.VTK_DOUBLE,
                    replace=True,
                )
                target_vtk = data.GetArray(target_name)

            target_np = vtknp.vtk_to_numpy(target_vtk)
            target_np[:] = flux_values
            target_vtk.Modified()
            return target_np

        results = {}

        point_data = mesh.GetPointData()
        if point_data.GetArray("A") is not None:
            node_positions = vtknp.vtk_to_numpy(mesh.GetPoints().GetData())
            results["node"] = _compute_and_store(
                data=point_data,
                positions=node_positions,
                source_name="A",
                target_name="PSI",
                target_location="node",
            )

        cell_data = mesh.GetCellData()
        if cell_data.GetArray("Ae") is not None:
            centroids_array = cell_data.GetArray("cell_centroids")
            cell_positions = np.asarray(vtknp.vtk_to_numpy(centroids_array), dtype=float)
            results["elem"] = _compute_and_store(
                data=cell_data,
                positions=cell_positions,
                source_name="Ae",
                target_name="PSIe",
                target_location="elem",
            )

        if not results:
            raise ValueError(
                "No source arrays found for poloidal flux. "
                "Expected PointData 'A' and/or CellData 'Ae'."
            )

        return results

    @staticmethod
    def build_toroidal_slices(
        input_mesh: "vtk.vtkUnstructuredGrid",
        target_mesh: "vtk.vtkUnstructuredGrid",
        location: str = "node",
        tolerance: float = 1e-6,
        extrusion_angle: float = 360.0,
        n_divisions: Optional[int] = None,
        coordinate_system: Optional["CoordinateSystem"] = None,
    ) -> np.ndarray:
        """
        Build a toroidal projection map from a 2D input mesh to a target mesh.

        Each input point/cell centroid is associated with target points/cell centroids that
        share the same radial distance from the rotation axis and the same axial coordinate.
        For 3D target meshes, the target is assumed to be a toroidal extrusion of the input
        mesh over ``extrusion_angle`` degrees with ``n_divisions`` elements in the toroidal
        direction. Node maps therefore contain ``n_divisions + 1`` matches per source node,
        while element maps contain ``n_divisions`` matches per source element. For 2D target
        meshes there is exactly one match per source ID.

        The resulting map is returned as a 2D NumPy array where the first column holds input IDs
        and the remaining columns hold corresponding target IDs. For 3D targets, each target
        column corresponds to one toroidal coordinate station.
        If selection arrays exist on the input/target mesh for the chosen location, only selected
        IDs are used.

        Args:
            input_mesh: 2D unstructured grid used to generate a surface of revolution.
            target_mesh: Target unstructured grid (typically 3D).
            location: "node" to map point IDs or "elem" to map cell IDs.
            tolerance: Absolute tolerance for matching radii and axial coordinates.
            extrusion_angle: Total toroidal extrusion angle in degrees.
            n_divisions: Number of toroidal element divisions in the extrusion.
            coordinate_system: Cylindrical coordinate system used for local coordinates.

        Returns:
            np.ndarray: Integer array with shape (num_input, 1 + num_matches), where
                column 0 stores source entity IDs and the remaining columns store target IDs.

        Raises:
            ValueError: If inputs are invalid or the input mesh is not 2D.
        """
        ensure_vtk_unstructured_grid(input_mesh, arg_name="input_mesh")
        ensure_vtk_unstructured_grid(target_mesh, arg_name="target_mesh")
        if location not in ["node", "elem"]:
            raise ValueError("location must be 'node' or 'elem'.")
        if tolerance <= 0:
            raise ValueError("tolerance must be positive.")
        if extrusion_angle <= 0:
            raise ValueError("extrusion_angle must be positive.")

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
                    (num_matches,), candidate_ids[0], dtype=candidate_ids.dtype
                )
            else:
                trimmed = np.empty((num_matches,), dtype=candidate_ids.dtype)
                used = np.zeros(candidate_idx.size, dtype=bool)
                poloidal_dist = np.hypot(dr[candidate_idx], dz[candidate_idx])
                for theta_idx, theta_value in enumerate(expected_theta):
                    theta_dist = angular_distance(target_theta[candidate_idx], theta_value)
                    score = theta_dist + poloidal_dist / max(tolerance_current, np.finfo(float).eps)
                    if np.any(~used):
                        score = np.where(used, np.inf, score)
                    pick = int(np.argmin(score))
                    trimmed[theta_idx] = candidate_ids[pick]
                    used[pick] = True
            mapping[row_idx, 1 : 1 + trimmed.size] = trimmed[:num_matches].astype(
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
        ensure_vtk_unstructured_grid(mesh)
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
        if cp_array is None:
            raise ValueError(f"Array '{array_name}' not found in PointData.")
        cp_np = vtknp.vtk_to_numpy(cp_array)

        # Reserve CP ids used by the full contiguous CP family: CP_0, CP_1, ..., CP_N.
        existing_values: set[int] = set()
        cp_idx = 0
        while True:
            cp_family_array = point_data.GetArray(f"CP_{cp_idx}")
            if cp_family_array is None:
                break
            cp_family_np = vtknp.vtk_to_numpy(cp_family_array)
            used = cp_family_np[np.isfinite(cp_family_np) & (cp_family_np != 0)]
            if used.size > 0:
                existing_values.update(used.astype(int).tolist())
            cp_idx += 1

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

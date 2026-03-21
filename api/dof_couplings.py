from typing import List

import numpy as np
import vtk
import vtk.util.numpy_support as vtknp


class DofCouplings:
    @staticmethod
    def _collect_cp_arrays(mesh: vtk.vtkUnstructuredGrid) -> List[vtk.vtkDataArray]:
        point_data = mesh.GetPointData()
        if point_data is None:
            return []

        arrays: List[vtk.vtkDataArray] = []
        for i in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(i)
            if array is None:
                continue
            name = array.GetName()
            if not name or not name.startswith("CP_"):
                continue
            if array.GetNumberOfComponents() != 1:
                continue
            arrays.append(array)
        return arrays

    @staticmethod
    def build_cp_pairs(
        mesh: vtk.vtkUnstructuredGrid,
    ) -> np.ndarray:
        cp_arrays = DofCouplings._collect_cp_arrays(mesh)
        if not cp_arrays:
            return np.empty((0, 2), dtype=np.int64)

        cp_value_blocks: List[np.ndarray] = []
        node_id_blocks: List[np.ndarray] = []

        for cp_array in cp_arrays:
            cp_values = vtknp.vtk_to_numpy(cp_array)
            if cp_values.size == 0:
                continue
            cp_values = np.asarray(cp_values).reshape(-1)
            active_mask = cp_values > 0
            if not np.any(active_mask):
                continue

            cp_value_blocks.append(cp_values[active_mask].astype(np.int64, copy=False))
            node_id_blocks.append(np.nonzero(active_mask)[0].astype(np.int64, copy=False))

        if not cp_value_blocks:
            return np.empty((0, 2), dtype=np.int64)

        cp_values_all = np.concatenate(cp_value_blocks)
        node_ids_all = np.concatenate(node_id_blocks)

        # Deduplicate (cp_id, node_id) pairs while preserving deterministic ordering.
        pairs_matrix = np.column_stack((cp_values_all, node_ids_all))
        unique_pairs = np.unique(pairs_matrix, axis=0)
        order = np.lexsort((unique_pairs[:, 1], unique_pairs[:, 0]))
        unique_pairs = unique_pairs[order]

        cp_ids = unique_pairs[:, 0]
        node_ids = unique_pairs[:, 1]

        pair_blocks: List[np.ndarray] = []
        unique_cp_ids, start_indices = np.unique(cp_ids, return_index=True)
        boundaries = np.append(start_indices, len(cp_ids))

        for i, cp_id in enumerate(unique_cp_ids):
            start = boundaries[i]
            end = boundaries[i + 1]
            cp_node_ids = node_ids[start:end]
            if cp_node_ids.size == 0:
                continue

            master = int(cp_node_ids.min())
            slave_ids = cp_node_ids[cp_node_ids != master]
            if slave_ids.size == 0:
                continue

            masters = np.full(slave_ids.shape, master, dtype=np.int64)
            pair_blocks.append(np.column_stack((masters, slave_ids.astype(np.int64, copy=False))))

        if not pair_blocks:
            return np.empty((0, 2), dtype=np.int64)

        all_pairs = np.vstack(pair_blocks)
        all_pairs = np.unique(all_pairs, axis=0)
        return np.ascontiguousarray(all_pairs, dtype=np.int64)

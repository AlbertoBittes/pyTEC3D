# VTK Tree–Cotree Gauge via Midside Nodes (general 2nd‑order cells)
#
# This module builds a tree–cotree decomposition on a vtkUnstructuredGrid that
# contains second‑order elements (edges identified by their midside nodes).
#
# Key ideas:
# - Every edge is uniquely identified by its midside node (point id).
# - Constrained corner nodes are excluded from the spanning forest traversal.
# - You can also provide an initial boolean mask over midside nodes to force
#   some edges into the cotree (e.g., boundary/Dirichlet edges).
# - The result is written into two vtkIntArray attached to PointData
#   (size = number of points), with 1 at midside nodes for tree/cotree.
#
# Supported cell types:
#   - VTK_QUADRATIC_EDGE (3)
#   - VTK_QUADRATIC_TRIANGLE (6)
#   - VTK_QUADRATIC_QUAD (8)
#   - VTK_QUADRATIC_TETRA (10)
#   - VTK_QUADRATIC_WEDGE (15)
#   - VTK_QUADRATIC_PYRAMID (13)
#   - VTK_QUADRATIC_HEXAHEDRON (20)
#
# Usage example:
#   tree_arr, cotree_arr, info = build_tree_cotree_for_midside_edges(
#       ugrid,
#       constrained_point_ids={...},            # corner nodes to exclude entirely
#       forced_cotree_mask=your_numpy_bool_mask  # length = ugrid.GetNumberOfPoints()
#   )
#   ugrid.GetPointData().AddArray(tree_arr)
#   ugrid.GetPointData().AddArray(cotree_arr)
#
# Notes:
# - If forced_cotree_mask is None, only constrained nodes are excluded.
# - Edges that touch any constrained node are automatically cotree.
# - Among the remaining edges, we build a spanning forest using Union‑Find.
#
# (c) Provide as-is; tailor as needed for your mesh conventions.

import vtk
import numpy as np
import vtk.util.numpy_support as vtknp

# Edge tables for 2nd-order VTK cells
EDGE_TABLES = {
    vtk.VTK_QUADRATIC_EDGE: (
        [(0, 1)],
        [2],
    ),
    vtk.VTK_QUADRATIC_TRIANGLE: (
        [(0, 1), (1, 2), (2, 0)],
        [3, 4, 5],
    ),
    vtk.VTK_QUADRATIC_QUAD: (
        [(0, 1), (1, 2), (2, 3), (3, 0)],
        [4, 5, 6, 7],
    ),
    vtk.VTK_QUADRATIC_TETRA: (
        [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
        [4, 5, 6, 7, 8, 9],
    ),
    vtk.VTK_QUADRATIC_WEDGE: (
        [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
        [6, 7, 8, 9, 10, 11, 12, 13, 14],
    ),
    vtk.VTK_QUADRATIC_PYRAMID: (
        [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
        [5, 6, 7, 8, 9, 10, 11, 12],
    ),
    vtk.VTK_QUADRATIC_HEXAHEDRON: (
        [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ],
        [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
    ),
}


class DSU:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True


# ------------------------------
# Core builder
# ------------------------------
class TreeGaugeHelper:
    @staticmethod
    def build_tree_cotree_for_midside_edges(
        ugrid,
        constrained_point_ids=None,
        forced_cotree_mask=None,
        tree_name="TREE",
        cotree_name="COTREE",
    ):
        """
        Build a tree–cotree decomposition for second‑order cells, identifying edges by
        their midside node ids. Results are vtkIntArray on PointData (1 at midside nodes).

        Parameters
        ----------
        ugrid : vtkUnstructuredGrid
            Mesh with second‑order elements (quadratic). Edges are implicit.
        constrained_point_ids : set[int], optional
            Corner point ids to exclude entirely from traversal (e.g., Dirichlet nodes).
        forced_cotree_mask : sequence of bool/int, optional
            Length = number of points. Nonzero/True at midside points forces that edge into cotree.
        tree_name, cotree_name : str
            Names for the output arrays (PointData).

        Returns
        -------
        tree_arr : vtkIntArray  (PointData)
        cotree_arr : vtkIntArray  (PointData)
        info : dict
            Diagnostic info: { 'num_edges', 'num_tree', 'num_cotree', 'skipped_unknown_cells', ... }
        """
        npts = ugrid.GetNumberOfPoints()
        constrained = set(constrained_point_ids or [])

        # Prepare forced cotree mask at point level
        if forced_cotree_mask is None:
            forced_cotree_mask = np.zeros(npts, dtype=np.int32)
        elif forced_cotree_mask.shape[0] != npts:
            raise ValueError(
                "forced_cotree_mask must have length equal to number of points"
            )

        # Use np.nan for non-midside nodes
        tree_np = np.full(npts, np.nan, dtype=np.float32)
        cotree_np = np.full(npts, np.nan, dtype=np.float32)

        # Collect unique edges using midside node id as key
        edge_map = {}
        unknown_cells = 0

        for ci in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(ci)
            ctype = cell.GetCellType()
            if ctype not in EDGE_TABLES:
                unknown_cells += 1
                continue
            corner_pairs, mid_locals = EDGE_TABLES[ctype]
            ids = cell.GetPointIds()
            for (v0_local, v1_local), mid_local in zip(corner_pairs, mid_locals):
                v0 = ids.GetId(v0_local)
                v1 = ids.GetId(v1_local)
                mid = ids.GetId(mid_local)
                if mid not in edge_map:
                    edge_map[mid] = (v0, v1)
                else:
                    ev0, ev1 = edge_map[mid]
                    if not ((ev0 == v0 and ev1 == v1) or (ev0 == v1 and ev1 == v0)):
                        raise RuntimeError(
                            f"Inconsistent edge definition for midside {mid}: "
                            f"existing ({ev0},{ev1}) vs new ({v0},{v1})"
                        )

        # Forced cotree
        cotree_np[forced_cotree_mask == 1] = 1.0

        if len(constrained) > 0:
            # Constrained edges → cotree
            for mid, (a, b) in edge_map.items():
                if cotree_np[mid] == 1.0:
                    continue
                if (a in constrained) or (b in constrained):
                    cotree_np[mid] = 1.0

        # Union–Find spanning tree
        dsu = DSU()
        for mid, (a, b) in edge_map.items():
            if cotree_np[mid] == 1.0:
                continue
            if (a in constrained) or (b in constrained):
                continue
            if dsu.union(a, b):
                tree_np[mid] = 1.0
            else:
                cotree_np[mid] = 1.0

        tree_arr = vtknp.numpy_to_vtk(tree_np, deep=1, array_type=vtk.VTK_FLOAT)
        tree_arr.SetName(tree_name)
        cotree_arr = vtknp.numpy_to_vtk(cotree_np, deep=1, array_type=vtk.VTK_FLOAT)
        cotree_arr.SetName(cotree_name)

        info = {
            "num_points": npts,
            "num_cells": ugrid.GetNumberOfCells(),
            "num_edges": len(edge_map),
            "num_tree": int(tree_np.sum()),
            "num_cotree": int(cotree_np.sum()),
            "skipped_unknown_cells": unknown_cells,
        }
        return tree_arr, cotree_arr, info

    # ------------------------------
    # Convenience: attach arrays to the grid
    # ------------------------------
    @staticmethod
    def build_tree_cotree_on_grid(
        ugrid,
        constrained_point_ids=None,
        forced_cotree_mask=None,
        tree_name="TREE",
        cotree_name="COTREE",
    ):
        tree_arr, cotree_arr, info = (
            TreeGaugeHelper.build_tree_cotree_for_midside_edges(
                ugrid,
                constrained_point_ids=constrained_point_ids,
                forced_cotree_mask=forced_cotree_mask,
                tree_name=tree_name,
                cotree_name=cotree_name,
            )
        )
        ugrid.GetPointData().AddArray(tree_arr)
        ugrid.GetPointData().AddArray(cotree_arr)
        return info

    # ------------------------------
    # Optional: helper to create a forced_cotree_mask from a set of midside ids
    # ------------------------------
    @staticmethod
    def forced_mask_from_midside_set(n_points, midside_ids):
        mask = np.zeros(n_points, dtype=np.int32)
        for mid in midside_ids:
            if 0 <= mid < n_points:
                mask[mid] = 1
        return mask

    @staticmethod
    def add_midside_nodes_mask_array(ugrid, array_name="MidsideNodes"):
        """
        Adds an array to the unstructured grid with 1 on midside nodes, 0 on corner nodes.
        The mask is determined by the cell types and VTK quadratic element conventions.
        """
        npts = ugrid.GetNumberOfPoints()
        mask = np.zeros(npts, dtype=np.int32)

        # Scan all cells and mark midside nodes
        for ci in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(ci)
            ctype = cell.GetCellType()
            if ctype not in EDGE_TABLES:
                continue
            _, midside_locals = EDGE_TABLES[ctype]
            ids = cell.GetPointIds()
            for mid_local in midside_locals:
                mid = ids.GetId(mid_local)
                mask[mid] = 1

        arr = vtknp.numpy_to_vtk(mask, deep=1, array_type=vtk.VTK_INT)
        arr.SetName(array_name)
        ugrid.GetPointData().AddArray(arr)
        return arr

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

from typing import Dict, Iterable, Optional, Sequence, Set, Tuple, TYPE_CHECKING

import numpy as np
import vtk
import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

from .logging_utils import get_logger
from .validation_utils import ensure_vtk_unstructured_grid


LOGGER = get_logger("apec_m.tree_gauge")

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import vtkmodules.vtkCommonDataModel as vtk_dm  # type: ignore[import-not-found]

# Edge tables for 2nd-order VTK cells
EDGE_TABLES: Dict[int, Tuple[Sequence[Tuple[int, int]], Sequence[int]]] = {
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


class _DSU:
    """Disjoint-set union structure used to build a spanning forest."""
    def __init__(self) -> None:
        """Initialize empty parent/rank maps.

        Returns:
            None
        """
        self.parent: Dict[int, int] = {}
        self.rank: Dict[int, int] = {}

    def find(self, x: int) -> int:
        """Find representative with path compression.

        Args:
            x: Node id.

        Returns:
            int: Representative id.
        """
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0
            return x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> bool:
        """Union by rank.

        Args:
            a: First node id.
            b: Second node id.

        Returns:
            bool: True when the sets were merged, False when already connected.
        """
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
class TreeGauge:
    @staticmethod
    def build_tree_cotree_for_midside_edges(
        ugrid: "vtk.vtkUnstructuredGrid",
        constrained_point_ids: Optional[Iterable[int]] = None,
        forced_cotree_mask: Optional[np.ndarray] = None,
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Tuple["vtk.vtkIntArray", "vtk.vtkIntArray", Dict[str, int]]:
        """
        Build a tree–cotree decomposition for second‑order cells, identifying edges by
        their midside node ids. Results are vtkIntArray on PointData (1 at midside nodes).

        Args:
            ugrid: Mesh with second‑order elements (quadratic). Edges are implicit.
            constrained_point_ids: Corner point ids to exclude entirely from traversal (e.g., Dirichlet nodes).
            forced_cotree_mask: Length = number of points. Nonzero/True at midside points forces that edge into cotree.
            tree_name: Name for the tree output array (PointData).
            cotree_name: Name for the cotree output array (PointData).

        Returns:
            Tuple[vtkIntArray, vtkIntArray, Dict[str, int]]: (tree_arr, cotree_arr, info) where info contains
            diagnostic counts like num_edges, num_tree, num_cotree, skipped_unknown_cells.
        """
        ensure_vtk_unstructured_grid(ugrid, arg_name="ugrid")
        npts = ugrid.GetNumberOfPoints()
        constrained = set(constrained_point_ids or [])

        # Prepare forced cotree mask at point level
        if forced_cotree_mask is None:
            forced_cotree_mask = np.zeros(npts, dtype=np.int32)
        elif forced_cotree_mask.shape[0] != npts:
            raise ValueError(
                "forced_cotree_mask must have length equal to number of points"
            )

        # Use 0 for non-midside nodes, 1 for tree/cotree
        tree_np = np.zeros(npts, dtype=np.int32)
        cotree_np = np.zeros(npts, dtype=np.int32)

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
        cotree_np[forced_cotree_mask == 1] = 1

        if len(constrained) > 0:
            # Constrained edges → cotree
            for mid, (a, b) in edge_map.items():
                if cotree_np[mid] == 1:
                    continue
                if (a in constrained) or (b in constrained):
                    cotree_np[mid] = 1

        # Union–Find spanning tree
        dsu = _DSU()
        for mid, (a, b) in edge_map.items():
            if cotree_np[mid] == 1:
                continue
            if (a in constrained) or (b in constrained):
                continue
            if dsu.union(a, b):
                tree_np[mid] = 1
            else:
                cotree_np[mid] = 1

        tree_arr = vtknp.numpy_to_vtk(tree_np, deep=1, array_type=vtk.VTK_INT)
        tree_arr.SetName(tree_name)
        cotree_arr = vtknp.numpy_to_vtk(cotree_np, deep=1, array_type=vtk.VTK_INT)
        cotree_arr.SetName(cotree_name)

        info = {
            "num_points": npts,
            "num_cells": ugrid.GetNumberOfCells(),
            "num_edges": len(edge_map),
            "num_tree": int(np.sum(tree_np)),
            "num_cotree": int(np.sum(cotree_np)),
            "skipped_unknown_cells": unknown_cells,
        }
        return tree_arr, cotree_arr, info

    @staticmethod
    def build_tree_cotree_on_grid(
        ugrid: "vtk.vtkUnstructuredGrid",
        constrained_point_ids: Optional[Iterable[int]] = None,
        forced_cotree_mask: Optional[np.ndarray] = None,
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Dict[str, int]:
        """
        Convenience wrapper that builds tree/cotree arrays and attaches them to the grid.

        Args:
            ugrid: Mesh to annotate.
            constrained_point_ids: Corner point ids to exclude entirely from traversal.
            forced_cotree_mask: Length = number of points; truthy at midside points that must belong to the cotree.
            tree_name: Name for the tree output array.
            cotree_name: Name for the cotree output array.

        Returns:
            Dict[str, int]: Diagnostic info from `build_tree_cotree_for_midside_edges`.
        """
        tree_arr, cotree_arr, info = TreeGauge.build_tree_cotree_for_midside_edges(
            ugrid,
            constrained_point_ids=constrained_point_ids,
            forced_cotree_mask=forced_cotree_mask,
            tree_name=tree_name,
            cotree_name=cotree_name,
        )
        ugrid.GetPointData().AddArray(tree_arr)
        ugrid.GetPointData().AddArray(cotree_arr)
        return info

    @staticmethod
    def tree_gauge(
        ugrid: "vtk.vtkUnstructuredGrid",
        array_name: str = "Az",
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Dict[str, int]:
        """
        Build tree/cotree masks and zero out the tree-marked entries in a PointData array.

        This mirrors the APEC_M.tree_gauge behavior but operates directly on the mesh.

        Args:
            ugrid: Target mesh to analyze and annotate.
            array_name: PointData array used to derive constrained midside edges.
            tree_name: Name for the output tree mask array.
            cotree_name: Name for the output cotree mask array.

        Returns:
            Dict[str, int]: Diagnostic info from the tree/cotree build.
        """
        arr = ugrid.GetPointData().GetArray(array_name)
        if arr is None:
            raise ValueError(f"Array '{array_name}' not found in PointData.")
        arr_np = vtknp.vtk_to_numpy(arr)

        constrained_corners = TreeGauge.constrained_corners_from_midside_array(
            ugrid, array_name
        )

        info = TreeGauge.build_tree_cotree_on_grid(
            ugrid,
            constrained_point_ids=constrained_corners,
            tree_name=tree_name,
            cotree_name=cotree_name,
        )
        LOGGER.info("Tree/cotree gauge built. Info: %s", info)

        tree_arr = ugrid.GetPointData().GetArray(tree_name)
        if tree_arr is not None:
            tree_np = vtknp.vtk_to_numpy(tree_arr)
            arr_np[tree_np == 1] = 0
            arr.Modified()

        return info

    # ------------------------------
    # Optional: helper to create a forced_cotree_mask from a set of midside ids
    # ------------------------------
    @staticmethod
    def forced_mask_from_midside_set(
        n_points: int, midside_ids: Iterable[int]
    ) -> np.ndarray:
        """
        Create a cotree mask (length = number of points) with 1 at specified midside ids.

        Args:
            n_points: Total number of points in the mesh.
            midside_ids: Iterable of midside point ids to force into the cotree.

        Returns:
            np.ndarray: Integer mask with ones at midside ids, zeros elsewhere.
        """
        mask = np.zeros(n_points, dtype=np.int32)
        for mid in midside_ids:
            if 0 <= mid < n_points:
                mask[mid] = 1
        return mask

    @staticmethod
    def add_midside_nodes_mask_array(
        ugrid: "vtk.vtkUnstructuredGrid", array_name: str = "MidsideNodes"
    ) -> "vtk.vtkIntArray":
        """
        Add a PointData array with 1 on midside nodes, 0 on corner nodes.

        The mask is determined by the cell types and VTK quadratic element conventions.

        Args:
            ugrid: Target unstructured grid.
            array_name: Name for the mask array (default "MidsideNodes").

        Returns:
            vtkIntArray: Added mask array.
        """
        ensure_vtk_unstructured_grid(ugrid, arg_name="ugrid")
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

    @staticmethod
    def constrained_corners_from_midside_array(
        ugrid: "vtk.vtkUnstructuredGrid", array_name: str
    ) -> Set[int]:
        """
        Given an array in PointData (array_name) with NaN at unconstrained midside nodes or corner nodes,
        and a value at constrained midside nodes, return a set of corner node ids to be used as constrained_point_ids
        for build_tree_cotree_for_midside_edges.

        The set will include:
        - All corner nodes connected to any constrained midside node.

        Args:
            ugrid: Unstructured grid containing quadratic elements.
            array_name: Name of the PointData array with constraints on midside nodes.

        Returns:
            Set[int]: Corner point ids connected to constrained midside nodes.
        """
        ensure_vtk_unstructured_grid(ugrid, arg_name="ugrid")
        arr = ugrid.GetPointData().GetArray(array_name)
        if arr is None:
            raise ValueError(f"Array '{array_name}' not found in PointData.")
        arr_np = vtknp.vtk_to_numpy(arr)

        midside_mask_arr = ugrid.GetPointData().GetArray("MidsideNodes")
        if midside_mask_arr is None:
            raise ValueError("Mesh missing 'MidsideNodes' array.")
        midside_mask_np = vtknp.vtk_to_numpy(midside_mask_arr)

        # Find constrained midside node ids
        constrained_midside_ids = np.where(
            (midside_mask_np == 1) & (~np.isnan(arr_np))
        )[0]

        constrained_set = set()

        # For each constrained midside node, add its connected corner nodes
        for ci in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(ci)
            ctype = cell.GetCellType()
            if ctype not in EDGE_TABLES:
                continue
            corner_pairs, midside_locals = EDGE_TABLES[ctype]
            ids = cell.GetPointIds()
            for (v0_local, v1_local), mid_local in zip(corner_pairs, midside_locals):
                mid = ids.GetId(mid_local)
                if mid in constrained_midside_ids:
                    v0 = ids.GetId(v0_local)
                    v1 = ids.GetId(v1_local)
                    constrained_set.add(v0)
                    constrained_set.add(v1)

        return constrained_set

# VTK Tree–Cotree Gauge via Midside Nodes (general 2nd‑order cells)
#
# This module builds a tree–cotree decomposition on a vtkUnstructuredGrid that
# contains second‑order elements (edges identified by their midside nodes).
#
# Key ideas:
# - Every edge is uniquely identified by its midside node (point id).
# - The tree is built with boundary-first BFS over non-conducting nodes.
# - Start-surface nodes provide BFS seeds (one or more per disconnected region).
# - Periodic slave/master midside couplings are collapsed to canonical edges.
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
#       non_conducting_point_ids={...},
#       start_surface_point_ids={...},
#       periodic_midside_coupling=periodic_map,
#   )
#   ugrid.GetPointData().AddArray(tree_arr)
#   ugrid.GetPointData().AddArray(cotree_arr)
#
# Notes:
# - Among the remaining eligible edges, a boundary-first BFS builds the tree.
#
# (c) Provide as-is; tailor as needed for your mesh conventions.

from collections import deque
from typing import Dict, Iterable, Mapping, Optional, Sequence, Set, Tuple, TYPE_CHECKING

import numpy as np
import vtk
import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

from .logging_utils import get_logger
from .selector import Selector
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


# ------------------------------
# Core builder
# ------------------------------
class TreeGauge:
    @staticmethod
    def _normalize_periodic_midside_mapping(
        periodic_midside_coupling: Optional[object],
        ugrid: Optional["vtk.vtkUnstructuredGrid"] = None,
        midside_mask: Optional[np.ndarray] = None,
    ) -> Dict[int, int]:
        """Normalize periodic midside coupling as slave->master map."""
        if periodic_midside_coupling is None:
            return {}

        if isinstance(periodic_midside_coupling, str):
            cp_array_name = periodic_midside_coupling
            if not cp_array_name.startswith("CP_"):
                raise ValueError(
                    "periodic_midside_coupling string must be a CP_* PointData array name."
                )
            if ugrid is None or midside_mask is None:
                raise ValueError(
                    "ugrid and midside_mask are required when periodic_midside_coupling is a string."
                )

            cp_arr = ugrid.GetPointData().GetArray(cp_array_name)
            if cp_arr is None:
                raise ValueError(
                    f"PointData array '{cp_array_name}' was not found on the mesh."
                )
            cp_np = vtknp.vtk_to_numpy(cp_arr)
            if cp_np.shape[0] != ugrid.GetNumberOfPoints():
                raise ValueError(
                    f"'{cp_array_name}' length must equal number of points in the mesh."
                )

            midside_ids = np.where(midside_mask == 1)[0].astype(int)
            groups: Dict[object, list[int]] = {}
            for mid in midside_ids:
                value = cp_np[mid]
                if isinstance(value, (float, np.floating)) and np.isnan(value):
                    continue
                key = value.item() if isinstance(value, np.generic) else value
                if key == 0:
                    continue
                groups.setdefault(key, []).append(mid)

            mapping: Dict[int, int] = {}
            for mids in groups.values():
                if len(mids) < 2:
                    continue
                master = int(min(mids))
                for slave in mids:
                    slave_i = int(slave)
                    if slave_i != master:
                        mapping[slave_i] = master
            return mapping

        mapping: Dict[int, int] = {}
        if isinstance(periodic_midside_coupling, Mapping):
            for slave, master in periodic_midside_coupling.items():
                mapping[int(slave)] = int(master)
            return mapping

        if isinstance(periodic_midside_coupling, np.ndarray):
            if periodic_midside_coupling.ndim != 2 or periodic_midside_coupling.shape[1] != 2:
                raise ValueError(
                    "periodic_midside_coupling ndarray must have shape (N, 2)"
                )
            for row in periodic_midside_coupling:
                mapping[int(row[0])] = int(row[1])
            return mapping

        raise TypeError(
            "periodic_midside_coupling must be dict-like, a (N, 2) numpy array, a CP_* array name string, or None"
        )

    @staticmethod
    def _canonical_midside_id(mid: int, slave_to_master: Dict[int, int]) -> int:
        """Follow slave->master links until a canonical midside id is reached."""
        seen: Set[int] = set()
        cur = int(mid)
        while cur in slave_to_master:
            if cur in seen:
                raise ValueError(
                    f"Cycle detected in periodic_midside_coupling at midside id {cur}."
                )
            seen.add(cur)
            nxt = int(slave_to_master[cur])
            if nxt == cur:
                break
            cur = nxt
        return cur

    @staticmethod
    def _canonical_corner_id(corner: int, slave_to_master: Dict[int, int]) -> int:
        """Follow slave->master links until a canonical corner id is reached."""
        seen: Set[int] = set()
        cur = int(corner)
        while cur in slave_to_master:
            if cur in seen:
                raise ValueError(
                    f"Cycle detected in periodic corner mapping at point id {cur}."
                )
            seen.add(cur)
            nxt = int(slave_to_master[cur])
            if nxt == cur:
                break
            cur = nxt
        return cur

    @staticmethod
    def _normalize_periodic_corner_mapping(
        periodic_midside_coupling: Optional[object],
        ugrid: "vtk.vtkUnstructuredGrid",
        midside_mask: np.ndarray,
        raw_edge_map: Dict[int, Tuple[int, int]],
        canonical_to_raw_mids: Dict[int, Set[int]],
    ) -> Dict[int, int]:
        """Build corner-node slave->master periodic map for topology canonicalization."""
        # Preferred path: use CP_* values on corner nodes directly.
        if isinstance(periodic_midside_coupling, str):
            cp_arr = ugrid.GetPointData().GetArray(periodic_midside_coupling)
            if cp_arr is not None:
                cp_np = vtknp.vtk_to_numpy(cp_arr)
                corner_ids = np.where(midside_mask == 0)[0].astype(int)
                groups: Dict[object, list[int]] = {}
                for pid in corner_ids:
                    value = cp_np[pid]
                    if isinstance(value, (float, np.floating)) and np.isnan(value):
                        continue
                    key = value.item() if isinstance(value, np.generic) else value
                    if key == 0:
                        continue
                    groups.setdefault(key, []).append(int(pid))

                mapping: Dict[int, int] = {}
                for corners in groups.values():
                    if len(corners) < 2:
                        continue
                    master = int(min(corners))
                    for slave in corners:
                        if slave != master:
                            mapping[int(slave)] = master
                return mapping

        # Fallback: infer corner correspondences from periodic-collapsed midside edge groups.
        mapping: Dict[int, int] = {}
        for raw_mids in canonical_to_raw_mids.values():
            mids = sorted(int(m) for m in raw_mids if int(m) in raw_edge_map)
            if len(mids) < 2:
                continue

            rep_raw = mids[0]
            rep_a, rep_b = raw_edge_map[rep_raw]
            for raw_mid in mids[1:]:
                a, b = raw_edge_map[raw_mid]

                options = [
                    ((a, rep_a), (b, rep_b)),
                    ((a, rep_b), (b, rep_a)),
                ]

                def conflict_count(option: Tuple[Tuple[int, int], Tuple[int, int]]) -> int:
                    conflicts = 0
                    for slave, master in option:
                        if slave == master:
                            continue
                        if slave in mapping and mapping[slave] != master:
                            conflicts += 1
                    return conflicts

                option0_conflicts = conflict_count(options[0])
                option1_conflicts = conflict_count(options[1])
                chosen = options[0] if option0_conflicts <= option1_conflicts else options[1]

                for slave, master in chosen:
                    if slave == master:
                        continue
                    if slave in mapping and mapping[slave] != master:
                        continue
                    mapping[slave] = master

        return mapping

    @staticmethod
    def build_tree_cotree_for_midside_edges(
        ugrid: "vtk.vtkUnstructuredGrid",
        non_conducting_point_ids: Optional[Iterable[int]] = None,
        start_surface_point_ids: Optional[Iterable[int]] = None,
        periodic_midside_coupling: Optional[object] = None,
        spanning_method: str = "kruskal",
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Tuple["vtk.vtkIntArray", "vtk.vtkIntArray", Dict[str, int]]:
        """
        Build a tree–cotree decomposition for second‑order cells, identifying edges by
        their midside node ids. Results are vtkIntArray on PointData (1 at midside nodes).

        Args:
            ugrid: Mesh with second‑order elements (quadratic). Edges are implicit.
            non_conducting_point_ids: Point ids for non-conducting regions where gauge tree constraints are built.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to nodes attached to currently selected elements (Selector.nsle).
            start_surface_point_ids: Point ids used as BFS roots (boundary-first traversal).
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                Multiple disconnected regions are supported by providing roots in each region.
                If omitted, defaults to currently selected nodes.
            periodic_midside_coupling: Periodic edge coupling as either a dict-like slave->master map
                or an ndarray with shape (N, 2), where each row is [slave_id, master_id],
                or a CP_* PointData array name. For CP_* arrays, midside edges sharing the
                same CP value are coupled automatically. Non-midside pairs are ignored.
            spanning_method: Tree construction algorithm. "kruskal" (default) uses
                a weighted spanning-forest with priority levels
                (Dirichlet boundary edges: 0, periodic-cut edges: 1, interior edges: 2+quality).
                "bfs" runs boundary-first BFS from start_surface_point_ids.
            tree_name: Name for the tree output array (PointData).
            cotree_name: Name for the cotree output array (PointData).

        Returns:
            Tuple[vtkIntArray, vtkIntArray, Dict[str, int]]: (tree_arr, cotree_arr, info) where info contains
            diagnostic counts like num_edges, num_tree, num_cotree, skipped_unknown_cells.
        """
        ensure_vtk_unstructured_grid(ugrid, arg_name="ugrid")
        npts = ugrid.GetNumberOfPoints()
        point_sel_array = ugrid.GetPointData().GetArray("selection_point")
        point_sel_original = vtknp.vtk_to_numpy(point_sel_array).copy()

        if start_surface_point_ids is None:
            start_surface_input: Set[int] = set(
                np.where(point_sel_original == 1)[0].astype(int).tolist()
            )
        else:
            start_surface_input = set(int(pid) for pid in start_surface_point_ids)

        if non_conducting_point_ids is None:
            # Build default non-conducting set from selected elements.
            Selector.nsle(ugrid, mode="s")
            point_sel_current = vtknp.vtk_to_numpy(point_sel_array)
            non_conducting_input: Set[int] = set(
                np.where(point_sel_current == 1)[0].astype(int).tolist()
            )
            # Preserve caller selection state.
            point_sel_current[:] = point_sel_original
            point_sel_array.Modified()
        else:
            non_conducting_input = set(int(pid) for pid in non_conducting_point_ids)

        midside_arr = ugrid.GetPointData().GetArray("MidsideNodes")
        midside_mask = vtknp.vtk_to_numpy(midside_arr)

        # BFS graph nodes are corner nodes (edge endpoints). Drop midside ids from user sets.
        raw_non_conducting = set(i for i in non_conducting_input if midside_mask[i] == 0)
        raw_start_surface = set(i for i in start_surface_input if midside_mask[i] == 0)

        periodic_raw = TreeGauge._normalize_periodic_midside_mapping(
            periodic_midside_coupling,
            ugrid=ugrid,
            midside_mask=midside_mask,
        )
        slave_to_master: Dict[int, int] = {}
        for slave, master in periodic_raw.items():
            if midside_mask[slave] != 1 or midside_mask[master] != 1:
                continue
            slave_to_master[slave] = master

        # Use 0 for non-midside nodes, 1 for tree/cotree
        tree_np = np.zeros(npts, dtype=np.int32)
        cotree_np = np.zeros(npts, dtype=np.int32)

        # Collect unique edges using canonical midside id as key.
        # The same canonical edge may have multiple raw midside ids due to periodic coupling.
        edge_map: Dict[int, Tuple[int, int]] = {}
        raw_edge_map: Dict[int, Tuple[int, int]] = {}
        canonical_to_raw_mids: Dict[int, Set[int]] = {}
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
                raw_mid = ids.GetId(mid_local)
                if raw_mid not in raw_edge_map:
                    raw_edge_map[raw_mid] = (v0, v1)
                mid = TreeGauge._canonical_midside_id(raw_mid, slave_to_master)

                if mid not in canonical_to_raw_mids:
                    canonical_to_raw_mids[mid] = set()
                canonical_to_raw_mids[mid].add(raw_mid)

                if mid not in edge_map:
                    edge_map[mid] = (v0, v1)
                else:
                    # For periodic-coupled midsides, corner ids can differ across
                    # opposite surfaces. Keep the first representative endpoints.
                    pass

        corner_slave_to_master = TreeGauge._normalize_periodic_corner_mapping(
            periodic_midside_coupling,
            ugrid,
            midside_mask,
            raw_edge_map,
            canonical_to_raw_mids,
        )

        # Canonicalize corner-node ids before topology/BFS construction.
        non_conducting = {
            TreeGauge._canonical_corner_id(pid, corner_slave_to_master)
            for pid in raw_non_conducting
        }
        start_surface = {
            TreeGauge._canonical_corner_id(pid, corner_slave_to_master)
            for pid in raw_start_surface
        }
        for mid, (a, b) in list(edge_map.items()):
            edge_map[mid] = (
                TreeGauge._canonical_corner_id(a, corner_slave_to_master),
                TreeGauge._canonical_corner_id(b, corner_slave_to_master),
            )

        # Eligible edges are those where constraints can be computed.
        eligible_edges: Dict[int, Tuple[int, int]] = {}
        for mid, (a, b) in edge_map.items():
            if (a not in non_conducting) or (b not in non_conducting):
                # In conducting regions no tree-gauge constraints are needed.
                for raw_mid in canonical_to_raw_mids[mid]:
                    cotree_np[raw_mid] = 1
                continue
            eligible_edges[mid] = (a, b)

        num_edges_before_periodic_collapse = sum(
            1
            for a, b in raw_edge_map.values()
            if (a in raw_non_conducting) and (b in raw_non_conducting)
        )

        # Build adjacency for boundary-first BFS.
        adjacency: Dict[int, list[Tuple[int, int]]] = {}
        for mid, (a, b) in eligible_edges.items():
            adjacency.setdefault(a, []).append((b, mid))
            adjacency.setdefault(b, []).append((a, mid))

        # Z2 seam marker: canonical edges created by collapsing multiple periodic raw midsides.
        periodic_edge_parity: Dict[int, int] = {
            mid: 1 if len(canonical_to_raw_mids.get(mid, set())) > 1 else 0
            for mid in eligible_edges
        }

        method = spanning_method.strip().lower()
        if method not in {"kruskal", "bfs"}:
            raise ValueError("spanning_method must be either 'kruskal' or 'bfs'.")

        def run_bfs() -> Dict[int, int]:
            edge_state_local: Dict[int, int] = {}  # 1=tree, 2=cotree
            visited: Set[int] = set()
            queue: deque[int] = deque()

            def run_bfs_from_seed(seed: int) -> None:
                if seed in visited:
                    return
                if seed not in adjacency:
                    visited.add(seed)
                    return

                visited.add(seed)
                queue.append(seed)
                while queue:
                    node = queue.popleft()
                    for nbr, mid in adjacency.get(node, []):
                        if mid in edge_state_local:
                            continue
                        if nbr not in visited:
                            edge_state_local[mid] = 1
                            visited.add(nbr)
                            queue.append(nbr)
                        else:
                            edge_state_local[mid] = 2

            # Primary boundary-first traversal from user-provided start surfaces.
            for seed in sorted(start_surface):
                if seed in non_conducting:
                    run_bfs_from_seed(seed)

            for node in sorted(adjacency.keys()):
                if node not in visited:
                    run_bfs_from_seed(node)

            # Any remaining eligible edge not classified by BFS is conservatively cotree.
            for mid in eligible_edges:
                if mid not in edge_state_local:
                    edge_state_local[mid] = 2

            return edge_state_local

        def run_kruskal() -> Dict[int, int]:
            edge_state_local: Dict[int, int] = {}  # 1=tree, 2=cotree
            # Kruskal spanning forest with multilevel priorities:
            # 0: Dirichlet boundary edges, 1: periodic-cut edges, 2+quality: interior edges.
            nodes = sorted(adjacency.keys())
            parent_uf: Dict[int, int] = {n: n for n in nodes}
            rank_uf: Dict[int, int] = {n: 0 for n in nodes}

            def uf_find(x: int) -> int:
                root = x
                while parent_uf[root] != root:
                    root = parent_uf[root]
                while parent_uf[x] != x:
                    nxt = parent_uf[x]
                    parent_uf[x] = root
                    x = nxt
                return root

            def uf_union(a: int, b: int) -> bool:
                ra = uf_find(a)
                rb = uf_find(b)
                if ra == rb:
                    return False
                if rank_uf[ra] < rank_uf[rb]:
                    parent_uf[ra] = rb
                elif rank_uf[ra] > rank_uf[rb]:
                    parent_uf[rb] = ra
                else:
                    parent_uf[rb] = ra
                    rank_uf[ra] += 1
                return True

            def edge_weight(mid: int, a: int, b: int) -> tuple[int, float, int]:
                if (a in start_surface) and (b in start_surface):
                    return (0, 0.0, mid)
                if periodic_edge_parity.get(mid, 0) == 1:
                    return (1, 0.0, mid)
                pa = np.asarray(ugrid.GetPoint(a), dtype=float)
                pb = np.asarray(ugrid.GetPoint(b), dtype=float)
                quality = float(np.linalg.norm(pa - pb))
                return (2, quality, mid)

            weighted_edges: list[tuple[tuple[int, float, int], int, int, int]] = []
            for mid, (a, b) in eligible_edges.items():
                weighted_edges.append((edge_weight(mid, a, b), mid, a, b))
            weighted_edges.sort(key=lambda item: item[0])

            for _, mid, a, b in weighted_edges:
                if uf_union(a, b):
                    edge_state_local[mid] = 1
                else:
                    edge_state_local[mid] = 2

            return edge_state_local

        edge_state = run_bfs() if method == "bfs" else run_kruskal()

        # Compute parity-to-root and component labels from the current tree set.
        parity_to_root: Dict[int, int] = {}
        node_component: Dict[int, int] = {}
        tree_adjacency: Dict[int, list[Tuple[int, int]]] = {}
        for mid, state in edge_state.items():
            if state != 1:
                continue
            a, b = eligible_edges[mid]
            tree_adjacency.setdefault(a, []).append((b, mid))
            tree_adjacency.setdefault(b, []).append((a, mid))

        component_id = 0
        seen_nodes: Set[int] = set()
        for root in sorted(adjacency.keys()):
            if root in seen_nodes:
                continue
            stack: list[int] = [root]
            seen_nodes.add(root)
            node_component[root] = component_id
            parity_to_root[root] = 0
            while stack:
                node = stack.pop()
                for nbr, mid in tree_adjacency.get(node, []):
                    if nbr in seen_nodes:
                        continue
                    seen_nodes.add(nbr)
                    node_component[nbr] = component_id
                    parity_to_root[nbr] = parity_to_root[node] ^ periodic_edge_parity.get(mid, 0)
                    stack.append(nbr)
            component_id += 1

        # Kill remaining periodic global modes by anchoring one odd-periodic
        # fundamental-cycle cotree edge per connected component.
        odd_periodic_cotree: list[int] = []
        odd_by_component: Dict[int, list[int]] = {}
        for mid, state in edge_state.items():
            if state != 2:
                continue
            a, b = eligible_edges[mid]
            parity = parity_to_root.get(a, 0) ^ parity_to_root.get(b, 0) ^ periodic_edge_parity.get(mid, 0)
            if parity == 1:
                odd_periodic_cotree.append(mid)
                comp = node_component.get(a, -1)
                odd_by_component.setdefault(comp, []).append(mid)

        anchored_global_mode_mids: list[int] = []
        for mids in odd_by_component.values():
            if not mids:
                continue
            anchored_mid = min(mids)
            edge_state[anchored_mid] = 1
            anchored_global_mode_mids.append(anchored_mid)

        anchored_global_mode_mid: Optional[int] = (
            min(anchored_global_mode_mids) if anchored_global_mode_mids else None
        )

        # Push canonical decisions back to the raw midside ids.
        for mid, state in edge_state.items():
            for raw_mid in canonical_to_raw_mids[mid]:
                # Stop the tree constraint from hitting the dependent DoF
                if raw_mid in slave_to_master:
                    continue

                if state == 1:
                    tree_np[raw_mid] = 1
                else:
                    cotree_np[raw_mid] = 1

        # Keep sets disjoint.
        tree_np[cotree_np == 1] = 0

        tree_arr = vtknp.numpy_to_vtk(tree_np, deep=1, array_type=vtk.VTK_INT)
        tree_arr.SetName(tree_name)
        cotree_arr = vtknp.numpy_to_vtk(cotree_np, deep=1, array_type=vtk.VTK_INT)
        cotree_arr.SetName(cotree_name)

        num_edges_after_periodic_collapse = len(eligible_edges)

        info = {
            "num_points": npts,
            "num_cells": ugrid.GetNumberOfCells(),
            "num_edges": num_edges_before_periodic_collapse,
            "num_edges_after_periodic_collapse": num_edges_after_periodic_collapse,
            "num_tree": int(np.sum(tree_np)),
            "num_cotree": int(np.sum(cotree_np)),
            "skipped_unknown_cells": unknown_cells,
            "spanning_method": method,
            "num_periodic_couplings": len(slave_to_master),
            "num_odd_periodic_fundamental_cycles": len(odd_periodic_cotree),
            "num_odd_periodic_components": len(odd_by_component),
            "num_global_mode_anchors": len(anchored_global_mode_mids),
            "anchored_global_mode_mids": anchored_global_mode_mids,
            "anchored_global_mode_mid": -1 if anchored_global_mode_mid is None else int(anchored_global_mode_mid),
        }
        return tree_arr, cotree_arr, info

    @staticmethod
    def build_tree_cotree_on_grid(
        ugrid: "vtk.vtkUnstructuredGrid",
        non_conducting_point_ids: Optional[Iterable[int]] = None,
        start_surface_point_ids: Optional[Iterable[int]] = None,
        periodic_midside_coupling: Optional[object] = None,
        spanning_method: str = "kruskal",
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Dict[str, int]:
        """
        Convenience wrapper that builds tree/cotree arrays and attaches them to the grid.

        Args:
            ugrid: Mesh to annotate.
            non_conducting_point_ids: Point ids where tree constraints are computed.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to nodes attached to currently selected elements (Selector.nsle).
            start_surface_point_ids: Point ids used as BFS roots.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to currently selected nodes.
            periodic_midside_coupling: Periodic midside edge coupling (slave->master),
                or CP_* PointData array name.
            spanning_method: "kruskal" (default) or "bfs".
            tree_name: Name for the tree output array.
            cotree_name: Name for the cotree output array.

        Returns:
            Dict[str, int]: Diagnostic info from `build_tree_cotree_for_midside_edges`.
        """
        tree_arr, cotree_arr, info = TreeGauge.build_tree_cotree_for_midside_edges(
            ugrid,
            non_conducting_point_ids=non_conducting_point_ids,
            start_surface_point_ids=start_surface_point_ids,
            periodic_midside_coupling=periodic_midside_coupling,
            spanning_method=spanning_method,
            tree_name=tree_name,
            cotree_name=cotree_name,
        )
        # Replace existing arrays with the same names to avoid stale duplicates.
        ugrid.GetPointData().RemoveArray(tree_name)
        ugrid.GetPointData().RemoveArray(cotree_name)
        ugrid.GetPointData().AddArray(tree_arr)
        ugrid.GetPointData().AddArray(cotree_arr)
        return info

    @staticmethod
    def tree_gauge(
        ugrid: "vtk.vtkUnstructuredGrid",
        non_conducting_point_ids: Optional[Iterable[int]] = None,
        start_surface_point_ids: Optional[Iterable[int]] = None,
        periodic_midside_coupling: Optional[object] = None,
        spanning_method: str = "kruskal",
        array_name: str = "_Az",
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
    ) -> Dict[str, int]:
        """
        Build tree/cotree masks and zero out the tree-marked entries in a PointData array.

        Args:
            ugrid: Target mesh to analyze and annotate.
            array_name: PointData array to be gauge-fixed (tree entries set to zero).
            tree_name: Name for the output tree mask array.
            cotree_name: Name for the output cotree mask array.
            non_conducting_point_ids: Point ids where tree constraints are computed.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to nodes attached to currently selected elements (Selector.nsle).
            start_surface_point_ids: BFS root point ids.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to currently selected nodes.
            periodic_midside_coupling: Periodic midside edge coupling (slave->master)
                as dict-like, ndarray (N,2), or CP_* PointData array name.
            spanning_method: "kruskal" (default) or "bfs".

        Returns:
            Dict[str, int]: Diagnostic info from the tree/cotree build.
        """
        arr = ugrid.GetPointData().GetArray(array_name)
        if arr is None:
            raise ValueError(f"Array '{array_name}' not found in PointData.")
        arr_np = vtknp.vtk_to_numpy(arr)

        non_conducting_points = (
            None
            if non_conducting_point_ids is None
            else set(int(i) for i in non_conducting_point_ids)
        )
        start_surface_points = (
            None
            if start_surface_point_ids is None
            else set(int(i) for i in start_surface_point_ids)
        )

        LOGGER.info(
            "Tree gauge input: non-conducting points=%d, start-surface points=%d",
            -1 if non_conducting_points is None else len(non_conducting_points),
            -1 if start_surface_points is None else len(start_surface_points),
        )

        info = TreeGauge.build_tree_cotree_on_grid(
            ugrid,
            non_conducting_point_ids=non_conducting_points,
            start_surface_point_ids=start_surface_points,
            periodic_midside_coupling=periodic_midside_coupling,
            spanning_method=spanning_method,
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

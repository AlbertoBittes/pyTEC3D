# VTK Tree–Cotree Gauge via Midside Nodes (general 2nd-order cells)
#
# Changes in this revision (relative to the previous fixed version):
#
#   [CHG-1] Removed `match_periodic_points_by_geometry` — geometric matching
#           of periodic node pairs is now handled by a separate module owned
#           by the caller. This module only consumes the resulting mapping.
#
#   [CHG-2] `periodic_midside_coupling` and `periodic_corner_coupling` are
#           replaced by a single `periodic_point_coupling` argument. The
#           caller now supplies ONE dict / (N,2) ndarray / CP_* array name
#           containing BOTH corner and midside slave->master pairs mixed
#           together; this module splits them internally by point type
#           (using the MidsideNodes mask) into the midside-edge map and the
#           corner-node map it needs internally.
#
#   [CHG-3] `tree_gauge()` no longer zeros out `array_name` entries at
#           periodic-slave points, even if they carry the tree bit. Those
#           points are tied to their master (A_slave = ± A_master, applied
#           elsewhere by the caller's periodic constraint code) rather than
#           independently fixed to zero. Non-slave tree points are zeroed as
#           before.
#
# Carried over from the previous revision:
#
#   [FIX-3] BFS root selection is true multi-source: ALL start-surface seeds
#           belonging to a connected component act as simultaneous roots.
#           This means you can freely pass point ids from several different
#           physical surfaces (e.g. multiple PEC boundaries bounding one air
#           region) in a single `start_surface_point_ids` set — each will
#           anchor its own local branch of the tree rather than only the
#           first-encountered one winning, which keeps tree paths short.
#
#   [FIX-4] Diagnostics: `info` reports connected-component count and
#           validates the spanning-forest identity
#           (#tree_edges == #nodes - #components); mismatches are logged as
#           errors so a broken graph is never silent.
#
#   [FIX-5] Skipped (non-quadratic / unknown) cells emit a logged warning
#           with the offending cell type ids instead of a silent count only.
#
# Known limitation (unchanged, documented on purpose): if the non-conducting
# region is not simply connected (e.g. air threading a conductor loop), plain
# tree-cotree is still correct and sufficient for an A / A-V formulation
# driven by a prescribed (divergence-free) source current — the source term
# supplies the circulations that a scalar-potential formulation would need
# explicit cohomology "cuts" for. Cuts are NOT implemented here and are not
# needed unless you introduce a magnetic scalar potential (T-Omega / reduced
# scalar potential) somewhere in the formulation.
#
# Supported cell types:
#   - VTK_QUADRATIC_EDGE
#   - VTK_QUADRATIC_TRIANGLE
#   - VTK_QUADRATIC_QUAD
#   - VTK_QUADRATIC_TETRA
#   - VTK_QUADRATIC_WEDGE
#   - VTK_QUADRATIC_PYRAMID
#   - VTK_QUADRATIC_HEXAHEDRON
#
# Usage example:
#   tree_arr, cotree_arr, info = build_tree_cotree_for_midside_edges(
#       ugrid,
#       non_conducting_point_ids={...},
#       start_surface_point_ids={...},          # can span multiple surfaces
#       periodic_point_coupling=periodic_map,   # corner + midside ids mixed
#   )
#   ugrid.GetPointData().AddArray(tree_arr)
#   ugrid.GetPointData().AddArray(cotree_arr)
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
    def _canonical_id(point_id: int, slave_to_master: Dict[int, int]) -> int:
        """Follow slave->master links until a canonical point id is reached.
        Shared helper for both midside and corner canonicalization."""
        seen: Set[int] = set()
        cur = int(point_id)
        while cur in slave_to_master:
            if cur in seen:
                raise ValueError(
                    f"Cycle detected in periodic point coupling at point id {cur}."
                )
            seen.add(cur)
            nxt = int(slave_to_master[cur])
            if nxt == cur:
                break
            cur = nxt
        return cur

    # [CHG-2] Unified periodic coupling normalization: accepts corner and
    # midside slave->master pairs mixed together in a single input, and
    # splits them into two separate maps by point type.
    @staticmethod
    def _normalize_periodic_point_coupling(
        periodic_point_coupling: Optional[object],
        ugrid: "vtk.vtkUnstructuredGrid",
        midside_mask: np.ndarray,
    ) -> Tuple[Dict[int, int], Dict[int, int]]:
        """
        Normalize a unified periodic point coupling input into two separate
        slave->master maps, one for midside (edge) ids and one for corner
        (node) ids.

        Accepted formats for `periodic_point_coupling`:
          - None: both maps are empty.
          - str starting with "CP_": a PointData array name. Points sharing
            the same nonzero CP value are grouped, with the lowest point id
            in each group chosen as master. Corner and midside points are
            grouped independently (a group never mixes the two types, since
            grouping is performed separately per node type using each
            point's own CP value from the same array).
          - Mapping[int, int] or ndarray of shape (N, 2): explicit
            slave->master pairs with corner and midside ids mixed together.
            Each pair is routed to the corner or midside map according to
            the type of its endpoints. A pair whose slave and master are of
            different types (one corner, one midside) is dropped with a
            logged warning, since that indicates a mis-specified mapping.

        Returns:
            Tuple[Dict[int, int], Dict[int, int]]:
                (midside_slave_to_master, corner_slave_to_master)
        """
        if periodic_point_coupling is None:
            return {}, {}

        if isinstance(periodic_point_coupling, str):
            cp_array_name = periodic_point_coupling
            if not cp_array_name.startswith("CP_"):
                raise ValueError(
                    "periodic_point_coupling string must be a CP_* PointData array name."
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

            def build_groups(mask_value: int) -> Dict[int, int]:
                ids = np.where(midside_mask == mask_value)[0].astype(int)
                groups: Dict[object, list[int]] = {}
                for pid in ids:
                    value = cp_np[pid]
                    if isinstance(value, (float, np.floating)) and np.isnan(value):
                        continue
                    key = value.item() if isinstance(value, np.generic) else value
                    if key == 0:
                        continue
                    groups.setdefault(key, []).append(int(pid))
                mapping: Dict[int, int] = {}
                for members in groups.values():
                    if len(members) < 2:
                        continue
                    master = int(min(members))
                    for slave in members:
                        if slave != master:
                            mapping[int(slave)] = master
                return mapping

            midside_map = build_groups(1)
            corner_map = build_groups(0)
            return midside_map, corner_map

        # dict / ndarray path: pairs may mix corner and midside ids together.
        raw_pairs: list[Tuple[int, int]] = []
        if isinstance(periodic_point_coupling, Mapping):
            raw_pairs = [(int(s), int(m)) for s, m in periodic_point_coupling.items()]
        elif isinstance(periodic_point_coupling, np.ndarray):
            if periodic_point_coupling.ndim != 2 or periodic_point_coupling.shape[1] != 2:
                raise ValueError("periodic_point_coupling ndarray must have shape (N, 2)")
            raw_pairs = [(int(r[0]), int(r[1])) for r in periodic_point_coupling]
        else:
            raise TypeError(
                "periodic_point_coupling must be dict-like, a (N, 2) numpy array, "
                "a CP_* array name string, or None"
            )

        midside_map = {}
        corner_map = {}
        for slave, master in raw_pairs:
            slave_is_mid = bool(midside_mask[slave])
            master_is_mid = bool(midside_mask[master])
            if slave_is_mid != master_is_mid:
                LOGGER.warning(
                    "periodic_point_coupling pair (%d -> %d) mixes a corner id "
                    "with a midside id; skipping this pair.", slave, master,
                )
                continue
            if slave_is_mid:
                midside_map[slave] = master
            else:
                corner_map[slave] = master
        return midside_map, corner_map

    @staticmethod
    def build_tree_cotree_for_midside_edges(
        ugrid: "vtk.vtkUnstructuredGrid",
        non_conducting_point_ids: Optional[Iterable[int]] = None,
        start_surface_point_ids: Optional[Iterable[int]] = None,
        periodic_point_coupling: Optional[object] = None,  # [CHG-2]
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
        periodic_slave_name: str = "PERIODIC_SLAVE",
    ) -> Tuple["vtk.vtkIntArray", "vtk.vtkIntArray", Dict[str, int]]:
        """
        Build a tree-cotree decomposition for second-order cells, identifying edges by
        their midside node ids. Results are vtkIntArray on PointData (1 at midside nodes).

        Args:
            ugrid: Mesh with second-order elements (quadratic). Edges are implicit.
            non_conducting_point_ids: Point ids for non-conducting regions where gauge tree constraints are built.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to nodes attached to currently selected elements (Selector.nsle).
            start_surface_point_ids: Point ids used as BFS roots (boundary-first traversal).
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                You can freely combine point ids from several distinct physical
                surfaces here (e.g. two different PEC boundaries bounding the
                same air region) — every seed that lands in a given connected
                component acts as a simultaneous root for that component, so
                using multiple surfaces only shortens the resulting tree, it
                never conflicts. If omitted, defaults to currently selected nodes.
            periodic_point_coupling: [CHG-2] Unified periodic coupling input
                containing BOTH corner and midside slave->master pairs mixed
                together. Accepts a dict-like mapping, an (N, 2) ndarray, or a
                CP_* PointData array name (grouped by matching CP value,
                corner and midside groups kept separate automatically).
            tree_name: Name for the tree output array (PointData).
            cotree_name: Name for the cotree output array (PointData).
            periodic_slave_name: Name for a PointData array marking
                periodic-slave midside nodes (1 = slave), so downstream
                constraint code can distinguish "tied to master" from
                "tree-fixed-to-zero" even though both may carry the tree bit.

        Returns:
            Tuple[vtkIntArray, vtkIntArray, Dict[str, int]]: (tree_arr, cotree_arr, info) where info contains
            diagnostic counts like num_edges, num_tree, num_cotree, skipped_unknown_cells,
            num_components, and a tree_identity_ok flag.
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
            Selector.nsle(ugrid, mode="s")
            point_sel_current = vtknp.vtk_to_numpy(point_sel_array)
            non_conducting_input: Set[int] = set(
                np.where(point_sel_current == 1)[0].astype(int).tolist()
            )
            point_sel_current[:] = point_sel_original
            point_sel_array.Modified()
        else:
            non_conducting_input = set(int(pid) for pid in non_conducting_point_ids)

        midside_arr = ugrid.GetPointData().GetArray("MidsideNodes")
        midside_mask = vtknp.vtk_to_numpy(midside_arr)

        raw_non_conducting = set(i for i in non_conducting_input if midside_mask[i] == 0)
        raw_start_surface = set(i for i in start_surface_input if midside_mask[i] == 0)

        # [CHG-2] Single normalization call splits corner/midside pairs.
        periodic_midside_raw, corner_slave_to_master = TreeGauge._normalize_periodic_point_coupling(
            periodic_point_coupling, ugrid, midside_mask
        )
        slave_to_master: Dict[int, int] = {}
        for slave, master in periodic_midside_raw.items():
            if midside_mask[slave] != 1 or midside_mask[master] != 1:
                continue
            slave_to_master[slave] = master

        tree_np = np.zeros(npts, dtype=np.int32)
        cotree_np = np.zeros(npts, dtype=np.int32)
        periodic_slave_np = np.zeros(npts, dtype=np.int32)

        edge_map: Dict[int, Tuple[int, int]] = {}
        raw_edge_map: Dict[int, Tuple[int, int]] = {}
        canonical_to_raw_mids: Dict[int, Set[int]] = {}
        unknown_cell_types: Dict[int, int] = {}

        for ci in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(ci)
            ctype = cell.GetCellType()
            if ctype not in EDGE_TABLES:
                unknown_cell_types[ctype] = unknown_cell_types.get(ctype, 0) + 1
                continue
            corner_pairs, mid_locals = EDGE_TABLES[ctype]
            ids = cell.GetPointIds()
            for (v0_local, v1_local), mid_local in zip(corner_pairs, mid_locals):
                v0 = ids.GetId(v0_local)
                v1 = ids.GetId(v1_local)
                raw_mid = ids.GetId(mid_local)
                if raw_mid not in raw_edge_map:
                    raw_edge_map[raw_mid] = (v0, v1)
                mid = TreeGauge._canonical_id(raw_mid, slave_to_master)

                canonical_to_raw_mids.setdefault(mid, set()).add(raw_mid)

                if mid not in edge_map:
                    edge_map[mid] = (v0, v1)

        if unknown_cell_types:
            total_skipped = sum(unknown_cell_types.values())
            LOGGER.warning(
                "Skipped %d cell(s) with unsupported (non-quadratic) cell types "
                "while building the tree-cotree graph: %s. Edges belonging to "
                "these cells are NOT included in the gauge graph — verify this "
                "is intended (e.g. mixed-order mesh, boundary tag elements).",
                total_skipped, unknown_cell_types,
            )

        non_conducting = {
            TreeGauge._canonical_id(pid, corner_slave_to_master)
            for pid in raw_non_conducting
        }
        start_surface = {
            TreeGauge._canonical_id(pid, corner_slave_to_master)
            for pid in raw_start_surface
        }
        for mid, (a, b) in list(edge_map.items()):
            edge_map[mid] = (
                TreeGauge._canonical_id(a, corner_slave_to_master),
                TreeGauge._canonical_id(b, corner_slave_to_master),
            )

        eligible_edges: Dict[int, Tuple[int, int]] = {}
        for mid, (a, b) in edge_map.items():
            if (a not in non_conducting) or (b not in non_conducting):
                for raw_mid in canonical_to_raw_mids[mid]:
                    cotree_np[raw_mid] = 1
                continue
            eligible_edges[mid] = (a, b)

        num_edges_before_periodic_collapse = sum(
            1
            for a, b in raw_edge_map.values()
            if (a in raw_non_conducting) and (b in raw_non_conducting)
        )

        adjacency: Dict[int, list[Tuple[int, int]]] = {}
        for mid, (a, b) in eligible_edges.items():
            adjacency.setdefault(a, []).append((b, mid))
            adjacency.setdefault(b, []).append((a, mid))

        edge_state: Dict[int, int] = {}  # 1=tree, 2=cotree
        visited: Set[int] = set()
        queue: deque[int] = deque()

        # Union-Find over graph nodes. Each seed starts as its OWN singleton
        # set (they are NOT pre-merged with each other) so that an edge
        # directly bridging two still-separate seed-rooted fragments is
        # correctly recognized as a genuine tree edge instead of always
        # being discarded as a "revisit". See dsu_bug note in repo memory.
        dsu_parent: Dict[int, int] = {}

        def dsu_find(x: int) -> int:
            dsu_parent.setdefault(x, x)
            while dsu_parent[x] != x:
                dsu_parent[x] = dsu_parent[dsu_parent[x]]
                x = dsu_parent[x]
            return x

        def dsu_union(a: int, b: int) -> bool:
            ra, rb = dsu_find(a), dsu_find(b)
            if ra == rb:
                return False
            dsu_parent[rb] = ra
            return True

        # True multi-source BFS: all seeds belonging to the same connected
        # component act as simultaneous roots. Pass point ids from as many
        # distinct physical surfaces as you like in start_surface_point_ids;
        # each becomes its own local root within its component, so the tree
        # stays shallow instead of only the first-encountered seed winning.
        def run_multi_source_bfs(seeds: Iterable[int]) -> None:
            local_seeds = [s for s in seeds if s not in visited]
            if not local_seeds:
                return
            for s in local_seeds:
                visited.add(s)
                dsu_find(s)
                if s in adjacency:
                    queue.append(s)
            while queue:
                node = queue.popleft()
                for nbr, mid in adjacency.get(node, []):
                    if mid in edge_state:
                        continue
                    if nbr not in visited:
                        edge_state[mid] = 1
                        visited.add(nbr)
                        dsu_union(node, nbr)
                        queue.append(nbr)
                    else:
                        # nbr already visited: this bridges two distinct
                        # multi-source fragments (tree edge, merges them) or
                        # closes a real cycle within one fragment (cotree).
                        # Without this DSU check every seed-adjacent edge
                        # would be miscounted as a cycle, undercounting the
                        # tree by (#seeds - #true_components) edges.
                        if dsu_union(node, nbr):
                            edge_state[mid] = 1
                        else:
                            edge_state[mid] = 2

        # Primary pass: all boundary/reference seeds fire together.
        run_multi_source_bfs(s for s in sorted(start_surface) if s in non_conducting)

        # Any remaining components (no start-surface seed reached them) get
        # an arbitrary single root each.
        for node in sorted(adjacency.keys()):
            if node not in visited:
                run_multi_source_bfs([node])

        for mid in eligible_edges:
            if mid not in edge_state:
                edge_state[mid] = 2  # safety net; should not trigger

        # Propagate classification to periodic-slave raw mids too, so they
        # carry the same tree/cotree state as their master. Downstream
        # constraint code must check PERIODIC_SLAVE before deciding whether
        # to zero a tree point or tie it to its master instead (see
        # tree_gauge() below for the reference implementation of that check).
        for mid, state in edge_state.items():
            for raw_mid in canonical_to_raw_mids[mid]:
                if state == 1:
                    tree_np[raw_mid] = 1
                else:
                    cotree_np[raw_mid] = 1
                if raw_mid in slave_to_master:
                    periodic_slave_np[raw_mid] = 1

        tree_np[cotree_np == 1] = 0

        tree_arr = vtknp.numpy_to_vtk(tree_np, deep=1, array_type=vtk.VTK_INT)
        tree_arr.SetName(tree_name)
        cotree_arr = vtknp.numpy_to_vtk(cotree_np, deep=1, array_type=vtk.VTK_INT)
        cotree_arr.SetName(cotree_name)
        periodic_slave_arr = vtknp.numpy_to_vtk(
            periodic_slave_np, deep=1, array_type=vtk.VTK_INT
        )
        periodic_slave_arr.SetName(periodic_slave_name)
        ugrid.GetPointData().RemoveArray(periodic_slave_name)
        ugrid.GetPointData().AddArray(periodic_slave_arr)

        num_edges_after_periodic_collapse = len(eligible_edges)

        # Spanning-forest identity check: for each connected component,
        # #tree_edges_in_component == #nodes_in_component - 1. Count
        # connected components directly (cheap, same complexity as the BFS
        # already performed).
        num_components = 0
        visited_for_count: Set[int] = set()
        for node in adjacency.keys():
            if node in visited_for_count:
                continue
            num_components += 1
            stack = [node]
            visited_for_count.add(node)
            while stack:
                cur = stack.pop()
                for nbr, _mid in adjacency.get(cur, []):
                    if nbr not in visited_for_count:
                        visited_for_count.add(nbr)
                        stack.append(nbr)

        num_tree_canonical = sum(1 for s in edge_state.values() if s == 1)
        num_nodes_in_graph = len(adjacency)
        expected_tree = num_nodes_in_graph - num_components if num_nodes_in_graph else 0
        tree_identity_ok = (num_tree_canonical == expected_tree)
        if not tree_identity_ok:
            LOGGER.error(
                "Tree-cotree identity check FAILED: got %d tree edges, expected "
                "%d (= %d nodes - %d components). This indicates a bug in the "
                "gauge graph construction (e.g. wrong non-conducting/PEC "
                "classification or a periodic-mapping error) — do not trust "
                "this gauge until resolved.",
                num_tree_canonical, expected_tree, num_nodes_in_graph, num_components,
            )

        info = {
            "num_points": npts,
            "num_cells": ugrid.GetNumberOfCells(),
            "num_edges": num_edges_before_periodic_collapse,
            "num_edges_after_periodic_collapse": num_edges_after_periodic_collapse,
            "num_tree": int(np.sum(tree_np)),
            "num_cotree": int(np.sum(cotree_np)),
            "num_periodic_slave_marked": int(np.sum(periodic_slave_np)),
            "skipped_unknown_cells": sum(unknown_cell_types.values()),
            "num_periodic_couplings": len(slave_to_master),
            "num_components": num_components,
            "num_tree_canonical": num_tree_canonical,
            "expected_tree_canonical": expected_tree,
            "tree_identity_ok": tree_identity_ok,
        }
        return tree_arr, cotree_arr, info

    @staticmethod
    def build_tree_cotree_on_grid(
        ugrid: "vtk.vtkUnstructuredGrid",
        non_conducting_point_ids: Optional[Iterable[int]] = None,
        start_surface_point_ids: Optional[Iterable[int]] = None,
        periodic_point_coupling: Optional[object] = None,  # [CHG-2]
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
        periodic_slave_name: str = "PERIODIC_SLAVE",
    ) -> Dict[str, int]:
        """
        Convenience wrapper that builds tree/cotree arrays and attaches them to the grid.
        """
        tree_arr, cotree_arr, info = TreeGauge.build_tree_cotree_for_midside_edges(
            ugrid,
            non_conducting_point_ids=non_conducting_point_ids,
            start_surface_point_ids=start_surface_point_ids,
            periodic_point_coupling=periodic_point_coupling,
            tree_name=tree_name,
            cotree_name=cotree_name,
            periodic_slave_name=periodic_slave_name,
        )
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
        periodic_point_coupling: Optional[object] = None,  # [CHG-2]
        array_name: str = "_Az",
        tree_name: str = "TREE",
        cotree_name: str = "COTREE",
        periodic_slave_name: str = "PERIODIC_SLAVE",
    ) -> Dict[str, int]:
        """
        Build tree/cotree masks and zero out the tree-marked entries in a
        PointData array — EXCEPT at periodic-slave points. [CHG-3]

        Periodic-slave points are tied to their master (A_slave = ± A_master)
        by the caller's own periodic constraint code, so this method leaves
        their `array_name` value untouched even if they carry the tree bit;
        zeroing them here would be redundant at best and wrong at worst if
        the periodic sign/offset is nonzero. Non-slave tree points are
        zeroed as before.

        Args:
            ugrid: Target mesh to analyze and annotate.
            array_name: PointData array to be gauge-fixed (tree entries set to zero).
            tree_name: Name for the output tree mask array.
            cotree_name: Name for the output cotree mask array.
            non_conducting_point_ids: Point ids where tree constraints are computed.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to nodes attached to currently selected elements (Selector.nsle).
            start_surface_point_ids: BFS root point ids; may span multiple surfaces.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                If omitted, defaults to currently selected nodes.
            periodic_point_coupling: Unified periodic coupling (corner and
                midside pairs mixed) as dict-like, ndarray (N,2), or CP_*
                PointData array name.
            periodic_slave_name: Name of the PointData array marking
                periodic-slave midside nodes.

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
            periodic_point_coupling=periodic_point_coupling,
            tree_name=tree_name,
            cotree_name=cotree_name,
            periodic_slave_name=periodic_slave_name,
        )
        if not info.get("tree_identity_ok", True):
            LOGGER.error(
                "tree_gauge(): tree-cotree identity check failed — see prior "
                "error log. Proceeding, but the resulting gauge is suspect."
            )
        LOGGER.info("Tree/cotree gauge built. Info: %s", info)

        tree_arr = ugrid.GetPointData().GetArray(tree_name)
        if tree_arr is not None:
            tree_np = vtknp.vtk_to_numpy(tree_arr)
            periodic_slave_arr = ugrid.GetPointData().GetArray(periodic_slave_name)
            if periodic_slave_arr is not None:
                periodic_slave_np = vtknp.vtk_to_numpy(periodic_slave_arr)
            else:
                periodic_slave_np = np.zeros_like(tree_np)
            # [CHG-3] Only zero tree points that are NOT periodic slaves.
            zero_mask = (tree_np == 1) & (periodic_slave_np == 0)
            arr_np[zero_mask] = 0
            arr.Modified()

        return info

    @staticmethod
    def add_midside_nodes_mask_array(
        ugrid: "vtk.vtkUnstructuredGrid", array_name: str = "MidsideNodes"
    ) -> "vtk.vtkIntArray":
        """
        Add a PointData array with 1 on midside nodes, 0 on corner nodes.
        """
        ensure_vtk_unstructured_grid(ugrid, arg_name="ugrid")
        npts = ugrid.GetNumberOfPoints()
        mask = np.zeros(npts, dtype=np.int32)

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
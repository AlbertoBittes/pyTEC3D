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
#   [CHG-4] Tree construction uses a multi-source shortest-path tree
#           (Dijkstra) on the non-conducting corner-node graph, with edge
#           cost = edge length only.
#
# Carried over from the previous revision:
#
#   [FIX-3] Root selection is true multi-source: all start-surface seeds in
#           a connected component act as simultaneous roots for that
#           component; components with no seed fall back to an arbitrary root.
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
#       start_surface_point_ids={...},          # multi-source seeds
#       periodic_point_coupling=periodic_map,   # corner + midside ids mixed
#   )
#   ugrid.GetPointData().AddArray(tree_arr)
#   ugrid.GetPointData().AddArray(cotree_arr)
#
# (c) Provide as-is; tailor as needed for your mesh conventions.

from typing import Dict, Iterable, Mapping, Optional, Sequence, Set, Tuple, TYPE_CHECKING, cast

import heapq
import numpy as np
import vtk
import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

from .logging_utils import get_logger
from .selector import Selector
from .validation_utils import ensure_vtk_unstructured_grid


LOGGER = get_logger("pytec3d.tree_gauge")

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
    """Build and apply tree/cotree gauges on VTK meshes."""

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
                ids = np.flatnonzero(midside_mask == mask_value)
                if ids.size == 0:
                    return {}

                vals = cp_np[ids]
                if np.issubdtype(vals.dtype, np.floating):
                    valid = np.isfinite(vals) & (vals != 0)
                else:
                    valid = vals != 0

                ids = ids[valid]
                vals = vals[valid]
                if ids.size < 2:
                    return {}

                order = np.argsort(vals, kind="mergesort")
                vals_sorted = vals[order]
                ids_sorted = ids[order]

                is_new = np.empty(vals_sorted.size, dtype=bool)
                is_new[0] = True
                is_new[1:] = vals_sorted[1:] != vals_sorted[:-1]
                starts = np.flatnonzero(is_new)
                ends = np.empty_like(starts)
                ends[:-1] = starts[1:]
                ends[-1] = vals_sorted.size

                mapping: Dict[int, int] = {}
                for s, e in zip(starts, ends):
                    if (e - s) < 2:
                        continue
                    members = ids_sorted[s:e]
                    master = int(members.min())
                    slaves = members[members != master]
                    if slaves.size:
                        mapping.update(zip(slaves.tolist(), [master] * int(slaves.size)))
                return mapping

            midside_map = build_groups(1)
            corner_map = build_groups(0)
            return midside_map, corner_map

        # dict / ndarray path: pairs may mix corner and midside ids together.
        raw_pairs: np.ndarray
        if isinstance(periodic_point_coupling, Mapping):
            raw_pairs = np.array(
                [(int(s), int(m)) for s, m in periodic_point_coupling.items()],
                dtype=np.int64,
            )
        elif isinstance(periodic_point_coupling, np.ndarray):
            coupling_arr = cast(np.ndarray, periodic_point_coupling)
            if coupling_arr.ndim != 2 or coupling_arr.shape[1] != 2:
                raise ValueError("periodic_point_coupling ndarray must have shape (N, 2)")
            raw_pairs = coupling_arr.astype(np.int64, copy=False)
        else:
            raise TypeError(
                "periodic_point_coupling must be dict-like, a (N, 2) numpy array, "
                "a CP_* array name string, or None"
            )

        if raw_pairs.size == 0:
            return {}, {}

        slave = raw_pairs[:, 0]
        master = raw_pairs[:, 1]
        slave_is_mid = midside_mask[slave] == 1
        master_is_mid = midside_mask[master] == 1
        same_type = slave_is_mid == master_is_mid

        skipped = int(np.count_nonzero(~same_type))
        if skipped:
            LOGGER.warning(
                "periodic_point_coupling skipped %d mixed corner/midside pair(s).",
                skipped,
            )

        midside_pairs = same_type & slave_is_mid
        corner_pairs = same_type & (~slave_is_mid)

        midside_map = dict(
            zip(slave[midside_pairs].tolist(), master[midside_pairs].tolist())
        )
        corner_map = dict(
            zip(slave[corner_pairs].tolist(), master[corner_pairs].tolist())
        )
        return midside_map, corner_map

    @staticmethod
    def _canonicalize_ids(
        ids: np.ndarray,
        slave_to_master: Dict[int, int],
        cache: Dict[int, int],
    ) -> np.ndarray:
        """Map point ids to canonical masters using periodic slave-to-master links.

        Args:
            ids: Point ids to canonicalize.
            slave_to_master: Mapping from slave ids to their direct masters.
            cache: Mutable memoization map of id -> canonical id.

        Returns:
            np.ndarray: Canonicalized ids with dtype int64.

        Raises:
            ValueError: If a cycle is detected in slave_to_master links.
        """
        if ids.size == 0 or not slave_to_master:
            return ids.astype(np.int64, copy=False)

        out = np.empty(ids.size, dtype=np.int64)
        for i in range(ids.size):
            pid = int(ids[i])
            cached = cache.get(pid)
            if cached is not None:
                out[i] = cached
                continue

            path = []
            seen = set()
            cur = pid
            while True:
                cached_cur = cache.get(cur)
                if cached_cur is not None:
                    root = cached_cur
                    break

                nxt = slave_to_master.get(cur)
                if nxt is None or int(nxt) == cur:
                    root = cur
                    break

                if cur in seen:
                    raise ValueError(
                        f"Cycle detected in periodic point coupling at point id {cur}."
                    )
                seen.add(cur)
                path.append(cur)
                cur = int(nxt)

            cache[pid] = root
            for p in path:
                cache[p] = root
            out[i] = root

        return out

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
            start_surface_point_ids: Multi-source seeds for the
                shortest-path (Dijkstra) builder.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                Components with no seed fall back to an arbitrary local root.
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
            start_surface_ids = np.flatnonzero(point_sel_original == 1).astype(np.int64)
        else:
            start_surface_ids = np.fromiter(
                (int(pid) for pid in start_surface_point_ids), dtype=np.int64
            )

        if non_conducting_point_ids is None:
            Selector.nsle(ugrid, mode="s")
            point_sel_current = vtknp.vtk_to_numpy(point_sel_array)
            non_conducting_ids = np.flatnonzero(point_sel_current == 1).astype(np.int64)
            point_sel_current[:] = point_sel_original
            point_sel_array.Modified()
        else:
            non_conducting_ids = np.fromiter(
                (int(pid) for pid in non_conducting_point_ids), dtype=np.int64
            )

        midside_arr = ugrid.GetPointData().GetArray("MidsideNodes")
        if midside_arr is None:
            raise ValueError("PointData array 'MidsideNodes' was not found on the mesh.")
        midside_mask = vtknp.vtk_to_numpy(midside_arr)

        raw_non_conducting = non_conducting_ids[midside_mask[non_conducting_ids] == 0]
        raw_start_surface = start_surface_ids[midside_mask[start_surface_ids] == 0]

        # [CHG-2] Single normalization call splits corner/midside pairs.
        periodic_midside_raw, corner_slave_to_master = TreeGauge._normalize_periodic_point_coupling(
            periodic_point_coupling, ugrid, midside_mask
        )
        slave_to_master = {
            int(slave): int(master)
            for slave, master in periodic_midside_raw.items()
            if midside_mask[slave] == 1 and midside_mask[master] == 1
        }

        tree_np = np.zeros(npts, dtype=np.int32)
        cotree_np = np.zeros(npts, dtype=np.int32)
        periodic_slave_np = np.zeros(npts, dtype=np.int32)

        points_data = ugrid.GetPoints()
        if points_data is None:
            raise ValueError("Mesh has no points data.")
        points_np = vtknp.vtk_to_numpy(points_data.GetData())

        edge_map: Dict[int, Tuple[int, int]] = {}
        raw_edge_map: Dict[int, Tuple[int, int]] = {}
        raw_to_canonical_mid: Dict[int, int] = {}
        unknown_cell_types: Dict[int, int] = {}
        corner_cache: Dict[int, int] = {}
        midside_cache: Dict[int, int] = {}

        def canonical_corner(pid: int) -> int:
            cached = corner_cache.get(pid)
            if cached is not None:
                return cached
            out = TreeGauge._canonical_id(pid, corner_slave_to_master)
            corner_cache[pid] = out
            return out

        def canonical_mid(pid: int) -> int:
            cached = midside_cache.get(pid)
            if cached is not None:
                return cached
            out = TreeGauge._canonical_id(pid, slave_to_master)
            midside_cache[pid] = out
            return out

        for ci in range(ugrid.GetNumberOfCells()):
            cell = ugrid.GetCell(ci)
            ctype = cell.GetCellType()
            if ctype not in EDGE_TABLES:
                unknown_cell_types[ctype] = unknown_cell_types.get(ctype, 0) + 1
                continue
            corner_pairs, mid_locals = EDGE_TABLES[ctype]
            ids = cell.GetPointIds()

            for (v0_local, v1_local), mid_local in zip(corner_pairs, mid_locals):
                v0 = int(ids.GetId(v0_local))
                v1 = int(ids.GetId(v1_local))
                raw_mid = int(ids.GetId(mid_local))
                if raw_mid not in raw_edge_map:
                    raw_edge_map[raw_mid] = (v0, v1)
                mid = canonical_mid(raw_mid)
                raw_to_canonical_mid[raw_mid] = mid

                if mid not in edge_map:
                    edge_map[mid] = (canonical_corner(v0), canonical_corner(v1))

        if unknown_cell_types:
            total_skipped = sum(unknown_cell_types.values())
            LOGGER.warning(
                "Skipped %d cell(s) with unsupported (non-quadratic) cell types "
                "while building the tree-cotree graph: %s. Edges belonging to "
                "these cells are NOT included in the gauge graph — verify this "
                "is intended (e.g. mixed-order mesh, boundary tag elements).",
                total_skipped, unknown_cell_types,
            )

        non_conducting = np.unique(
            TreeGauge._canonicalize_ids(raw_non_conducting, corner_slave_to_master, corner_cache)
        )
        start_surface = np.unique(
            TreeGauge._canonicalize_ids(raw_start_surface, corner_slave_to_master, corner_cache)
        )

        if non_conducting.size:
            non_cond_lookup = np.zeros(npts, dtype=bool)
            non_cond_lookup[non_conducting] = True
        else:
            non_cond_lookup = np.zeros(npts, dtype=bool)

        edge_mid_ids = np.fromiter(edge_map.keys(), dtype=np.int64)
        edge_uv = np.array(list(edge_map.values()), dtype=np.int64)

        if edge_mid_ids.size:
            eligible_mask = non_cond_lookup[edge_uv[:, 0]] & non_cond_lookup[edge_uv[:, 1]]
            eligible_mid_ids = edge_mid_ids[eligible_mask]
            eligible_uv = edge_uv[eligible_mask]
        else:
            eligible_mid_ids = np.empty(0, dtype=np.int64)
            eligible_uv = np.empty((0, 2), dtype=np.int64)

        if raw_edge_map:
            raw_uv = np.array(list(raw_edge_map.values()), dtype=np.int64)
            raw_non_cond_lookup = np.zeros(npts, dtype=bool)
            raw_non_cond_lookup[raw_non_conducting] = True
            num_edges_before_periodic_collapse = int(
                np.count_nonzero(
                    raw_non_cond_lookup[raw_uv[:, 0]] & raw_non_cond_lookup[raw_uv[:, 1]]
                )
            )
        else:
            num_edges_before_periodic_collapse = 0

        edge_state_lookup = np.zeros(npts, dtype=np.int8)  # 1=tree, 2=cotree
        edge_state_lookup[eligible_mid_ids] = 2

        num_components = 0
        num_roots_total = 0

        if eligible_mid_ids.size:
            u = eligible_uv[:, 0]
            v = eligible_uv[:, 1]
            lengths = np.linalg.norm(points_np[u] - points_np[v], axis=1)
            bad_len = (~np.isfinite(lengths)) | (lengths <= 0.0)
            if np.any(bad_len):
                lengths[bad_len] = 1.0e-12

            node_ids = np.unique(np.concatenate((u, v)))
            u_idx = np.searchsorted(node_ids, u)
            v_idx = np.searchsorted(node_ids, v)

            src = np.concatenate((u_idx, v_idx))
            dst = np.concatenate((v_idx, u_idx))
            mid_dir = np.concatenate((eligible_mid_ids, eligible_mid_ids))
            w = np.concatenate((lengths, lengths))

            order = np.argsort(src, kind="mergesort")
            src = src[order]
            dst = dst[order]
            mid_dir = mid_dir[order]
            w = w[order]

            n_nodes = int(node_ids.size)
            indptr = np.zeros(n_nodes + 1, dtype=np.int64)
            np.add.at(indptr, src + 1, 1)
            np.cumsum(indptr, out=indptr)

            start_lookup = np.zeros(npts, dtype=bool)
            if start_surface.size:
                start_lookup[start_surface] = True

            visited = np.zeros(n_nodes, dtype=bool)
            tol = 1.0e-14

            for root in range(n_nodes):
                if visited[root]:
                    continue

                num_components += 1
                stack = [root]
                visited[root] = True
                component = []

                while stack:
                    cur = stack.pop()
                    component.append(cur)
                    row_start = int(indptr[cur])
                    row_end = int(indptr[cur + 1])
                    nbrs = dst[row_start:row_end]
                    for nbr in nbrs:
                        if not visited[nbr]:
                            visited[nbr] = True
                            stack.append(int(nbr))

                component_np = np.array(component, dtype=np.int64)
                component_global = node_ids[component_np]
                seeds = component_np[start_lookup[component_global]]
                if seeds.size == 0:
                    seeds = np.array([root], dtype=np.int64)

                num_roots_total += int(seeds.size)

                dist = np.full(n_nodes, np.inf, dtype=np.float64)
                parent_node = np.full(n_nodes, -1, dtype=np.int64)
                parent_mid = np.full(n_nodes, -1, dtype=np.int64)
                heap: list[Tuple[float, int]] = []
                for s in seeds:
                    s_int = int(s)
                    dist[s_int] = 0.0
                    heapq.heappush(heap, (0.0, s_int))

                while heap:
                    dcur, node = heapq.heappop(heap)
                    if dcur > dist[node] + tol:
                        continue

                    row_start = int(indptr[node])
                    row_end = int(indptr[node + 1])
                    for j in range(row_start, row_end):
                        nbr = int(dst[j])
                        nd = dcur + float(w[j])
                        if nd < dist[nbr] - tol:
                            dist[nbr] = nd
                            parent_node[nbr] = node
                            parent_mid[nbr] = int(mid_dir[j])
                            heapq.heappush(heap, (nd, nbr))
                        elif abs(nd - dist[nbr]) <= tol and parent_mid[nbr] >= 0:
                            old_mid = int(parent_mid[nbr])
                            old_parent = int(parent_node[nbr])
                            new_mid = int(mid_dir[j])
                            if (new_mid < old_mid) or (
                                (new_mid == old_mid) and (node < old_parent)
                            ):
                                parent_node[nbr] = node
                                parent_mid[nbr] = new_mid

                tree_mids = np.unique(parent_mid[component_np][parent_mid[component_np] >= 0])
                edge_state_lookup[tree_mids] = 1

        raw_mids = np.fromiter(raw_to_canonical_mid.keys(), dtype=np.int64)
        raw_canonical = np.fromiter(raw_to_canonical_mid.values(), dtype=np.int64)
        if raw_mids.size:
            raw_is_tree = edge_state_lookup[raw_canonical] == 1
            tree_np[raw_mids[raw_is_tree]] = 1
            cotree_np[raw_mids[~raw_is_tree]] = 1

        if slave_to_master:
            periodic_slave_np[np.fromiter(slave_to_master.keys(), dtype=np.int64)] = 1

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

        num_edges_after_periodic_collapse = int(eligible_mid_ids.size)

        num_tree_canonical = int(np.count_nonzero(edge_state_lookup == 1))
        num_nodes_in_graph = 0
        if eligible_mid_ids.size:
            num_nodes_in_graph = int(np.unique(eligible_uv).size)

        expected_tree = num_nodes_in_graph - num_roots_total if num_nodes_in_graph else 0
        tree_identity_ok = (num_tree_canonical == expected_tree)
        if not tree_identity_ok:
            LOGGER.error(
                "Tree-cotree identity check FAILED: got %d tree edges, expected "
                "%d (= %d nodes - %d roots across %d components). This indicates "
                "a bug in the gauge graph construction (e.g. wrong "
                "non-conducting/PEC classification or a periodic-mapping error) "
                "— do not trust this gauge until resolved.",
                num_tree_canonical,
                expected_tree,
                num_nodes_in_graph,
                num_roots_total,
                num_components,
            )

        info = {
            "tree_build_algorithm": "multi_source_dijkstra_edge_length",
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
            "num_roots": num_roots_total,
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
            start_surface_point_ids: Multi-source seeds for the
                shortest-path (Dijkstra) builder.
                Mixed corner+midside inputs are accepted; midside ids are ignored.
                Components with no seed fall back to an arbitrary local root.
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
            "Tree gauge input: non-conducting points=%d, start-surface points=%d "
            "(used as multi-source seeds for edge-length Dijkstra)",
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

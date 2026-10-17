"""Tree-cotree gauge of edge (H(curl)) fields on second-order meshes.

An edge of a quadratic cell is identified by its midside point, so the gauge is written on point arrays: 1 at the
midside points of the tree (or cotree) edges. The tree is a spanning forest of the edges of the non-conducting region:
a multi-source shortest-path tree (Dijkstra, edge length as cost) grown from the start-surface points, one arbitrary
root for each connected component without such points. Fixing the edge values of the tree to zero gauges the
vector potential.

Periodic couplings join slave points to their masters: the slave edges and corners are merged with their masters
before the tree is built, and slave midside points are marked (``PERIODIC_SLAVE``) so that the gauge leaves them to
the periodic constraint instead of fixing them to zero.

Known limitation: no cohomology cuts. They are not needed for an A or A-V formulation driven by a divergence-free
source current, even in a non-simply-connected non-conducting region; they would be for a scalar potential.

Supported cells: quadratic edge, triangle, quad, tetrahedron, wedge, pyramid and hexahedron.
"""

from __future__ import annotations

import heapq
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonCore import VTK_INT, vtkIntArray
from vtkmodules.vtkCommonDataModel import (
    VTK_QUADRATIC_EDGE,
    VTK_QUADRATIC_HEXAHEDRON,
    VTK_QUADRATIC_PYRAMID,
    VTK_QUADRATIC_QUAD,
    VTK_QUADRATIC_TETRA,
    VTK_QUADRATIC_TRIANGLE,
    VTK_QUADRATIC_WEDGE,
    vtkUnstructuredGrid,
)

from . import mesh as meshes
from . import selection
from ._logging import get_logger

LOGGER = get_logger("pytec3d.tree_gauge")

#: Point array marking the midside points (1) of the quadratic cells.
MIDSIDE_ARRAY = "MidsideNodes"

#: For each quadratic cell type: the local corner pairs of its edges and the local midside point of each edge.
EDGE_TABLES: dict[int, tuple[tuple[tuple[int, int], ...], tuple[int, ...]]] = {
    VTK_QUADRATIC_EDGE: (((0, 1),), (2,)),
    VTK_QUADRATIC_TRIANGLE: (((0, 1), (1, 2), (2, 0)), (3, 4, 5)),
    VTK_QUADRATIC_QUAD: (((0, 1), (1, 2), (2, 3), (3, 0)), (4, 5, 6, 7)),
    VTK_QUADRATIC_TETRA: (((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)), (4, 5, 6, 7, 8, 9)),
    VTK_QUADRATIC_WEDGE: (
        ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)),
        (6, 7, 8, 9, 10, 11, 12, 13, 14),
    ),
    VTK_QUADRATIC_PYRAMID: (
        ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)),
        (5, 6, 7, 8, 9, 10, 11, 12),
    ),
    VTK_QUADRATIC_HEXAHEDRON: (
        ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)),
        (8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
    ),
}

_MAX_EDGES_PER_CELL = 16


# ---------------------------------------------------------------------------------------------------- edges


def _cell_edges(mesh: vtkUnstructuredGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[int, int]]:
    """``(v0, v1, midside)`` of every edge of every quadratic cell, in cell order, and the skipped cell types.

    Returns:
        Three arrays (corner, corner, midside point id), one entry per (cell, local edge), and the count of cells of
        each unsupported type.
    """
    offsets, connectivity = meshes.cell_connectivity(mesh)
    types = meshes.cell_types(mesh)
    keys, v0, v1, mid = [], [], [], []
    for cell_type, (pairs, mids) in EDGE_TABLES.items():
        cells = np.flatnonzero(types == cell_type)
        if cells.size == 0:
            continue
        starts = offsets[cells][:, np.newaxis]
        local = np.asarray(pairs)
        keys.append((cells[:, np.newaxis] * _MAX_EDGES_PER_CELL + np.arange(len(mids))).ravel())
        v0.append(connectivity[starts + local[:, 0]].ravel())
        v1.append(connectivity[starts + local[:, 1]].ravel())
        mid.append(connectivity[starts + np.asarray(mids)].ravel())
    unsupported_types, counts = np.unique(types[~np.isin(types, list(EDGE_TABLES))], return_counts=True)
    skipped = {int(t): int(c) for t, c in zip(unsupported_types, counts, strict=True)}
    if not keys:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty, empty, skipped
    order = np.argsort(np.concatenate(keys), kind="stable")
    return np.concatenate(v0)[order], np.concatenate(v1)[order], np.concatenate(mid)[order], skipped


def add_midside_nodes_mask(mesh: vtkUnstructuredGrid, name: str = MIDSIDE_ARRAY) -> vtkIntArray:
    """Add (or replace) a point array with 1 at the midside points of the quadratic cells, 0 elsewhere."""
    meshes.check_grid(mesh)
    mask = np.zeros(mesh.GetNumberOfPoints(), dtype=np.int32)
    _, _, midside, _ = _cell_edges(mesh)
    mask[midside] = 1
    array = vtknp.numpy_to_vtk(mask, deep=True, array_type=VTK_INT)
    array.SetName(name)
    mesh.GetPointData().AddArray(array)
    return array


# ---------------------------------------------------------------------------------------------------- periodic maps


def _follow(point_id: int, slave_to_master: Mapping[int, int]) -> int:
    """The master at the end of the slave -> master chain of a point."""
    seen: set[int] = set()
    current = int(point_id)
    while current in slave_to_master:
        if current in seen:
            raise ValueError(f"Cycle in the periodic point coupling at point {current}.")
        seen.add(current)
        following = int(slave_to_master[current])
        if following == current:
            break
        current = following
    return current


def _canonical_map(n_points: int, slave_to_master: Mapping[int, int]) -> np.ndarray:
    """Point id -> its final master (itself when not a slave)."""
    canonical = np.arange(n_points, dtype=np.int64)
    for slave in slave_to_master:
        canonical[int(slave)] = _follow(slave, slave_to_master)
    return canonical


def _groups_to_masters(ids: np.ndarray, values: np.ndarray) -> dict[int, int]:
    """Points with the same nonzero CP value form a group whose master is the lowest id; slave -> master."""
    finite = np.isfinite(values) if np.issubdtype(values.dtype, np.floating) else np.ones(values.shape, dtype=bool)
    keep = finite & (values != 0)
    ids, values = ids[keep], values[keep]
    if ids.size < 2:
        return {}
    order = np.argsort(values, kind="mergesort")
    ids, values = ids[order], values[order]
    starts = np.flatnonzero(np.concatenate(([True], values[1:] != values[:-1])))
    ends = np.append(starts[1:], values.size)
    mapping: dict[int, int] = {}
    for start, end in zip(starts, ends, strict=True):
        if end - start < 2:
            continue
        members = ids[start:end]
        master = int(members.min())
        mapping.update(dict.fromkeys(members[members != master].tolist(), master))
    return mapping


def _periodic_maps(
    coupling: Any, mesh: vtkUnstructuredGrid, midside: np.ndarray
) -> tuple[dict[int, int], dict[int, int]]:
    """Split a periodic point coupling into (midside slave -> master, corner slave -> master).

    ``coupling`` is None, the name of a CP_* point array (points sharing a nonzero value are coupled, the lowest id
    being the master), a mapping slave -> master, or an (N, 2) array of (slave, master) rows. Pairs mixing a corner
    and a midside point are dropped with a warning.
    """
    if coupling is None:
        return {}, {}
    if isinstance(coupling, str):
        if not coupling.startswith("CP_"):
            raise ValueError("periodic_point_coupling must name a CP_* point array.")
        values = meshes.array_values(mesh, coupling, "node")
        return (
            _groups_to_masters(np.flatnonzero(midside == 1), values[midside == 1]),
            _groups_to_masters(np.flatnonzero(midside == 0), values[midside == 0]),
        )
    if isinstance(coupling, Mapping):
        pairs = np.array([(int(s), int(m)) for s, m in coupling.items()], dtype=np.int64).reshape(-1, 2)
    elif isinstance(coupling, np.ndarray) and coupling.ndim == 2 and coupling.shape[1] == 2:
        pairs = coupling.astype(np.int64, copy=False)
    else:
        raise TypeError(
            "periodic_point_coupling must be None, a CP_* array name, a slave -> master mapping or an (N, 2) array."
        )
    if pairs.size == 0:
        return {}, {}
    slave, master = pairs[:, 0], pairs[:, 1]
    slave_mid, master_mid = midside[slave] == 1, midside[master] == 1
    same = slave_mid == master_mid
    if not np.all(same):
        LOGGER.warning(
            "periodic_point_coupling: %d pair(s) mixing a corner and a midside point skipped.", int((~same).sum())
        )
    mid_pairs, corner_pairs = same & slave_mid, same & ~slave_mid
    return (
        dict(zip(slave[mid_pairs].tolist(), master[mid_pairs].tolist(), strict=True)),
        dict(zip(slave[corner_pairs].tolist(), master[corner_pairs].tolist(), strict=True)),
    )


# ---------------------------------------------------------------------------------------------------- tree


def _shortest_path_forest(
    edge_mids: np.ndarray, edge_uv: np.ndarray, coordinates: np.ndarray, n_points: int, start_surface: np.ndarray
) -> tuple[np.ndarray, int, int]:
    """Tree edges (their midside ids) of the multi-source shortest-path forest of the edge graph.

    Returns:
        ``(tree_mids, components, roots)``.
    """
    u, v = edge_uv[:, 0], edge_uv[:, 1]
    lengths = np.linalg.norm(coordinates[u] - coordinates[v], axis=1)
    lengths[~np.isfinite(lengths) | (lengths <= 0.0)] = 1.0e-12

    nodes = np.unique(np.concatenate((u, v)))
    u_index, v_index = np.searchsorted(nodes, u), np.searchsorted(nodes, v)
    source = np.concatenate((u_index, v_index))
    order = np.argsort(source, kind="mergesort")
    source = source[order]
    target = np.concatenate((v_index, u_index))[order]
    edge_mid = np.concatenate((edge_mids, edge_mids))[order]
    weight = np.concatenate((lengths, lengths))[order]

    n_nodes = int(nodes.size)
    row_start = np.zeros(n_nodes + 1, dtype=np.int64)
    np.add.at(row_start, source + 1, 1)
    np.cumsum(row_start, out=row_start)
    target_list, mid_list, weight_list = target.tolist(), edge_mid.tolist(), weight.tolist()

    is_start = np.zeros(n_points, dtype=bool)
    is_start[start_surface] = True
    visited = np.zeros(n_nodes, dtype=bool)
    tolerance = 1.0e-14
    tree: list[np.ndarray] = []
    components = roots = 0

    for root in range(n_nodes):
        if visited[root]:
            continue
        components += 1
        stack, component = [root], []
        visited[root] = True
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in target_list[row_start[current]:row_start[current + 1]]:
                if not visited[neighbour]:
                    visited[neighbour] = True
                    stack.append(neighbour)

        members = np.array(component, dtype=np.int64)
        seeds = members[is_start[nodes[members]]]
        if seeds.size == 0:
            seeds = np.array([root], dtype=np.int64)
        roots += int(seeds.size)

        distance = np.full(n_nodes, np.inf)
        parent = np.full(n_nodes, -1, dtype=np.int64)
        parent_mid = np.full(n_nodes, -1, dtype=np.int64)
        heap: list[tuple[float, int]] = []
        for seed in seeds.tolist():
            distance[seed] = 0.0
            heapq.heappush(heap, (0.0, seed))
        while heap:
            current_distance, node = heapq.heappop(heap)
            if current_distance > distance[node] + tolerance:
                continue
            for j in range(row_start[node], row_start[node + 1]):
                neighbour = target_list[j]
                candidate = current_distance + weight_list[j]
                if candidate < distance[neighbour] - tolerance:
                    distance[neighbour] = candidate
                    parent[neighbour] = node
                    parent_mid[neighbour] = mid_list[j]
                    heapq.heappush(heap, (candidate, neighbour))
                elif abs(candidate - distance[neighbour]) <= tolerance and parent_mid[neighbour] >= 0:
                    # Equal paths: keep the lowest edge id, then the lowest parent, so the tree is reproducible.
                    old_mid, new_mid = int(parent_mid[neighbour]), mid_list[j]
                    if new_mid < old_mid or (new_mid == old_mid and node < int(parent[neighbour])):
                        parent[neighbour] = node
                        parent_mid[neighbour] = new_mid
        member_mids = parent_mid[members]
        tree.append(np.unique(member_mids[member_mids >= 0]))
    tree_mids = np.concatenate(tree) if tree else np.empty(0, dtype=np.int64)
    return tree_mids, components, roots


def build_tree_cotree(
    mesh: vtkUnstructuredGrid,
    non_conducting_point_ids: Iterable[int] | None = None,
    start_surface_point_ids: Iterable[int] | None = None,
    periodic_point_coupling: Any = None,
    tree_name: str = "TREE",
    cotree_name: str = "COTREE",
    periodic_slave_name: str = "PERIODIC_SLAVE",
) -> tuple[vtkIntArray, vtkIntArray, dict[str, Any]]:
    """Build the tree-cotree decomposition of the edges of the non-conducting region.

    Args:
        mesh: Mesh of quadratic cells, with the ``MidsideNodes`` point array (:func:`add_midside_nodes_mask`).
        non_conducting_point_ids: Points of the non-conducting region (midside points are ignored). Default: the
            points of the selected cells.
        start_surface_point_ids: Seeds of the shortest-path tree (midside points are ignored). Default: the
            selected points.
        periodic_point_coupling: None, a CP_* point array name, a slave -> master mapping or an (N, 2) array of
            (slave, master) rows; corner and midside pairs together.
        tree_name: Name of the returned tree array.
        cotree_name: Name of the returned cotree array.
        periodic_slave_name: Point array, added to the mesh, marking the periodic slave midside points.

    Returns:
        ``(tree, cotree, info)``: point arrays (1 at the midside points of the tree / cotree edges, not added to the
        mesh) and diagnostics: edge, tree, cotree, component and root counts, and ``tree_identity_ok``, False when
        #tree edges != #nodes - #roots (a broken graph).

    Raises:
        MissingArrayError: If ``MidsideNodes`` or a named CP array is missing.
        ValueError: If the periodic coupling has a cycle.
    """
    meshes.check_grid(mesh)
    n_points = mesh.GetNumberOfPoints()
    selected_points = selection.selected_ids(mesh, "node")
    start_ids = selected_points if start_surface_point_ids is None else np.fromiter(
        (int(i) for i in start_surface_point_ids), dtype=np.int64
    )
    if non_conducting_point_ids is None:
        current = selection.selection_mask(mesh, "node")
        selection.nsle(mesh, mode="s")
        region_ids = selection.selected_ids(mesh, "node")
        selection.apply_mask(mesh, current, "node", "s")
    else:
        region_ids = np.fromiter((int(i) for i in non_conducting_point_ids), dtype=np.int64)

    midside = meshes.array_values(mesh, MIDSIDE_ARRAY, "node")
    region_corners = region_ids[midside[region_ids] == 0]
    start_corners = start_ids[midside[start_ids] == 0]

    midside_masters, corner_masters = _periodic_maps(periodic_point_coupling, mesh, midside)
    midside_masters = {s: m for s, m in midside_masters.items() if midside[s] == 1 and midside[m] == 1}
    canonical_mid = _canonical_map(n_points, midside_masters)
    canonical_corner = _canonical_map(n_points, corner_masters)

    v0, v1, raw_mids, skipped = _cell_edges(mesh)
    if skipped:
        LOGGER.warning(
            "Tree-cotree: %d cell(s) of unsupported (non-quadratic) types skipped: %s. Their edges are not in the "
            "gauge graph; check that this is intended (mixed-order mesh, boundary elements).",
            sum(skipped.values()), skipped,
        )

    # One edge per canonical midside point, taken where it first appears (cell order), as corner pairs.
    mids = canonical_mid[raw_mids]
    _, first = np.unique(mids, return_index=True)
    first.sort()
    edge_mids = mids[first]
    edge_uv = np.column_stack((canonical_corner[v0[first]], canonical_corner[v1[first]]))

    in_region = np.zeros(n_points, dtype=bool)
    in_region[np.unique(canonical_corner[region_corners])] = True
    eligible = in_region[edge_uv[:, 0]] & in_region[edge_uv[:, 1]] if edge_mids.size else np.zeros(0, dtype=bool)
    eligible_mids, eligible_uv = edge_mids[eligible], edge_uv[eligible]

    raw_in_region = np.zeros(n_points, dtype=bool)
    raw_in_region[region_corners] = True
    _, raw_first = np.unique(raw_mids, return_index=True)
    edges_before_collapse = int(np.count_nonzero(raw_in_region[v0[raw_first]] & raw_in_region[v1[raw_first]]))

    state = np.zeros(n_points, dtype=np.int8)  # 1 tree, 2 cotree, by canonical midside point
    state[eligible_mids] = 2
    components = roots = 0
    if eligible_mids.size:
        start_surface = np.unique(canonical_corner[start_corners])
        coordinates = meshes.points(mesh)
        tree_mids, components, roots = _shortest_path_forest(
            eligible_mids, eligible_uv, coordinates, n_points, start_surface
        )
        state[tree_mids] = 1

    tree = np.zeros(n_points, dtype=np.int32)
    cotree = np.zeros(n_points, dtype=np.int32)
    unique_raw = np.unique(raw_mids)
    raw_is_tree = state[canonical_mid[unique_raw]] == 1
    tree[unique_raw[raw_is_tree]] = 1
    cotree[unique_raw[~raw_is_tree]] = 1
    tree[cotree == 1] = 0

    periodic_slave = np.zeros(n_points, dtype=np.int32)
    if midside_masters:
        periodic_slave[np.fromiter(midside_masters, dtype=np.int64)] = 1
    meshes.remove_array(mesh, periodic_slave_name, "node")
    slave_array = vtknp.numpy_to_vtk(periodic_slave, deep=True, array_type=VTK_INT)
    slave_array.SetName(periodic_slave_name)
    mesh.GetPointData().AddArray(slave_array)

    tree_canonical = int(np.count_nonzero(state == 1))
    graph_nodes = int(np.unique(eligible_uv).size) if eligible_mids.size else 0
    expected_tree = graph_nodes - roots if graph_nodes else 0
    identity_ok = tree_canonical == expected_tree
    if not identity_ok:
        LOGGER.error(
            "Tree-cotree identity check failed: %d tree edges, expected %d (%d nodes - %d roots in %d components). "
            "The gauge graph is broken (wrong non-conducting region or periodic coupling): do not trust this gauge.",
            tree_canonical, expected_tree, graph_nodes, roots, components,
        )

    tree_array = vtknp.numpy_to_vtk(tree, deep=True, array_type=VTK_INT)
    tree_array.SetName(tree_name)
    cotree_array = vtknp.numpy_to_vtk(cotree, deep=True, array_type=VTK_INT)
    cotree_array.SetName(cotree_name)
    info: dict[str, Any] = {
        "tree_build_algorithm": "multi_source_dijkstra_edge_length",
        "num_points": n_points,
        "num_cells": mesh.GetNumberOfCells(),
        "num_edges": edges_before_collapse,
        "num_edges_after_periodic_collapse": int(eligible_mids.size),
        "num_tree": int(tree.sum()),
        "num_cotree": int(cotree.sum()),
        "num_periodic_slave_marked": int(periodic_slave.sum()),
        "skipped_unknown_cells": sum(skipped.values()),
        "num_periodic_couplings": len(midside_masters),
        "num_components": components,
        "num_roots": roots,
        "num_tree_canonical": tree_canonical,
        "expected_tree_canonical": expected_tree,
        "tree_identity_ok": identity_ok,
    }
    return tree_array, cotree_array, info


def tree_gauge(
    mesh: vtkUnstructuredGrid,
    non_conducting_point_ids: Iterable[int] | None = None,
    start_surface_point_ids: Iterable[int] | None = None,
    periodic_point_coupling: Any = None,
    array_name: str = "_A",
    component: int = 2,
    tree_name: str = "TREE",
    cotree_name: str = "COTREE",
    periodic_slave_name: str = "PERIODIC_SLAVE",
) -> dict[str, Any]:
    """Gauge an edge field: build the tree-cotree arrays and fix the tree edges to zero in a Dirichlet array.

    The tree, cotree and periodic-slave arrays are added to the mesh (replacing older ones). In ``array_name`` the
    entries of the tree midside points are set to 0, except at periodic slaves, which the periodic constraint ties to
    their masters instead.

    Args:
        mesh: Mesh of quadratic cells, modified in place.
        non_conducting_point_ids: See :func:`build_tree_cotree`.
        start_surface_point_ids: See :func:`build_tree_cotree`.
        periodic_point_coupling: See :func:`build_tree_cotree`.
        array_name: Dirichlet point array of the edge field.
        component: Component of ``array_name`` holding the edge values (the solver reads component 2 of ``_F``).
        tree_name: Name of the tree point array.
        cotree_name: Name of the cotree point array.
        periodic_slave_name: Name of the periodic-slave point array.

    Returns:
        The diagnostics of :func:`build_tree_cotree`.

    Raises:
        MissingArrayError: If ``array_name`` or ``MidsideNodes`` is missing.
    """
    target = meshes.get_array(mesh, array_name, "node")
    if non_conducting_point_ids is not None:
        non_conducting_point_ids = sorted({int(i) for i in non_conducting_point_ids})
    if start_surface_point_ids is not None:
        start_surface_point_ids = sorted({int(i) for i in start_surface_point_ids})
    LOGGER.info(
        "Tree gauge: %s non-conducting points, %s start-surface points.",
        "selected cells'" if non_conducting_point_ids is None else len(non_conducting_point_ids),
        "selected" if start_surface_point_ids is None else len(start_surface_point_ids),
    )
    tree, cotree, info = build_tree_cotree(
        mesh,
        non_conducting_point_ids=non_conducting_point_ids,
        start_surface_point_ids=start_surface_point_ids,
        periodic_point_coupling=periodic_point_coupling,
        tree_name=tree_name,
        cotree_name=cotree_name,
        periodic_slave_name=periodic_slave_name,
    )
    point_data = mesh.GetPointData()
    for array in (tree, cotree):
        point_data.RemoveArray(array.GetName())
        point_data.AddArray(array)
    LOGGER.info("Tree-cotree gauge built: %s", info)

    fix = (vtknp.vtk_to_numpy(tree) == 1) & (meshes.array_values(mesh, periodic_slave_name, "node") == 0)
    values = vtknp.vtk_to_numpy(target)
    if values.ndim > 1:
        values[fix, component] = 0
    else:
        values[fix] = 0
    target.Modified()
    return info

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import numpy as np


@dataclass(frozen=True)
class CouplingForest:
    """Describes the parent -> child edges that couple a single CP set."""

    array_name: str
    cp_id: int
    pairs: np.ndarray
    roots: np.ndarray
    isolated: np.ndarray

    def has_edges(self) -> bool:
        """True when the forest contains at least one coupling edge."""
        return self.pairs.size > 0


class DOFCouplingHelper:
    """Construct coupling forests from CP_* point data arrays in a VTK mesh."""

    CP_PREFIX = "CP_"

    @staticmethod
    def build_coupling_forests(grid) -> List[CouplingForest]:
        """Return shallow spanning forests for every CP_* coupling set in *grid*.

        Parameters
        ----------
        grid : vtk.vtkUnstructuredGrid
            Mesh containing CP_* point data arrays.

        Returns
        -------
        List[CouplingForest]
            One forest per CP_* array / CP id combination. Forests are acyclic by
            construction and promote sparsity by pairing adjacent nodes first.
        """
        if grid is None:
            raise ValueError("grid must be a vtk.vtkUnstructuredGrid instance.")

        npoints = grid.GetNumberOfPoints()
        if npoints == 0:
            return []

        cp_sets = DOFCouplingHelper._collect_cp_sets(grid, npoints)
        if not cp_sets:
            return []

        adjacency = DOFCouplingHelper._build_node_adjacency(grid, npoints)
        coords = DOFCouplingHelper._collect_point_coordinates(grid, npoints)

        forests: List[CouplingForest] = []
        for array_name, groups in cp_sets.items():
            for cp_id, nodes in groups.items():
                forests.append(
                    DOFCouplingHelper._build_forest_for_nodes(
                        array_name=array_name,
                        cp_id=cp_id,
                        nodes=nodes,
                        adjacency=adjacency,
                        coords=coords,
                    )
                )
        return forests

    @staticmethod
    def print_coupling_forests(forests: Iterable[CouplingForest]) -> None:
        """Pretty-print the couplings produced by build_coupling_forests."""
        forests = list(forests)
        if not forests:
            print("No coupling forests computed.")
            return

        for forest in forests:
            header = f"[{forest.array_name}] cp={forest.cp_id}"
            if forest.roots:
                header += f" roots={', '.join(str(root) for root in forest.roots)}"
            print(header)

            if forest.pairs:
                for parent, child in forest.pairs:
                    print(f"  {parent} -> {child}")
            else:
                print("  (no coupling edges)")

            if forest.isolated:
                print(f"  isolated: {', '.join(str(idx) for idx in forest.isolated)}")

    @staticmethod
    def _build_forest_for_nodes(
        array_name: str,
        cp_id: int,
        nodes: Sequence[int],
        adjacency: Sequence[Set[int]],
        coords: Sequence[Tuple[float, float, float]],
    ) -> CouplingForest:
        nodes_set = set(nodes)
        remaining = set(nodes)
        edges: List[Tuple[int, int]] = []
        roots: List[int] = []
        isolated: List[int] = []

        while remaining:
            seed = DOFCouplingHelper._pick_seed_node(remaining, adjacency, nodes_set)
            remaining.remove(seed)

            component = DOFCouplingHelper._collect_component(seed, remaining, adjacency)
            component_set = set(component)

            if len(component_set) == 1:
                isolated.append(seed)
                continue

            root = DOFCouplingHelper._pick_root(component_set, adjacency)
            roots.append(root)
            component_edges, visited = DOFCouplingHelper._build_tree_edges(
                root, component_set, adjacency, coords
            )
            edges.extend(component_edges)

            leftovers = component_set - visited
            if leftovers:
                isolated.extend(sorted(leftovers))

        pairs_array = (
            np.array(edges, dtype=np.int64).reshape(-1, 2)
            if edges
            else np.zeros((0, 2), dtype=np.int64)
        )
        roots_array = np.array(roots, dtype=np.int64)
        isolated_array = np.array(isolated, dtype=np.int64)

        return CouplingForest(
            array_name=array_name,
            cp_id=cp_id,
            pairs=pairs_array,
            roots=roots_array,
            isolated=isolated_array,
        )

    @staticmethod
    def _collect_cp_sets(grid, npoints: int) -> Dict[str, Dict[int, List[int]]]:
        point_data = grid.GetPointData()
        cp_sets: Dict[str, Dict[int, List[int]]] = {}

        for idx in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(idx)
            if array is None:
                continue

            name = array.GetName()
            if not name or not name.startswith(DOFCouplingHelper.CP_PREFIX):
                continue

            suffix = name[len(DOFCouplingHelper.CP_PREFIX) :]
            if not suffix.isdigit():
                continue

            groups: Dict[int, List[int]] = defaultdict(list)
            for pid in range(npoints):
                value = DOFCouplingHelper._safe_get_scalar(array, pid)
                if value is None:
                    continue
                cp_id = int(round(value))
                if cp_id < 0:
                    continue
                groups[cp_id].append(pid)

            filtered = {cp_id: pts for cp_id, pts in groups.items() if len(pts) >= 2}
            if filtered:
                cp_sets[name] = filtered

        return cp_sets

    @staticmethod
    def _build_node_adjacency(grid, npoints: int) -> List[Set[int]]:
        adjacency: List[Set[int]] = [set() for _ in range(npoints)]

        for cell_id in range(grid.GetNumberOfCells()):
            cell = grid.GetCell(cell_id)
            point_ids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
            if len(point_ids) < 2:
                continue

            for i, u in enumerate(point_ids):
                neighbors = adjacency[u]
                for v in point_ids[i + 1 :]:
                    neighbors.add(v)
                    adjacency[v].add(u)

        return adjacency

    @staticmethod
    def _collect_point_coordinates(grid, npoints: int) -> List[Tuple[float, float, float]]:
        return [tuple(map(float, grid.GetPoint(pid))) for pid in range(npoints)]

    @staticmethod
    def _collect_component(
        seed: int,
        remaining: Set[int],
        adjacency: Sequence[Set[int]],
    ) -> List[int]:
        component = [seed]
        queue: deque[int] = deque([seed])

        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.append(neighbor)
                    queue.append(neighbor)

        return component

    @staticmethod
    def _build_tree_edges(
        root: int,
        component_set: Set[int],
        adjacency: Sequence[Set[int]],
        coords: Sequence[Tuple[float, float, float]],
    ) -> Tuple[List[Tuple[int, int]], Set[int]]:
        visited: Set[int] = {root}
        queue: deque[int] = deque([root])
        edges: List[Tuple[int, int]] = []

        while queue:
            parent = queue.popleft()
            neighbors = [
                child
                for child in adjacency[parent]
                if child in component_set and child not in visited
            ]
            neighbors.sort(
                key=lambda child: DOFCouplingHelper._squared_distance(
                    coords[parent], coords[child]
                )
            )

            for child in neighbors:
                visited.add(child)
                queue.append(child)
                edges.append((parent, child))

        return edges, visited

    @staticmethod
    def _pick_seed_node(
        candidates: Set[int],
        adjacency: Sequence[Set[int]],
        allowed: Set[int],
    ) -> int:
        best = None
        best_score = -1

        for node in candidates:
            score = len(adjacency[node] & allowed)
            if score > best_score or (score == best_score and (best is None or node < best)):
                best = node
                best_score = score

        if best is None:
            raise RuntimeError("Failed to pick a seed node for coupling forest construction.")

        return best

    @staticmethod
    def _pick_root(component_set: Set[int], adjacency: Sequence[Set[int]]) -> int:
        best = None
        best_score = -1

        for node in component_set:
            score = len(adjacency[node] & component_set)
            if score > best_score or (score == best_score and (best is None or node < best)):
                best = node
                best_score = score

        if best is None:
            raise RuntimeError("Failed to pick a root for coupling tree construction.")

        return best

    @staticmethod
    def _safe_get_scalar(array, index: int):
        tuple_value = array.GetTuple(index)
        if not tuple_value:
            return None
        return tuple_value[0]

    @staticmethod
    def _squared_distance(
        a: Tuple[float, float, float],
        b: Tuple[float, float, float],
    ) -> float:
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        return dx * dx + dy * dy + dz * dz

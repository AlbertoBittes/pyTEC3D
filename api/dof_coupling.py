from __future__ import annotations

from collections import defaultdict, deque
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    TYPE_CHECKING,
)
import warnings

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - import is for type hints only
    import pandas as pd
    import vtk

Adjacency = Sequence[Set[int]]
Coordinates = Sequence[Tuple[float, float, float]]
Forest = Dict[str, Any]


class DOFCoupling:
    """Construct coupling forests from CP_* point data arrays in a VTK mesh."""

    CP_PREFIX = "CP_"

    @staticmethod
    def build_dof_couplings(
        grid: "vtk.vtkUnstructuredGrid",
        coupled_variables: Optional["pd.DataFrame"] = None,
    ) -> List[Forest]:
        """
        Return shallow spanning forests for every CP_* coupling set in *grid*.

        Args:
            grid: Mesh containing CP_* point data arrays.
            coupled_variables: Optional DataFrame mapping CP ids to variable names.

        Returns:
            List[Forest]: One forest per CP_* array / CP id combination. Forests are
            acyclic by construction and promote sparsity by pairing adjacent nodes first.
            Each forest dictionary exposes the keys 'array_name' (string), 'cp_id' (int),
            'pairs' (shape [N, 2] int64 array), 'roots' (shape [R] int64 array), and
            'isolated' (shape [K] int64 array).

        Raises:
            ValueError: If `grid` is None.
        """
        if grid is None:
            raise ValueError("grid must be a vtk.vtkUnstructuredGrid instance.")

        npoints = grid.GetNumberOfPoints()
        if npoints == 0:
            return []

        cp_sets = DOFCoupling._collect_cp_sets(grid, npoints)
        if not cp_sets:
            return []

        coupled_variables_map = DOFCoupling._collect_coupled_variables(
            coupled_variables
        )
        adjacency = DOFCoupling._build_node_adjacency(grid, npoints)
        coords = DOFCoupling._collect_point_coordinates(grid, npoints)

        forests: List[Dict[str, Any]] = []
        for array_name, groups in cp_sets.items():
            for cp_id, nodes in groups.items():
                variables = coupled_variables_map.get(cp_id, [])
                if not variables:
                    warnings.warn(
                        f"CP set {cp_id} in '{array_name}' has no coupled variables defined.",
                        RuntimeWarning,
                    )
                forests.append(
                    DOFCoupling._build_forest_for_nodes(
                        array_name=array_name,
                        cp_id=cp_id,
                        nodes=nodes,
                        adjacency=adjacency,
                        coords=coords,
                        coupled_variables=variables,
                    )
                )
        return forests

    @staticmethod
    def print_dof_couplings(forests: Iterable[Forest]) -> None:
        """
        Pretty-print the couplings produced by build_dof_couplings.

        Args:
            forests: Iterable of forest dictionaries as returned by build_dof_couplings.

        Returns:
            None
        """
        forests = list(forests)
        if not forests:
            print("No coupling forests computed.")
            return

        for forest in forests:
            header = f"[{forest['array_name']}] cp={forest['cp_id']}"
            if forest["roots"].size:
                header += " roots=" + ", ".join(
                    str(root) for root in forest["roots"].tolist()
                )
            print(header)

            if "coupled_variables" in forest and forest["coupled_variables"].size:
                vars_str = ", ".join(
                    str(var) for var in forest["coupled_variables"].tolist()
                )
                print(f"  variables: {vars_str}")

            if forest["pairs"].size:
                for parent, child in forest["pairs"].tolist():
                    print(f"  {parent} -> {child}")
            else:
                print("  (no coupling edges)")

            if forest["isolated"].size:
                isolated_str = ", ".join(
                    str(idx) for idx in forest["isolated"].tolist()
                )
                print(f"  isolated: {isolated_str}")

    @staticmethod
    def _build_forest_for_nodes(
        array_name: str,
        cp_id: int,
        nodes: Sequence[int],
        adjacency: Adjacency,
        coords: Coordinates,
        coupled_variables: Sequence[str],
    ) -> Forest:
        """
        Build a coupling forest for a single CP id.

        Args:
            array_name: Name of the CP_* array the nodes came from.
            cp_id: Coupling set identifier extracted from the CP_* values.
            nodes: Node ids that belong to this coupling set.
            adjacency: Point connectivity graph for the mesh.
            coords: Point coordinates used for distance-based edge ordering.

        Returns:
            Forest: Dictionary containing 'array_name', 'cp_id', 'pairs', 'roots', and 'isolated'.
        """
        nodes_set = set(nodes)
        remaining = set(nodes)
        edges: List[Tuple[int, int]] = []
        roots: List[int] = []
        isolated: List[int] = []

        while remaining:
            seed = DOFCoupling._pick_seed_node(remaining, adjacency, nodes_set)
            remaining.remove(seed)

            component = DOFCoupling._collect_component(seed, remaining, adjacency)
            component_set = set(component)

            if len(component_set) == 1:
                isolated.append(seed)
                continue

            root = DOFCoupling._pick_root(component_set, adjacency)
            roots.append(root)
            component_edges, visited = DOFCoupling._build_tree_edges(
                root, component_set, adjacency, coords
            )
            edges.extend(component_edges)

            leftovers = component_set - visited
            if leftovers:
                isolated.extend(sorted(leftovers))

        pairs_array = (
            np.asarray(edges, dtype=np.int64).reshape(-1, 2)
            if edges
            else np.empty((0, 2), dtype=np.int64)
        )
        roots_array = np.asarray(roots, dtype=np.int64)
        isolated_array = np.asarray(isolated, dtype=np.int64)

        return {
            "array_name": array_name,
            "cp_id": int(cp_id),
            "pairs": pairs_array,
            "roots": roots_array,
            "isolated": isolated_array,
            "coupled_variables": np.asarray(coupled_variables, dtype=str),
        }

    @staticmethod
    def _collect_cp_sets(
        grid: "vtk.vtkUnstructuredGrid", npoints: int
    ) -> Dict[str, Dict[int, List[int]]]:
        """
        Collect CP_* point data arrays and group node ids by coupling id.

        Args:
            grid: Mesh containing CP_* point data arrays.
            npoints: Number of points in the mesh.

        Returns:
            Dict[str, Dict[int, List[int]]]: Mapping of array name -> cp_id -> node ids.
            Only cp_ids with at least two nodes are returned.
        """
        point_data = grid.GetPointData()
        cp_sets: Dict[str, Dict[int, List[int]]] = {}

        for idx in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(idx)
            if array is None:
                continue

            name = array.GetName()
            if not name or not name.startswith(DOFCoupling.CP_PREFIX):
                continue

            suffix = name[len(DOFCoupling.CP_PREFIX) :]
            if not suffix.isdigit():
                continue

            groups: Dict[int, List[int]] = defaultdict(list)
            for pid in range(npoints):
                value = DOFCoupling._safe_get_scalar(array, pid)
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
    def _collect_coupled_variables(
        coupled_variables: Optional["pd.DataFrame"],
    ) -> Dict[int, List[str]]:
        """
        Build a mapping of CP id to coupled variable names.

        Args:
            coupled_variables: DataFrame with 'cp_id' and 'variable' columns.

        Returns:
            Dict[int, List[str]]: CP id -> list of variable names.

        Raises:
            ValueError: If required columns are missing.
        """
        if coupled_variables is None:
            return {}

        if coupled_variables.empty:
            return {}

        required_columns = {"cp_id", "variable"}
        missing = required_columns - set(coupled_variables.columns)
        if missing:
            raise ValueError(
                f"Coupled variables dataframe is missing required columns: {sorted(missing)}."
            )

        mapping: Dict[int, List[str]] = defaultdict(list)
        for _, row in coupled_variables.iterrows():
            cp_value = row["cp_id"]
            variable = row["variable"]
            if cp_value is None or variable is None:
                continue
            try:
                cp_id = int(cp_value)
            except (TypeError, ValueError):
                continue
            var_name = str(variable).strip()
            if not var_name:
                continue
            mapping[cp_id].append(var_name)

        return mapping

    @staticmethod
    def _build_node_adjacency(
        grid: "vtk.vtkUnstructuredGrid", npoints: int
    ) -> Adjacency:
        """
        Build an undirected adjacency list of mesh nodes based on cell connectivity.

        Args:
            grid: Mesh whose cell connectivity defines adjacency.
            npoints: Number of points in the mesh.

        Returns:
            List[Set[int]]: Adjacency list indexed by point id.
        """
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
    def _collect_point_coordinates(
        grid: "vtk.vtkUnstructuredGrid", npoints: int
    ) -> Coordinates:
        """
        Extract point coordinates from the mesh for distance calculations.

        Args:
            grid: Mesh providing point coordinates.
            npoints: Number of points in the mesh.

        Returns:
            Coordinates: List of (x, y, z) tuples.
        """
        coords: List[Tuple[float, float, float]] = []
        for pid in range(npoints):
            x, y, z = grid.GetPoint(pid)
            coords.append((float(x), float(y), float(z)))
        return coords

    @staticmethod
    def _collect_component(
        seed: int,
        remaining: Set[int],
        adjacency: Adjacency,
    ) -> List[int]:
        """
        Collect all nodes reachable from a seed within the allowed set.

        Args:
            seed: Starting node id.
            remaining: Mutable set of candidate nodes (will be consumed).
            adjacency: Mesh adjacency list.

        Returns:
            List[int]: Nodes in the connected component including the seed.
        """
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
        adjacency: Adjacency,
        coords: Coordinates,
    ) -> Tuple[List[Tuple[int, int]], Set[int]]:
        """
        Build a shallow spanning tree for a component, favoring short edges.

        Args:
            root: Node id chosen as the tree root.
            component_set: Nodes belonging to the component.
            adjacency: Mesh adjacency list.
            coords: Point coordinates for distance ordering.

        Returns:
            Tuple[List[Tuple[int, int]], Set[int]]: Tree edges (parent, child) and visited nodes.
        """
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
                key=lambda child: DOFCoupling._squared_distance(
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
        adjacency: Adjacency,
        allowed: Set[int],
    ) -> int:
        """
        Select a seed node with the highest number of allowed neighbors.

        Args:
            candidates: Remaining candidate nodes to seed a component.
            adjacency: Mesh adjacency list.
            allowed: Nodes eligible to be connected.

        Returns:
            int: Chosen seed node id.

        Raises:
            RuntimeError: If no seed can be selected.
        """
        best = None
        best_score = -1

        for node in candidates:
            score = len(adjacency[node] & allowed)
            if score > best_score or (
                score == best_score and (best is None or node < best)
            ):
                best = node
                best_score = score

        if best is None:
            raise RuntimeError(
                "Failed to pick a seed node for coupling forest construction."
            )

        return best

    @staticmethod
    def _pick_root(component_set: Set[int], adjacency: Adjacency) -> int:
        """
        Choose a root node for a component, favoring highest degree then lowest id.

        Args:
            component_set: Nodes in the connected component.
            adjacency: Mesh adjacency list.

        Returns:
            int: Root node id.

        Raises:
            RuntimeError: If no root can be selected.
        """
        best = None
        best_score = -1

        for node in component_set:
            score = len(adjacency[node] & component_set)
            if score > best_score or (
                score == best_score and (best is None or node < best)
            ):
                best = node
                best_score = score

        if best is None:
            raise RuntimeError("Failed to pick a root for coupling tree construction.")

        return best

    @staticmethod
    def _safe_get_scalar(array: "vtk.vtkDataArray", index: int) -> Optional[float]:
        """
        Safely extract a scalar value from a VTK array tuple.

        Args:
            array: VTK data array.
            index: Tuple index to read.

        Returns:
            Optional[float]: Scalar value or None if tuple is empty.
        """
        tuple_value = array.GetTuple(index)
        if not tuple_value:
            return None
        return tuple_value[0]

    @staticmethod
    def _squared_distance(
        a: Tuple[float, float, float],
        b: Tuple[float, float, float],
    ) -> float:
        """
        Compute squared Euclidean distance between two points.

        Args:
            a: First point as (x, y, z).
            b: Second point as (x, y, z).

        Returns:
            float: Squared distance.
        """
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        return dx * dx + dy * dy + dz * dz

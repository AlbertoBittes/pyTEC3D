"""Node-to-node contact between bodies meshed with coincident interface nodes, through DOF couplings.

The two faces in contact carry distinct, coincident nodes (each body has its own). A **pair** is a node of the first
body (the master) and the coincident node of the second; ``n`` is the contact normal, from the first body to the
second. Three kinds of contact, from the constraint equations of the couplings (:doc:`couplings`):

* **bonded** (linear): the displacement of the pair is the same, all its components (``"U"``);
* **sliding** (linear, frictionless, no separation): only the normal component is tied (``"U[x]"``: the component
  along the axis of the nodal frame, ``node_angles``; the global axes when the frames are not rotated);
* **unilateral** (nonlinear, frictionless): the pair is tied in the normal direction while it is in compression and
  free when it would open, found by an active-set iteration (:func:`solve_unilateral`): solve with the active pairs
  tied; release those whose normal stress is tensile, activate those that interpenetrate; until the set does not
  change. The contact pressure of a pair is its normal stress, ``-n . sigma . n`` (the point array ``SIGMA``).

Couplings are fixed when a model is loaded: the active-set iteration reloads the model for each new set (a static
analysis: each iteration is a solve of its own). The CP ids of the pairs are taken above the largest id of the mesh,
in ``array_name`` (a CP_k array, created if missing).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import vtk
from vtkmodules.util import numpy_support as vtknp

from . import mesh as meshes
from ._logging import get_logger

LOGGER = get_logger("pytec3d.contact")

AXES = {"x": 0, "y": 1, "z": 2}


def interface_pairs(grid, first_ids, second_ids, tolerance: float = 1e-9) -> np.ndarray:
    """Pairs (first, second) of coincident nodes: each node of ``first_ids`` with the node of ``second_ids`` at its
    position (within ``tolerance``); shape (n, 2). Nodes without a match are skipped."""
    points = vtknp.vtk_to_numpy(grid.GetPoints().GetData())
    first_ids, second_ids = np.asarray(first_ids, dtype=np.int64), np.asarray(second_ids, dtype=np.int64)
    # Points on a grid of cell size tolerance: a match is in the same or a neighbouring cell.
    key = lambda p: tuple(np.floor(p / tolerance).astype(np.int64))
    cells: dict[tuple, list[int]] = {}
    for j in second_ids:
        cells.setdefault(key(points[j]), []).append(int(j))
    out = []
    for i in first_ids:
        base = np.floor(points[i] / tolerance).astype(np.int64)
        best, best_d = None, tolerance
        for offset in np.ndindex(3, 3, 3):
            for j in cells.get(tuple(base + np.array(offset) - 1), ()):
                d = float(np.linalg.norm(points[j] - points[i]))
                if d <= best_d and j != i:
                    best, best_d = j, d
        if best is not None:
            out.append((int(i), best))
    return np.array(out, dtype=np.int64).reshape(-1, 2)


def _cp_array(grid, array_name: str) -> np.ndarray:
    if grid.GetPointData().GetArray(array_name) is None:
        meshes.add_array(grid, array_name, dtype=vtk.VTK_INT)
    return meshes.array_values(grid, array_name, "node")


def _next_cp_id(grid) -> int:
    largest = 0
    data = grid.GetPointData()
    for i in range(data.GetNumberOfArrays()):
        name = data.GetArrayName(i)
        if name and name.startswith("CP_"):
            values = vtknp.vtk_to_numpy(data.GetArray(i))
            if values.size:
                largest = max(largest, int(values.max()))
    return largest + 1


def couple_pairs(model, pairs: np.ndarray, fields: list[str], array_name: str = "CP_1") -> dict[int, list[str]]:
    """Give each pair its own CP id in ``array_name`` and couple ``fields`` on it (``model.add_dof_couplings``);
    returns {cp_id: fields}. The model must be unlocked (no solver session)."""
    grid = model.mesh
    cp = _cp_array(grid, array_name)
    start = _next_cp_id(grid)
    couplings = {}
    for k, (a, b) in enumerate(np.asarray(pairs, dtype=np.int64)):
        cp[a] = cp[b] = start + k
        couplings[start + k] = list(fields)
    if couplings:
        model.add_dof_couplings(couplings)
    grid.GetPointData().GetArray(array_name).Modified()
    return couplings


def clear_pairs(model, cp_ids, array_name: str = "CP_1") -> None:
    """Remove the couplings ``cp_ids``: their nodes get 0 in ``array_name``, their entries leave ``dof_couplings``."""
    cp = _cp_array(model.mesh, array_name)
    ids = set(int(i) for i in cp_ids)
    cp[np.isin(cp, list(ids))] = 0
    couplings = model.data.get("dof_couplings", {})
    for key in list(couplings):
        if int(key) in ids:
            del couplings[key]
    model.data = {**model.data, "dof_couplings": couplings}


def bonded(model, pairs: np.ndarray, field: str = "U", array_name: str = "CP_1") -> dict[int, list[str]]:
    """Bonded contact: the pairs share all the components of ``field``."""
    return couple_pairs(model, pairs, [field], array_name)


def sliding(model, pairs: np.ndarray, axis: str = "z", field: str = "U", array_name: str = "CP_1") -> dict[int, list[str]]:
    """Frictionless contact without separation: the pairs share the component ``axis`` of ``field`` (in the nodal
    frame of the first node)."""
    return couple_pairs(model, pairs, [f"{field}[{axis}]"], array_name)


def solve_unilateral(solver, pairs: np.ndarray, axis: str = "z", sign: float = 1.0, field: str = "U",
                     array_name: str = "CP_1", max_iterations: int = 30, tolerance: float = 0.0,
                     active: np.ndarray | None = None) -> dict[str, Any]:
    """Frictionless unilateral contact of the pairs, normal ``n = sign e_axis`` (from the first body to the second),
    by an active-set iteration around static solves (the solver's model must not be loaded: each set is loaded).

    ``active``: the initial set (default: all the pairs, closed). A pair is released when its normal stress at the
    first node is tensile (``n . sigma . n > tolerance``), activated when the second node moves into the first body
    (``(u2 - u1) . n < -tolerance``, the length unit of the model). Returns ``active`` (bool per pair), ``pressure``
    (``-n . sigma . n`` at the first nodes, 0 where open), ``gap`` ((u2 - u1) . n), ``iterations``; raises
    RuntimeError when the set does not settle in ``max_iterations``."""
    model = solver.model
    grid = model.mesh
    pairs = np.asarray(pairs, dtype=np.int64)
    k = AXES[axis]
    active = np.ones(len(pairs), dtype=bool) if active is None else np.asarray(active, dtype=bool).copy()
    cp_ids: list[int] = []
    for iteration in range(1, max_iterations + 1):
        solver.close_session()
        clear_pairs(model, cp_ids, array_name)
        cp_ids = list(sliding(model, pairs[active], axis, field, array_name)) if active.any() else []
        solver.load_model_into_solver()
        solver.solve()
        u = meshes.array_values(grid, field, "node")
        sigma = meshes.array_values(grid, "SIGMA", "node").reshape(-1, 3, 3)
        normal_stress = sigma[pairs[:, 0], k, k]
        gap = sign * (u[pairs[:, 1], k] - u[pairs[:, 0], k])
        release = active & (normal_stress > tolerance)
        close = ~active & (gap < -tolerance)
        LOGGER.info("Contact iteration %d: %d active pair(s), %d to release, %d to close.", iteration,
                    int(active.sum()), int(release.sum()), int(close.sum()))
        if not release.any() and not close.any():
            pressure = np.where(active, -normal_stress, 0.0)
            return {"active": active, "pressure": pressure, "gap": gap, "iterations": iteration}
        active = (active & ~release) | close
    raise RuntimeError(f"Unilateral contact: the active set did not settle in {max_iterations} iterations.")

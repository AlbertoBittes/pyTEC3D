"""Build the meshes of the scalar-potential examples (core_3D_phi, ring_3D_tphi, rings_3D_tphi, coil_core_3D_tphi, wire_2D_tphi,
bulk_2D_tphi_axisym): structured quadratic grids (HEX20, QUAD8) with their TYPE and MATERIAL cell arrays, written as
<case>/<case>.vtu next to this folder. Only numpy and vtk.

    cd examples/files && python formulation_meshes.py          # all of them
    python formulation_meshes.py wire_2D_tphi                  # one of them

The examples ship these files; rerun this script only to change a geometry (the macros take the dimensions from their
own constants: keep them in step).
"""

from __future__ import annotations

import itertools
import os
import sys

import numpy as np
import vtk
from vtkmodules.util import numpy_support as vtknp

HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7))
HEX_CORNERS = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
QUAD_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
QUAD_CORNERS = ((0, 0), (1, 0), (1, 1), (0, 1))


def tensor_grid(axes: list[np.ndarray]) -> dict[str, np.ndarray]:
    """A quadratic grid on the tensor product of the coordinate lists ``axes`` (2 lists: QUAD8, 3: HEX20; the midside
    points halfway). Returns ``points`` (n, 3), ``cells`` (VTK order) and ``centroids``."""
    dim = len(axes)
    corners, edges = (HEX_CORNERS, HEX_EDGES) if dim == 3 else (QUAD_CORNERS, QUAD_EDGES)
    half = []  # coordinates at half indices: 2 i the grid lines, 2 i + 1 the midpoints
    for a in axes:
        h = np.empty(2 * len(a) - 1)
        h[0::2], h[1::2] = a, 0.5 * (a[1:] + a[:-1])
        half.append(h)
    index: dict[tuple[int, ...], int] = {}
    points: list[tuple[float, ...]] = []

    def node(key: tuple[int, ...]) -> int:
        if key not in index:
            index[key] = len(points)
            points.append(tuple(half[d][k] for d, k in enumerate(key)))
        return index[key]

    cells = []
    for cell in itertools.product(*(range(len(a) - 1) for a in axes)):
        row = [node(tuple(2 * c + 2 * o for c, o in zip(cell, offset))) for offset in corners]
        for p, q in edges:
            row.append(node(tuple(2 * c + corners[p][d] + corners[q][d] for d, c in enumerate(cell))))
        cells.append(row)
    xyz = np.zeros((len(points), 3))
    xyz[:, :dim] = np.array(points)
    cells = np.array(cells)
    return {"points": xyz, "cells": cells, "centroids": xyz[cells[:, :len(corners)]].mean(axis=1)}


def graded(inner: float, outer: float, h: float, cells: int) -> np.ndarray:
    """Coordinates from ``inner`` to ``outer`` (exclusive of inner) in ``cells`` cells growing geometrically from h."""
    lo, hi = 1.0, 2.0
    for _ in range(200):  # growth ratio q with h (q^cells - 1) / (q - 1) = outer - inner
        q = 0.5 * (lo + hi)
        lo, hi = (q, hi) if h * (q**cells - 1.0) / (q - 1.0) < outer - inner else (lo, q)
    return inner + h * np.cumsum(q ** np.arange(cells))


def ogrid_disk(m: int, n_in: int, n_out: int, length: float, radius: float, core: float) -> dict[str, np.ndarray]:
    """A QUAD8 O-grid of the disk of radius ``radius`` in the square [-L, L]^2: a square core of half side core * radius
    (m cells per side), n_in layers to the circle, n_out to the square. ``inside``: the disk cells."""
    keys: dict[tuple[int, int, int], int] = {}
    cells, inside = [], []

    def node(i: int, j: int, level: int) -> int:
        return keys.setdefault((i, j, level), len(keys))

    def add(params: list[tuple[int, int, int]], is_inside: bool) -> None:
        mids = [tuple((params[a][k] + params[b][k]) // 2 for k in range(3)) for a, b in QUAD_EDGES]
        cells.append([node(*p) for p in params + mids])
        inside.append(is_inside)

    for ci, cj in itertools.product(range(-m // 2, m // 2), repeat=2):
        add([(2 * ci + 2 * u, 2 * cj + 2 * v, 0) for u, v in QUAD_CORNERS], True)
    for axis in range(2):
        for side in (-1, 1):
            for c in range(-m // 2, m // 2):
                for layer in range(n_in + n_out):
                    params = []
                    for u, v in QUAD_CORNERS:
                        ij = [0, 0]
                        ij[axis], ij[1 - axis] = side * m, 2 * c + 2 * u
                        params.append((ij[0], ij[1], 2 * (layer + v)))
                    add(params, layer < n_in)
    key_list = sorted(keys, key=keys.get)
    ij = np.array([k[:2] for k in key_list], float) / m
    level = np.array([k[2] for k in key_list]) / 2.0
    norm = np.linalg.norm(ij, axis=1)
    direction = np.divide(ij, norm[:, None], out=np.zeros_like(ij), where=norm[:, None] > 0)
    t_in, t_out = np.clip(level / n_in, 0, 1)[:, None], np.clip((level - n_in) / n_out, 0, 1)[:, None]
    mapped = np.where(level[:, None] <= n_in, (1 - t_in) * core * radius * ij + t_in * radius * direction,
                      (1 - t_out) * radius * direction + t_out * length * ij)
    mapped[level == 0] = core * radius * ij[level == 0]
    xyz = np.zeros((len(key_list), 3))
    xyz[:, :2] = mapped
    cells = np.array(cells)
    p = xyz[cells[:, :4]]
    u, v = p[:, 1, :2] - p[:, 0, :2], p[:, 3, :2] - p[:, 0, :2]
    flip = u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0] < 0
    cells[flip] = cells[flip][:, [0, 3, 2, 1, 7, 6, 5, 4]]
    return {"points": xyz, "cells": cells, "centroids": xyz[cells[:, :4]].mean(axis=1), "inside": np.array(inside)}


def write(case: str, msh: dict[str, np.ndarray], types: np.ndarray) -> None:
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetData(vtknp.numpy_to_vtk(np.ascontiguousarray(msh["points"]), deep=True))
    grid.SetPoints(points)
    cell_type = vtk.VTK_QUADRATIC_HEXAHEDRON if msh["cells"].shape[1] == 20 else vtk.VTK_QUADRATIC_QUAD
    for row in msh["cells"]:
        ids = vtk.vtkIdList()
        for i in row:
            ids.InsertNextId(int(i))
        grid.InsertNextCell(cell_type, ids)
    for name in ("TYPE", "MATERIAL"):
        array = vtknp.numpy_to_vtk(np.ascontiguousarray(types, dtype=np.int32), deep=True, array_type=vtk.VTK_INT)
        array.SetName(name)
        grid.GetCellData().AddArray(array)
    folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", case)
    os.makedirs(folder, exist_ok=True)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(os.path.join(folder, case + ".vtu"))
    writer.SetInputData(grid)
    writer.SetDataModeToBinary()
    writer.SetCompressorTypeToZLib()
    writer.Write()
    counts = {int(t): int(c) for t, c in zip(*np.unique(types, return_counts=True))}
    print(f"{case}: {len(msh['points'])} points, {len(types)} cells, cells per TYPE {counts}")


def box_grid(lx: float, lz: float, d: float) -> dict[str, np.ndarray]:
    """The Cartesian grid |x|, |y| < lx, |z| < lz of cubes of side d."""
    xy = np.linspace(-lx, lx, round(2 * lx / d) + 1)
    return tensor_grid([xy, xy, np.linspace(-lz, lz, round(2 * lz / d) + 1)])


# ---------------------------------------------------------------------------------------------------------- cases
# The square ring / frame / winding a < max(|x|, |y|) < b, |z| < h on a 0.25 m grid, the box |x|, |y| < 2, |z| < 1.5 m.
LX, LZ, D = 2.0, 1.5, 0.25


def core_3D_phi() -> None:
    """A closed square iron frame (TYPE 2), a < max(|x|, |y|) < b, |z| < h (a = 0.5, b = 1, h = 0.5), with a coil
    (TYPE 3) of thickness t = 0.25 wound on the leg 0.5 < x < 1 over |y| < 0.5; air TYPE 1."""
    msh = box_grid(LX, LZ, D)
    x, y, z = msh["centroids"].T
    r = np.maximum(np.abs(x), np.abs(y))
    iron = (r > 0.5) & (r < 1.0) & (np.abs(z) < 0.5)
    coil = (np.abs(x - 0.75) < 0.5) & (np.abs(z) < 0.75) & (np.abs(y) < 0.5) & ~iron
    write("core_3D_phi", msh, np.where(iron, 2, np.where(coil, 3, 1)))


def ring_3D_tphi() -> None:
    """A square conducting ring (TYPE 2), 0.5 < max(|x|, |y|) < 1, |z| < 0.25, and the one-cell cut layer spanning its
    hole, max(|x|, |y|) < 0.5, 0 < z < 0.25 (TYPE 3); air TYPE 1."""
    msh = box_grid(LX, LZ, D)
    x, y, z = msh["centroids"].T
    r = np.maximum(np.abs(x), np.abs(y))
    ring = (r > 0.5) & (r < 1.0) & (np.abs(z) < 0.25)
    layer = (r < 0.5) & (z > 0.0) & (z < D)
    write("ring_3D_tphi", msh, np.where(ring, 2, np.where(layer, 3, 1)))


def rings_3D_tphi() -> None:
    """Two square rings stacked along z, 0.5 < max(|x|, |y|) < 1, |z -+ 0.5| < 0.25: the lower ring TYPE 2, the upper
    TYPE 3; their cut layers, max(|x|, |y|) < 0.5 and zc < z < zc + 0.25 (zc = -0.5, +0.5): TYPE 4 (lower), 5 (upper);
    air TYPE 1."""
    msh = box_grid(LX, LZ, D)
    x, y, z = msh["centroids"].T
    r = np.maximum(np.abs(x), np.abs(y))
    types = np.ones(len(r), dtype=int)
    for k, zc in ((0, -0.5), (1, 0.5)):
        types[(0.5 < r) & (r < 1.0) & (np.abs(z - zc) < 0.25)] = 2 + k
        types[(r < 0.5) & (zc < z) & (z < zc + D)] = 4 + k
    write("rings_3D_tphi", msh, types)


def coil_core_3D_tphi() -> None:
    """A square solenoid (winding TYPE 2), 0.5 < max(|x|, |y|) < 1, |z| < 0.5, around a conducting core (TYPE 3),
    max(|x|, |y|) < 0.25, |z| < 0.25; the rest of the hole, max(|x|, |y|) < 0.5, |z| < 0.5, TYPE 4; air TYPE 1."""
    msh = box_grid(LX, LZ, D)
    x, y, z = msh["centroids"].T
    r = np.maximum(np.abs(x), np.abs(y))
    winding = (r > 0.5) & (r < 1.0) & (np.abs(z) < 0.5)
    core = (r < 0.25) & (np.abs(z) < 0.25)
    hole = (r < 0.5) & (np.abs(z) < 0.5) & ~core
    write("coil_core_3D_tphi", msh, np.select([winding, core, hole], [2, 3, 4], 1))


def wire_2D_tphi() -> None:
    """A round wire of radius 1 (TYPE 2) in the square [-2, 2]^2, O-grid; the cut, a one-cell strip of air cells
    along the ray y = 0+, x > 1 (TYPE 3); air TYPE 1."""
    msh = ogrid_disk(12, 12, 6, 2.0, 1.0, 0.3)
    xyz, cells = msh["points"], msh["cells"]
    on_ray = np.isclose(xyz[:, 1], 0.0) & (xyz[:, 0] > 0.0)
    strip = ~msh["inside"] & (msh["centroids"][:, 1] > 0.0) & on_ray[cells[:, :4]].any(axis=1)
    write("wire_2D_tphi", msh, np.where(msh["inside"], 2, np.where(strip, 3, 1)))


def bulk_2D_tphi_axisym() -> None:
    """Axisymmetric (x = r, y = z) half section of a bulk disk (TYPE 2) of radius 10 mm and thickness 10 mm (|z| < 5
    mm), 0.5 mm cells, in air (TYPE 1) up to r = 60 mm, |z| < 60 mm (graded cells)."""
    h = 0.5e-3
    r = np.concatenate((np.linspace(0.0, 10e-3, 21), graded(10e-3, 60e-3, h, 16)))
    zp = np.concatenate((np.linspace(0.0, 5e-3, 11), graded(5e-3, 60e-3, h, 18)))
    z = np.concatenate((-zp[:0:-1], zp))
    msh = tensor_grid([r, z])
    c = msh["centroids"]
    write("bulk_2D_tphi_axisym", msh, np.where((c[:, 0] < 10e-3) & (np.abs(c[:, 1]) < 5e-3), 2, 1))


CASES = {f.__name__: f for f in (core_3D_phi, ring_3D_tphi, rings_3D_tphi, coil_core_3D_tphi, wire_2D_tphi, bulk_2D_tphi_axisym)}

if __name__ == "__main__":
    for name in sys.argv[1:] or CASES:
        CASES[name]()

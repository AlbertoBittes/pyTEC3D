"""End-to-end test of an installed pyTEC3D (API and solver): a 2D magnetostatic model solved through pytec3d.Solver.

A 2 m x 2 m air box (40 x 40 QUAD4 cells) with a 0.4 m x 0.4 m stranded conductor carrying Js = 1 MA/m^2 along z, and
A_z = 0 on the outer boundary: the same model as solver/examples/square_conductor_2d.cpp, whose results it checks. The
load is then doubled in place, and the solution must double (the solver reads the arrays at every solve).

Run: python test_square_conductor.py   (or mpiexec -n 2 python test_square_conductor.py); exit code 0 when it passes.
"""

import os
import sys
import tempfile

import numpy as np
import vtk
from vtkmodules.util import numpy_support as vtknp

from pytec3d import Solver, mesh

N = 40
EXPECTED_A_CENTRE = 6.579678e-02  # Wb/m
EXPECTED_B_MAX = 1.241108e-01  # T


def build_grid() -> tuple[vtk.vtkUnstructuredGrid, np.ndarray]:
    h = 2.0 / N
    xy = np.array([[-1 + i * h, -1 + j * h, 0.0] for j in range(N + 1) for i in range(N + 1)])
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetData(vtknp.numpy_to_vtk(xy, deep=True))
    grid.SetPoints(points)
    for j in range(N):
        for i in range(N):
            p = j * (N + 1) + i
            ids = vtk.vtkIdList()
            for q in (p, p + 1, p + N + 2, p + N + 1):
                ids.InsertNextId(q)
            grid.InsertNextCell(vtk.VTK_QUAD, ids)
    return grid, xy


def main() -> int:
    os.chdir(tempfile.mkdtemp(prefix="pytec3d_test_"))  # the solver writes its element matrices in ./emat
    grid, xy = build_grid()
    solver = Solver.from_grid(grid, mesh_name="square_conductor")
    identity = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    solver.model.data = {
        "element_types": {
            "1": {"element_type": "AV2D_nodal", "option1": "NON_CONDUCTING", "option2": 1, "option3": 0},
            "2": {"element_type": "AV2D_nodal", "option1": "STRANDED_CONDUCTOR", "option2": 1, "option3": 0},
        },
        "materials": {"1": {"mur": identity}},
        "constants": {"1": {"thickness": 1.0}},
    }
    solver.model.update_grid()
    centroids = mesh.cell_centroids(grid)
    inside = (np.abs(centroids[:, 0]) < 0.2) & (np.abs(centroids[:, 1]) < 0.2)
    mesh.array_values(grid, "TYPE", "elem")[:] = np.where(inside, 2, 1)
    js = mesh.array_values(grid, "Js", "elem")
    js[inside, 2] = 1e6
    boundary = np.isclose(np.abs(xy[:, 0]), 1) | np.isclose(np.abs(xy[:, 1]), 1)
    mesh.array_values(grid, "_A", "node")[boundary, 2] = 0.0

    solver.time_step.update(time=1.0, transient_effects=False)
    quiet = ("-ksp_monitor", "-ksp_view", "-ksp_converged_reason", "-log_view", "-ksp_monitor_true_residual",
             "-ksp_monitor_solution", "-mat_mumps_icntl_4")
    solver.solver_settings = {"petsc_options": dict.fromkeys(quiet)}
    solver.load_model_into_solver()
    solver.solve()

    a = mesh.array_values(grid, "A", "node")
    b = mesh.array_values(grid, "B", "node")
    centre = (N // 2) * (N + 1) + N // 2
    a_centre = float(a[centre, 2])
    b_max = float(np.hypot(b[:, 0], b[:, 1]).max())

    js[inside, 2] *= 2.0  # in place: the solver views this memory
    solver.solve()
    ratio = float(a[centre, 2]) / a_centre
    solver.close()

    failures = []
    if not np.isclose(a_centre, EXPECTED_A_CENTRE, rtol=1e-5):
        failures.append(f"A_z at the centre {a_centre:.6e}, expected {EXPECTED_A_CENTRE:.6e}")
    if not np.isclose(b_max, EXPECTED_B_MAX, rtol=1e-5):
        failures.append(f"max |B| {b_max:.6e}, expected {EXPECTED_B_MAX:.6e}")
    if not np.isclose(ratio, 2.0, rtol=1e-9):
        failures.append(f"A_z ratio with the doubled load {ratio:.12f}, expected 2")
    if solver.rank == 0:
        print(f"A_z centre {a_centre:.6e}  max |B| {b_max:.6e}  ratio {ratio:.12f}")
        print("FAILED:\n  " + "\n  ".join(failures) if failures else "PASSED")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

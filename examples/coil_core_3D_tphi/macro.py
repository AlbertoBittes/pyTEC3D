# =====================================================================================================================
# coil_core_3D_tphi - T-phi (TPHI3D_edge): a coil given as FIXED T around a conducting core
# =====================================================================================================================
#
# What this example shows
# -----------------------
# TPHI3D_edge has no source current and no circuit: a current-driven coil is T fixed on its edges. Any field T0 with
# curl T0 = J (the current density of the winding at unit current) that is zero outside the winding and the hole it
# encloses will do; its circulations along the edges, times I(t), are imposed on every edge of those cells. PHI does
# the rest. Here a square solenoid, current ramped then held, induces eddy currents in a conducting core block placed
# in its hole.
#   * T0 = f e_z, f = J (b - max(|x|, |y|)) in the winding (a < max(|x|, |y|) < b), J (b - a) in the hole column
#     (max(|x|, |y|) < a, over the height of the winding: "the disk filling the hole", here as thick as the coil),
#     0 elsewhere: curl T0 = J around z in the winding, 0 in the hole (T0 uniform there);
#   * a conductor in the hole has its own T, which cannot be fixed: T0 is moved out of the core by a gradient,
#     T0' = T0 - grad(g), g = J (b - a) z chi, chi = 1 on the nodes of the core cells: on the edges grad(g) is exact
#     (g(end) - g(start)), so T0' is zero on every edge of the core and carries the same current;
#   * the edge values: tree_gauge.edges() (each edge in the solver's orientation) and tree_gauge.edge_circulations();
#   * the winding is SOLID_CONDUCTOR (its T fixed: no eddy current there, its rho unused), the rest of the hole
#     CUT_LAYER (T without conduction), the core SOLID_CONDUCTOR with the usual gauges;
#   * the mean flux through a turn, Psi = int T0 . B dV per ampere-turn (the induced voltage is N dPsi / dt).
#
# Geometry and mesh (coil_core_3D_tphi.vtu, built by ../files/formulation_meshes.py): the box |x|, |y| < 2 m,
#   |z| < 1.5 m, 3072 HEX20 cells of 0.25 m. Winding TYPE / MATERIAL 2: 0.5 < max(|x|, |y|) < 1, |z| < 0.5; core
#   TYPE / MATERIAL 3: max(|x|, |y|) < 0.25, |z| < 0.25; the rest of the hole TYPE / MATERIAL 4: max(|x|, |y|) < 0.5,
#   |z| < 0.5; air TYPE 1.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order)
#   * 1: TPHI3D_edge / NON_CONDUCTING  / 1   air
#   * 2: TPHI3D_edge / SOLID_CONDUCTOR / 1   winding (T fixed)
#   * 3: TPHI3D_edge / SOLID_CONDUCTOR / 1   core (T free: eddy currents)
#   * 4: TPHI3D_edge / CUT_LAYER       / 1   the rest of the hole (T fixed, no conduction)
#   Materials: 1 air; 2 rho = 1 (not used); 3 rho = 1e-6 Ohm m; 4 air.
#
# Boundary conditions: the box is a flux wall (PHI free, one node fixed); core: n x T = 0 (T0' = 0 on its surface
#   edges), PHI = 0 at its inner nodes.
# Loads and time history: NI = 0.5 MA-turns (J = 1 MA/m^2) ramped from 0 over 0.1 s, then held; backward Euler,
#   dt = 0.02 s, 15 steps (the diffusion time of the core, mu0 c^2 / rho, is about 0.08 s).
# Results: output/coil_core_3D_tphi_out_<n>.vtu, coil_core_3D_tphi_steps.csv, coil_core_3D_tphi.log; printed: B_z at
#   the centre, the largest eddy current density in the core and Psi.
#
# How to run:  cd examples/coil_core_3D_tphi && python macro.py
# Documentation: "Element library - TPHI3D_edge (a coil as T fixed on its edges)", "Theory - scalar potentials".
# =====================================================================================================================

import os

import numpy as np

from pytec3d import Solver, mesh, selection, tree_gauge

A, B, H = 0.5, 1.0, 0.5                  # winding: a < max(|x|, |y|) < b, |z| < h
NI, T_RAMP, DT, STEPS = 5e5, 0.1, 0.02, 15
J_UNIT = 1.0 / ((B - A) * 2 * H)          # current density of the winding per ampere-turn (A/m^2)

# 1. Solver, mesh and model data. Everything below is set by selections: select points or cells, then
#    set_value_on_selection() writes one value, or one value per selected id (in increasing id order).
pytec3d = Solver()
filename = "coil_core_3D_tphi"
pytec3d.read_model(os.path.join(os.getcwd(), filename + ".vtu"))
pytec3d.model.read_json("model_data.json")
grid = pytec3d.model.mesh
xyz = mesh.points(grid)


# 2. T0 per ampere-turn and its edge values on the edges of the winding and of the hole (TYPE 2 and 4): the edges of
#    the selected cells, each with its midside point (sorted: the order of the selection) and its end points in the
#    solver's orientation.
def t0(p):
    r = np.clip(np.maximum(np.abs(p[:, 0]), np.abs(p[:, 1])), A, B)
    f = J_UNIT * (B - r) * (np.abs(p[:, 2]) < H)
    return np.column_stack((np.zeros(len(p)), np.zeros(len(p)), f))


selection.select(grid, "TYPE", [2, 4], location="elem")
mids, start, end = tree_gauge.edges(grid, selection.selected_ids(grid, "elem"))
values = tree_gauge.edge_circulations(grid, start, end, t0)
#    Out of the core: minus the edge differences of g = J (b - a) z chi, chi = 1 on the points of the core cells.
selection.select(grid, "TYPE", 3, location="elem")
selection.nsle(grid)
chi = selection.selection_mask(grid, "node").astype(float)
g = J_UNIT * (B - A) * xyz[:, 2] * chi
values -= g[end] - g[start]
core_mids, _, _ = tree_gauge.edges(grid, selection.selected_ids(grid, "elem"))
print("largest |T0'| on the core edges:", np.abs(values[np.isin(mids, core_mids)]).max())

# 3. Gauges of the core: n x T = 0 on its surface (its surface edges are hole edges: T0' = 0 there, fixed in the
#    loop), PHI = 0 at its inner corner points (the points of the core cells minus those of the other cells).
selection.select(grid, "TYPE", 3, location="elem")
selection.nsle(grid)
selection.select(grid, "TYPE", [1, 2, 4], location="elem")
selection.nsle(grid, mode="u")
selection.select(grid, "MidsideNodes", 0, mode="r")
selection.set_value_on_selection(grid, "_PHI", 0.0)

# 4. The box: a flux wall; PHI fixed on one of its points.
selection.select_all(grid)
selection.nsel_ext(grid)
selection.select(grid, "ID", int(selection.selected_ids(grid)[0]))
selection.set_value_on_selection(grid, "_PHI", 0.0)

# 5. Load the model; backward Euler; at each step _T = NI(t) T0' on the fixed edges (one value per selected midside
#    point: the selection of the midside points of the coil edges, kept through the loop).
pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 1.0
centroids, volumes = mesh.cell_centroids(grid), mesh.cell_sizes(grid)
t0_cells = t0(centroids)
core = mesh.array_values(grid, "TYPE", "elem") == 3
centre = int(np.argmin(np.linalg.norm(xyz, axis=1)))
selection.select(grid, "ID", mids)
for step in range(1, STEPS + 1):
    t = step * DT
    ni = NI * min(t / T_RAMP, 1.0)
    selection.set_value_on_selection(grid, "_T", ni * values, 2)
    pytec3d.time_step["time"] = t
    pytec3d.solve()
    pytec3d.write_output()
    b = mesh.array_values(grid, "B", "elem")
    jt = mesh.array_values(grid, "Jt", "elem")
    psi = float(np.einsum("ij,ij->i", t0_cells, b) @ volumes)
    print(f"t = {t:.2f} s: NI = {ni / 1e3:5.0f} kA, B_z(0) = {mesh.array_values(grid, 'B', 'node')[centre, 2]:.4f} T, "
          f"max |J| core = {np.linalg.norm(jt[core], axis=1).max() / 1e3:7.1f} kA/m^2, Psi = {psi:.4f} Wb")

pytec3d.save_model()

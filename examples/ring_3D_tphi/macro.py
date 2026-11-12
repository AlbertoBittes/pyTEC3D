# =====================================================================================================================
# ring_3D_tphi - T-phi (TPHI3D_edge): a conducting ring in a field switched on, its net current through a CUT
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The T-phi formulation, H = T - grad(PHI): T (edge, NEDELEC_ONE) only in the conductors and the cut layers, PHI (nodal)
# everywhere. A uniform axial field is switched on at t = 0 around a square conducting ring: the induced ring current
# then decays with the L / R time of the ring.
#   * the gauges of T-phi: n x T = 0 on the conductor surface (_T = 0 on its edges), PHI = 0 at the nodes inside the
#     conductor (first-order PHI: its gradients are in the edge space of T);
#   * the applied field: the total potential on the box, PHI = -H0 z;
#   * the CUT: outside the conductors H = -grad(PHI) has no circulation, so a loop around the ring could not see its
#     current. A one-cell layer of CUT_LAYER cells spanning the hole (0 < z < 0.25 inside the ring) carries the jump:
#     T = I grad(c), c from 0 (bottom face) to 1 (top face). Its edges in the two planes are zero; its edges crossing
#     the layer (vertical) are coupled on one DOF, the ring current I, unknown here (an induced current). CP ids 2 and
#     3 on the corner nodes of the two faces give the solver the orientation of each crossing edge.
#   * the ring current from the cell values of Jt over a section, and the E output.
#
# Geometry and mesh (ring_3D_tphi.vtu, built by ../files/formulation_meshes.py): the box |x|, |y| < 2 m, |z| < 1.5 m,
#   3072 HEX20 cells of 0.25 m. Ring TYPE / MATERIAL 2: 0.5 < max(|x|, |y|) < 1, |z| < 0.25; cut layer TYPE 3 (MATERIAL
#   3, air): max(|x|, |y|) < 0.5, 0 < z < 0.25; air TYPE 1.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order)
#   * 1: TPHI3D_edge / NON_CONDUCTING  / 1   air: PHI only
#   * 2: TPHI3D_edge / SOLID_CONDUCTOR / 1   ring: T and PHI
#   * 3: TPHI3D_edge / CUT_LAYER       / 1   the cut: T without conduction
#   Materials: 1 air; 2 rho = 1e-6 Ohm m; 3 air.
#   DOF couplings (set by the macro): 1, 2, 3 -> ["T"] (the crossing edges of the cut and the two faces of the layer).
#
# Boundary conditions: PHI = -H0 z on the box (B0 = mu0 H0 = 1 T along z from t = 0); the gauges above.
# Time history: backward Euler, dt = 0.01 s, 30 steps (the L / R time of the ring is about 0.07 s).
# Results: output/ring_3D_tphi_out_<n>.vtu (H, B, Jt, E, T, PHI), ring_3D_tphi_steps.csv, ring_3D_tphi.log;
#   printed: the ring current of each step.
#
# How to run:  cd examples/ring_3D_tphi && python macro.py
# Documentation: "Element library - TPHI3D_edge" and "Cuts", "Theory - scalar potentials".
# Variants (the tutorial "Cuts and fixed-T coils"): the current imposed (_T = I on the crossing edges, no coupling)
#   or the voltage (the coupling and VLTG = V on one crossing edge).
# =====================================================================================================================

import os

import numpy as np

from pytec3d import Solver, mesh, selection

MU0 = 4e-7 * np.pi
B0, DT, STEPS = 1.0, 0.01, 30

# 1. Solver, mesh and model data. Everything below is set by selections (as in MAPDL): select points or cells, then
#    set_value_on_selection() writes one value, or one value per selected id (in increasing id order), into an array.
pytec3d = Solver()
filename = "ring_3D_tphi"
pytec3d.read_model(os.path.join(os.getcwd(), filename + ".vtu"))
pytec3d.model.read_json("model_data.json")
grid = pytec3d.model.mesh
tol = 1e-6

# 2. Gauges of the ring. Its surface: the points of the ring cells that also belong to the other cells (nsle of the
#    ring, then reselect the points of the others); n x T = 0 on its edges: _T = 0 (component 2) on its midside points.
#    Inside: the points of the ring cells, minus those of the other cells; PHI = 0 at the corner points.
selection.select(grid, "TYPE", 2, location="elem")
selection.nsle(grid)
selection.select(grid, "TYPE", [1, 3], location="elem")
selection.nsle(grid, mode="r")                              # surface points of the ring
selection.select(grid, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(grid, "_T", 0.0, 2)
selection.select(grid, "TYPE", 2, location="elem")
selection.nsle(grid)
selection.select(grid, "TYPE", [1, 3], location="elem")
selection.nsle(grid, mode="u")                              # inner points of the ring
selection.select(grid, "MidsideNodes", 0, mode="r")
selection.set_value_on_selection(grid, "_PHI", 0.0)

# 3. The cut: the edges of the layer cells (TYPE 3). Crossing edges: one end on each face of the layer, here the
#    vertical ones (midside points at z = 0.125), all oriented along +z by the solver (from the lexicographically
#    smaller end point to the larger; tree_gauge.edges() gives the orientation of any edge). Planar edges (midside points
#    on z = 0 and z = 0.25): _T = 0. The crossing edges on the ring's inner surface belong to the cut: they leave the
#    n x T = 0 set (_T = NaN: free). All crossing edges are coupled on one DOF (CP_0 = 1), the ring current, unknown
#    here (an induced current); CP_0 = 2 and 3 on the corner points of the two faces give the solver the orientation of
#    each crossing edge (from the face 2, c = 0, to the face 3, c = 1).
selection.select(grid, "TYPE", 3, location="elem")
selection.nsle(grid)
selection.select(grid, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(grid, "_T", 0.0, 2)        # all the edges of the layer, then the crossing ones:
selection.select_pos(grid, 2, (0.125 - tol, 0.125 + tol), mode="r")
selection.set_value_on_selection(grid, "_T", np.nan, 2)
selection.set_value_on_selection(grid, "CP_0", 1)           # the coupled DOF: the ring current
for face, z in ((2, 0.0), (3, 0.25)):
    selection.select(grid, "TYPE", 3, location="elem")
    selection.nsle(grid)
    selection.select(grid, "MidsideNodes", 0, mode="r")
    selection.select_pos(grid, 2, (z - tol, z + tol), mode="r")
    selection.set_value_on_selection(grid, "CP_0", face)
pytec3d.model.add_dof_couplings({1: ["T"], 2: ["T"], 3: ["T"]})

# 4. The applied field: the total potential PHI = -H0 z on the box, from the first step on: one value per selected
#    point, computed from their z coordinates (in increasing id order, as selected_ids).
selection.select_all(grid)
selection.nsel_ext(grid)
z = mesh.points(grid)[selection.selected_ids(grid), 2]
selection.set_value_on_selection(grid, "_PHI", -B0 / MU0 * z)

# 5. Load the model; backward Euler.
pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 1.0

# The ring current: the flux of Jt_y through the section x > 0.5, 0 < y < 0.25 of the ring (cell values x 0.25 x 0.25:
# the cells of the ring with their centroid in that range).
selection.select(grid, "TYPE", 2, location="elem")
selection.select_pos(grid, 0, (0.0, 2.0), mode="r", location="elem")
selection.select_pos(grid, 1, (0.0, 0.25), mode="r", location="elem")
for step in range(1, STEPS + 1):
    pytec3d.time_step["time"] = step * DT
    pytec3d.solve()
    pytec3d.write_output()
    current = selection.get_selection_values(grid, "Jt", 1, location="elem").sum() * 0.25 * 0.25
    print(f"t = {step * DT:.2f} s: ring current {current / 1e3:8.2f} kA")

pytec3d.save_model()

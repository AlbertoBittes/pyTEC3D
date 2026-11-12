# =====================================================================================================================
# core_3D_phi - magnetostatics with the scalar potential (PHI3D_nodal): a coil on a closed, saturating iron frame
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A coil wound on one leg of a closed square iron frame, solved with PHI3D_nodal: H = Hs - grad(PHI), the scalar
# potential PHI (nodal) and the source field Hs (edge, NEDELEC_ONE) both unknowns of ONE static system.
#   * the source current: the cell array Js of the STRANDED_CONDUCTOR cells (any Js: curl Hs is its solenoidal part);
#   * the gauge of Hs: a tree of all the mesh edges grown from the box (tree_gauge(..., array_name="_Hs")) and
#     n x Hs = 0 on the box (_Hs = 0 on the box edges); needed in every PHI3D_nodal model;
#   * the iron is a closed loop around the current, the case where a reduced potential with a Biot-Savart source loses
#     its accuracy (cancellation) and MAPDL needs its generalized scalar potential: here Hs and grad(PHI) are in the
#     same edge space, they cancel exactly in the iron, and no cut is needed (PHI stays single valued: Hs carries the
#     current);
#   * a B-H curve on the iron (law bh_curve_h: input H, output mu), static load steps solved by Newton.
#
# Geometry and mesh (core_3D_phi.vtu, built by ../files/formulation_meshes.py): the box |x|, |y| < 2 m, |z| < 1.5 m,
#   3072 HEX20 cells of 0.25 m. Iron frame TYPE / MATERIAL 2: 0.5 < max(|x|, |y|) < 1, |z| < 0.5 (legs 0.5 x 1 m);
#   coil TYPE / MATERIAL 3: a 0.25 m thick winding around the leg 0.5 < x < 1, over |y| < 0.5; air TYPE 1.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order)
#   * 1: PHI3D_nodal / NON_CONDUCTING     / 1   air
#   * 2: PHI3D_nodal / SOFT_MAGNETIC      / 1   iron (mur replaced by the B-H curve)
#   * 3: PHI3D_nodal / STRANDED_CONDUCTOR / 1   the coil: source current Js
#   option2 = 1: PHI3D_nodal is first order only (the gradients of first-order PHI are the edge functions of Hs).
#   Materials: 1 air; 2 iron (mur = 1000, replaced by the curve set below); 3 coil (mur = 1).
#   Nonlinear laws: bh_curve_h on MATERIAL 2, added by this macro: B = mu0 H + Bs H / (H + Hk), Bs = 1.8 T, Hk = 300 A/m.
#
# Boundary conditions: the box is a flux wall (natural condition of PHI, n . B = 0), PHI = 0 on one node; Hs gauge.
# Loads: four static steps, NI = 0.5, 2, 10 and 50 kA-turns (J0 = NI / (0.25 m x 1 m) in the winding).
# Results: output/core_3D_phi_out_<n>.vtu (H, B, Jt = curl Hs, PHI, Hs), core_3D_phi_steps.csv, core_3D_phi.log;
#   printed: the mean B_y of the opposite leg (the flux of the frame) and the Newton iterations of each step.
#
# How to run:  cd examples/core_3D_phi && python macro.py
# Documentation: "Element library - PHI3D_nodal", "User manual - Nonlinear laws", "Theory - scalar potentials".
# =====================================================================================================================

import os

import numpy as np

from pytec3d import Solver, mesh, selection, tree_gauge

MU0 = 4e-7 * np.pi

# 1. Solver, mesh and model data. The cell centroids (cell array cell_centroids) are used by the position selections
#    of cells and to give Js.
pytec3d = Solver()
filename = "core_3D_phi"
pytec3d.read_model(os.path.join(os.getcwd(), filename + ".vtu"))
pytec3d.model.read_json("model_data.json")
grid = pytec3d.model.mesh
centroids = mesh.cell_centroids(grid)

# 2. The B-H curve of the iron, rows [H (A/m), B (T)], H and B increasing from 0; the curve goes on along its last
#    segment (slope about mu0 here: fully saturated).
H = np.concatenate(([0.0], np.geomspace(1.0, 1e8, 400)))
pytec3d.model.set_bh_curve(2, np.column_stack((H, MU0 * H + 1.8 * H / (H + 300.0))))

# 3. The gauge of the source field Hs: a spanning tree of the edges of all the cells, grown from the box, its edges
#    fixed to zero (_Hs = 0 on their midside nodes, component 2), and n x Hs = 0 on the box: _Hs = 0 on every edge
#    of the box (the edge DOFs live on the midside nodes; the value on the corner nodes is not used).
selection.select_all(grid)
selection.nsel_ext(grid)                                    # the points of the box
tree_gauge.tree_gauge(grid, array_name="_Hs")               # tree of the selected cells, rooted on the selected points
selection.set_value_on_selection(grid, "_Hs", 0.0, 2)       # n x Hs = 0 on the box

# 4. PHI: a free boundary is a flux wall; PHI is defined up to a constant: fix it on one node of the box.
box = selection.selected_ids(grid, "node")
selection.select(grid, "ID", int(box[0]))
selection.set_value_on_selection(grid, "_PHI", 0.0)

# 5. The source current density of the winding at 1 ampere-turn, per coil cell. The current circulates about +y
#    around the leg (centre x = 0.75, half width 0.25; half height 0.5): along -z on the side x > 1, +z on the side
#    x < 0.5, +x above the leg, -x below, J0 = 1 / (0.25 x 1.0) A/m^2. The corner cells take half of each direction:
#    that Js is not exactly divergence-free, and need not be: Hs only takes its solenoidal part. One vector per
#    selected cell (in increasing id order), all the components (component=None).
selection.select(grid, "TYPE", 3, location="elem")
c = centroids[selection.selected_ids(grid, "elem")]
u, w = c[:, 0] - 0.75, c[:, 2]
du, dw = np.abs(u) - 0.25, np.abs(w) - 0.5                  # how far a coil cell is beyond the leg in x and in z
along_z = np.where(np.isclose(du, dw), 0.5, (du > dw).astype(float))
j_unit = np.column_stack((np.sign(w) * (1.0 - along_z), np.zeros(len(c)), -np.sign(u) * along_z)) / (0.25 * 1.0)

# 6. Load the model; static steps (no time derivative: transient_effects False). The coil cells stay selected: at
#    each step Js = NI j_unit on them.
pytec3d.load_model_into_solver()
pytec3d.time_step["transient_effects"] = False

opposite = (mesh.array_values(grid, "TYPE", "elem") == 2) & (centroids[:, 0] < -0.5) & (np.abs(centroids[:, 1]) < 0.5)
volumes = mesh.cell_sizes(grid)
for step, ni in enumerate((500.0, 2e3, 1e4, 5e4), start=1):
    pytec3d.time_step["time"] = float(step)
    selection.set_value_on_selection(grid, "Js", ni * j_unit, None, location="elem")
    pytec3d.solve()
    pytec3d.write_output()
    b = mesh.array_values(grid, "B", "elem")
    flux_density = abs(volumes[opposite] @ b[opposite, 1]) / volumes[opposite].sum()
    print(f"step {step}: NI = {ni:8.0f} A-turns, mean |B_y| in the opposite leg = {flux_density:.3f} T, "
          f"Newton iterations {pytec3d.steps[-1]['newton_iterations']}")

pytec3d.save_model()

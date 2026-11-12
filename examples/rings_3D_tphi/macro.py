# =====================================================================================================================
# rings_3D_tphi - T-phi (TPHI3D_edge): several conductors with their own cuts, connected in SERIES and in PARALLEL
# =====================================================================================================================
#
# What this example shows
# -----------------------
# Two square conducting rings stacked along z, each with its own cut (a one-cell CUT_LAYER layer spanning its hole),
# driven by a voltage around the loop (the load VLTG on a cut DOF), no applied field. Two analyses:
#   * SERIES: the crossing edges of both cuts are coupled on ONE DOF (the same CP ids 1, 2, 3 for both layers): one
#     current through both rings, as two turns of a winding in series; VLTG on one crossing edge is the voltage of the
#     whole series loop;
#   * PARALLEL: each cut its own DOF (CP ids 1, 2, 3 for the lower ring, 4, 5, 6 for the upper one), both loaded with
#     the same VLTG: the rings share the voltage, each carries its own current.
# During the transient the rings are coupled by their mutual inductance (the current rises more slowly than in one
# ring alone); a last, very long step gives the DC state: series I = V / (R1 + R2), parallel I_k = V / R_k.
#   * the gauges of each ring (n x T = 0 on its surface, PHI = 0 inside), the cut of each ring (planar edges _T = 0,
#     crossing edges coupled, the corner nodes of the two faces of each layer labelled with the orientation ids);
#   * every value is set by selections and set_value_on_selection;
#   * rules (user manual, "Cuts"): one free cut DOF per hole, cuts do not cross; the coupling ids of a connection
#     are the circuit: a total current shared between parallel conductors or a resistor in series is not available.
#
# Geometry and mesh (rings_3D_tphi.vtu, built by ../files/formulation_meshes.py): the box |x|, |y| < 2 m, |z| < 1.5 m,
#   3072 HEX20 cells of 0.25 m. Rings 0.5 < max(|x|, |y|) < 1, |z -+ 0.5| < 0.25: lower TYPE / MATERIAL 2, upper 3;
#   their cut layers max(|x|, |y|) < 0.5, zc < z < zc + 0.25: lower TYPE 4, upper TYPE 5 (MATERIAL 4, 5: air); air 1.
#
# Model data (model_data.json)
#   Element types: 1 TPHI3D_edge / NON_CONDUCTING; 2, 3 TPHI3D_edge / SOLID_CONDUCTOR; 4, 5 TPHI3D_edge / CUT_LAYER;
#   all first order (option2 = 1). Materials: 2, 3 rho = 1e-6 Ohm m.
#   DOF couplings (set by the macro): series 1, 2, 3 -> ["T"]; parallel 1 ... 6 -> ["T"].
#
# Boundary conditions: PHI = 0 on the box; the gauges and cuts above.
# Loads and time history: VLTG = 0.1 V from t = 0; backward Euler, 10 steps of 0.02 s, then one step of 1000 s (DC).
# Results: series/ and parallel/ folders (output/, steps, log); printed: the current of each ring at each step and
#   the DC check.
#
# How to run:  cd examples/rings_3D_tphi && python macro.py
# Documentation: "Element library - Cuts (several conductors, series, parallel)", "DOF Couplings".
# =====================================================================================================================

import os
import shutil

import numpy as np

from pytec3d import Solver, mesh, selection, tree_gauge

V0, DT, STEPS, T_DC = 0.1, 0.02, 10, 1000.0
RINGS = ((2, 4, -0.5), (3, 5, 0.5))                  # (ring TYPE, its layer TYPE, the centre z of the ring)
tol = 1e-6


def run(connection):
    """One analysis in its own folder: ``connection`` "series" or "parallel"; returns the currents of the rings."""
    folder = os.path.join(os.getcwd(), connection)
    os.makedirs(folder, exist_ok=True)
    for name in ("rings_3D_tphi.vtu", "model_data.json"):
        shutil.copy(name, folder)
    here = os.getcwd()
    os.chdir(folder)
    pytec3d = Solver()
    pytec3d.read_model(os.path.join(folder, "rings_3D_tphi.vtu"))
    pytec3d.model.read_json("model_data.json")
    grid = pytec3d.model.mesh
    xyz = mesh.points(grid)
    others = [1, 2, 3, 4, 5]
    couplings = {}
    for k, (ring, layer, zc) in enumerate(RINGS):
        rest = [t for t in others if t != ring]
        # 1. Gauges of the ring: _T = 0 on the midside nodes of its surface (shared with other cells), _PHI = 0 at its
        #    inner corner nodes.
        selection.select(grid, "TYPE", ring, location="elem"); selection.nsle(grid)
        selection.select(grid, "TYPE", rest, location="elem"); selection.nsle(grid, mode="r")
        selection.select(grid, "MidsideNodes", 1, mode="r")
        selection.set_value_on_selection(grid, "_T", 0.0, 2)
        selection.select(grid, "TYPE", ring, location="elem"); selection.nsle(grid)
        selection.select(grid, "TYPE", rest, location="elem"); selection.nsle(grid, mode="u")
        selection.select(grid, "MidsideNodes", 0, mode="r")
        selection.set_value_on_selection(grid, "_PHI", 0.0)
        # 2. Its cut: the edges of its layer; planar edges _T = 0; crossing edges (vertical, midside at zc + 0.125,
        #    all oriented along +z) free and coupled. Series: the same ids for both cuts (one DOF); parallel: ids of
        #    their own (base 0 for the lower ring, 3 for the upper).
        base = 0 if connection == "series" else 3 * k
        selection.select(grid, "TYPE", layer, location="elem"); selection.nsle(grid)
        selection.select(grid, "MidsideNodes", 1, mode="r")
        selection.set_value_on_selection(grid, "_T", 0.0, 2)
        selection.select_pos(grid, 2, (zc + 0.125 - tol, zc + 0.125 + tol), mode="r")
        selection.set_value_on_selection(grid, "_T", np.nan, 2)
        selection.set_value_on_selection(grid, "CP_0", base + 1)
        crossing = selection.selected_ids(grid)
        for face, z in ((base + 2, zc), (base + 3, zc + 0.25)):      # the corner nodes of the two faces of the layer
            selection.select(grid, "TYPE", layer, location="elem"); selection.nsle(grid)
            selection.select(grid, "MidsideNodes", 0, mode="r")
            selection.select_pos(grid, 2, (z - tol, z + tol), mode="r")
            selection.set_value_on_selection(grid, "CP_0", face)
        couplings.update({base + 1: ["T"], base + 2: ["T"], base + 3: ["T"]})
        # 3. The voltage: VLTG on one crossing midside node of the cut DOF (series: once, for the whole loop).
        if connection == "parallel" or k == 0:
            selection.select(grid, "ID", int(crossing[0]))
            selection.set_value_on_selection(grid, "VLTG", V0)
    pytec3d.model.add_dof_couplings(couplings)
    # 4. The box: PHI = 0 (no applied field).
    selection.select_all(grid); selection.nsel_ext(grid)
    selection.set_value_on_selection(grid, "_PHI", 0.0)

    # 5. Load; backward Euler; the current of each ring from Jt_y across its section x > 0.5, 0 < y < 0.25.
    pytec3d.load_model_into_solver()
    pytec3d.time_step["theta"] = 1.0
    c = mesh.cell_centroids(grid)
    types = mesh.array_values(grid, "TYPE", "elem")
    sections = [(types == ring) & (c[:, 0] > 0) & (c[:, 1] > 0) & (c[:, 1] < 0.25) for ring, _, _ in RINGS]
    currents = []
    for t in list(DT * np.arange(1, STEPS + 1)) + [T_DC]:
        pytec3d.time_step["time"] = float(t)
        pytec3d.solve()
        pytec3d.write_output()
        jy = mesh.array_values(grid, "Jt", "elem")[:, 1]
        currents.append([jy[sec].sum() * 0.25 * 0.25 for sec in sections])
        print(f"{connection:8s} t = {t:7.2f} s: I lower = {currents[-1][0] / 1e3:8.3f} kA, "
              f"I upper = {currents[-1][1] / 1e3:8.3f} kA")
    pytec3d.close()
    os.chdir(here)
    return np.array(currents)


series = run("series")
parallel = run("parallel")
r1, r2 = V0 / parallel[-1, 0], V0 / parallel[-1, 1]                  # each ring's resistance, from the parallel DC
print(f"DC: R_lower = {r1 * 1e6:.4f} uOhm, R_upper = {r2 * 1e6:.4f} uOhm; series current "
      f"{series[-1, 0] / 1e3:.4f} kA = V / (R1 + R2) = {V0 / (r1 + r2) / 1e3:.4f} kA; "
      f"the same in both rings at every step: {np.allclose(series[:, 0], series[:, 1], rtol=1e-12)}")

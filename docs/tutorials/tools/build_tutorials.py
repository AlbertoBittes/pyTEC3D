"""Build the tutorial notebooks from the test cases (besides tutorial_coil_2D), with their outputs.

Usage, from docs/tutorials/tools, with the pytec3d environment (pyvista installed):
    python build_tutorials.py ..                  # all of them
    python build_tutorials.py .. edge_elements    # one of them
"""
import os
import sys

import nbbuild

SETUP = '''
import os, shutil, tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

import pytec3d
from pytec3d import Solver, selection
from pytec3d import mesh as meshes

pv.set_jupyter_backend("static")          # images in the notebook; "trame" gives interactive views
pytec3d.set_log_level("ERROR")            # warnings ("No loads" at every voltage-driven step...) stay in the job log

CASE = Path("../../tests/{case}").resolve()
WORK = Path(tempfile.mkdtemp(prefix="pytec3d_{case}_"))
for name in {files}:
    shutil.copy(CASE / name, WORK / name)
os.chdir(WORK)
print("working in a temporary folder:", WORK.name)

import contextlib, io
def quietly():
    """Keep the report tables that the solver prints (run provenance, run summary) out of the notebook: they are in
    the job log <model name>.log of the working folder."""
    return contextlib.redirect_stdout(io.StringIO())
'''


# ======================================================================================================================
def stranded_coil_iron():
    nb = nbbuild.Notebook()
    nb.md("""
# Stranded coil, circuit equation and saturating iron (2D axisymmetric)

This notebook follows `tests/coil_2D_avoltemf_ferr`: a 1000-turn **stranded coil** fed by a voltage source, a passive
copper ring above it, and a rectangular iron frame around both, described by a **B-H curve**. It shows:

* the circuit unknowns of a stranded coil, `VOLT` (terminal voltage) and `EMF` (induced voltage), and how to read the
  coil current from them, `i = (VOLT - EMF) / R`;
* a voltage-driven transient: the step response of an inductive circuit, `i(t) -> V / R` with the time constant
  `L / R`;
* the effect of the iron: a larger inductance, then saturation, which a nonlinear law handles with Newton's method;
* the report of the operating range of the B-H law at every step.

We run the same model twice, without and with the iron, and compare. Requirements: the `pytec3d` environment with
`pip install -e ".[notebooks]"` (matplotlib; pyvista is part of the conda environment). See the user manual (*Element Library: stranded coil*,
*Nonlinear Materials*) and the theory (*Stranded coils*, *Newton-Raphson*).
""")
    nb.code(SETUP.format(case="coil_2D_avoltemf_ferr", files='("coil_2D_avoltemf_ferr.vtu", "model_data.json")'))
    nb.md("""
## The model

Axisymmetric about the y axis (x = r): the coil (TYPE 23, `STRANDED_COIL`) at r = 0.5 m, z = -0.5 m, the passive ring
(TYPE 24, `SOLID_CONDUCTOR`: in 2D a massive conductor without terminal, a closed ring) at z = +0.5 m, and the iron
frame (TYPE 26). The coil real
constants are in CONSTANT 1:
""")
    nb.code('''
import json
data = json.load(open("model_data.json"))
coil = data["constants"]["1"]
print(json.dumps(coil))
R, N = coil["resistance"], coil["turns"]
print(f"R = {R} Ohm, N = {N:.0f} turns: a 1500 V source drives {1500 / R:.0f} A, i.e. {1500 / R * N / 1e3:.0f} kA-turns")
''')
    nb.md("""
## One function for the whole analysis

`run(with_iron)` builds the analysis as the test macro does, then applies a 1500 V step at t = 0 and follows it for
1 s with backward Euler (100 steps of 10 ms):

1. read the mesh and the model data;
2. the iron frame: either TYPE / MATERIAL 26 with the B-H curve, or air (TYPE / MATERIAL 21);
3. A_z = 0 on the external boundary;
4. the coil terminal: `CP_0 = 1` on the coil nodes, coupling `["EMF", "VOLT"]` (model data `dof_couplings["1"]`);
5. the passive ring: nothing to set (a 2D solid conductor has no electric unknown);
6. load the model, then the steps with `_VOLT = 1500` on the coil master node.

After each step it records the terminal quantities at the master node, the coil current, the Newton iterations and
the largest |B| in the iron (the *law range* of the step report). The report tables of each run go to its log file
(`air.log`, `iron.log`), not to the notebook.
""")
    nb.code('''
BH = np.array([[0.0, 0.0], [100.0, 0.372], [200.0, 0.744], [255.0, 0.936], [299.1, 1.102], [373.2, 1.287],
               [447.5, 1.415], [521.3, 1.488], [595.2, 1.553], [743.3, 1.641], [1113.8, 1.763], [1484.7, 1.838],
               [1855.2, 1.879], [2226.9, 1.922], [3661.0, 1.998], [5870.7, 2.092], [10289.3, 2.220],
               [20396.8, 2.435]])                                       # rows [H (A/m), B (T)]

def run(with_iron, volt=1500.0, dt=1e-2, steps=100):
    with contextlib.redirect_stdout(io.StringIO()):     # the report tables go to the log file only
        return _run(with_iron, volt, dt, steps)   # (quietly() of the setup cell does the same)

def _run(with_iron, volt, dt, steps):
    solver = Solver(log_file="iron.log" if with_iron else "air.log")
    solver.read_model("coil_2D_avoltemf_ferr.vtu")
    solver.model.read_json("model_data.json")
    mesh = solver.model.mesh

    # 2. iron frame (the mesh already has TYPE 26 on it) or air
    selection.select(mesh, "TYPE", 26, location="elem")
    if with_iron:
        solver.model.set_bh_curve(26, BH)
    else:
        selection.set_value_on_selection(mesh, "TYPE", 21, location="elem")
        selection.set_value_on_selection(mesh, "MATERIAL", 21, location="elem")

    # 3. A_z = 0 on the external boundary
    selection.select_all(mesh); selection.nsel_ext(mesh)
    selection.set_value_on_selection(mesh, "_A", 0.0, 2)

    # 4. coil terminal (5. the passive ring is a solid conductor: nothing to set)
    selection.select(mesh, "MATERIAL", 23, location="elem"); selection.nsle(mesh)
    selection.set_value_on_selection(mesh, "CP_0", 1)
    master = int(selection.selected_ids(mesh, "node")[0])

    # 6. load, then the voltage step
    solver.load_model_into_solver()
    solver.time_step["theta"] = 1.0
    selection.select(mesh, "ID", master)
    selection.set_value_on_selection(mesh, "_VOLT", volt)
    rows = []
    for n in range(1, steps + 1):
        solver.time_step["time"] = n * dt
        solver.solve()
        v = float(meshes.array_values(mesh, "VOLT", "node")[master])
        e = float(meshes.array_values(mesh, "EMF", "node")[master])
        report = solver.steps[-1]["report"]
        b_iron = max([r["max_input"] for r in report["law_ranges"]], default=0.0)
        rows.append((n * dt, v, e, (v - e) / R, solver.steps[-1]["newton_iterations"], b_iron))
    if with_iron:
        solver.write_output()          # keep the last state for the plots below
    grid = pv.wrap(mesh).copy()
    solver.close()
    return np.array(rows), grid

air, _ = run(with_iron=False)
iron, grid = run(with_iron=True)
print("final coil current: air %.2f A, iron %.2f A (V / R = %.0f A)" % (air[-1, 3], iron[-1, 3], 1500 / R))
print("Newton iterations per step (iron): min %d, max %d" % (iron[:, 4].min(), iron[:, 4].max()))
''')
    nb.md("""
## Step response

In air the current rises with a single time constant `L / R` (about 20 ms here). The iron frame guides the flux and
multiplies the inductance: the rise is an order of magnitude slower. The EMF is the voltage induced in the winding by
the change of the linked flux, `EMF = dPsi/dt`; it starts close to the source voltage and decays as the current
settles, `i = (VOLT - EMF) / R`. In the last part of the transient the iron reaches the knee of its B-H curve: its
differential permeability, and the inductance seen by the source, decrease, and Newton needs more than one iteration
per step (in the linear first segment of the curve one iteration is exact).
""")
    nb.code('''
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
for data, label in ((air, "air"), (iron, "iron")):
    ax[0].plot(data[:, 0] * 1e3, data[:, 3], label=label)
    ax[1].plot(data[:, 0] * 1e3, data[:, 2], label=label)
ax[0].axhline(1500 / R, color="k", lw=0.5, ls="--"); ax[0].set_ylabel("coil current i [A]")
ax[1].set_ylabel("EMF [V]")
ax[2].plot(iron[:, 0] * 1e3, iron[:, 5], "C1"); ax[2].set_ylabel("largest |B| in the iron [T]")
for a in ax:
    a.set_xlabel("t [ms]"); a.grid(alpha=0.3)
ax[0].legend(); plt.tight_layout(); plt.show()
''')
    nb.md("""
## Where the iron works on its B-H curve

The law `bh_curve` uses the table as a piecewise-linear curve, extended along its last segment. The largest |B| of
each step (the element mean, as reported by the solver) shows how far into saturation the frame goes.
""")
    nb.code('''
fig, ax = plt.subplots(figsize=(5.5, 3.6))
ax.plot(BH[:, 0], BH[:, 1], "k.-", label="B-H table")
h_op = np.interp(iron[:, 5], BH[:, 1], BH[:, 0])
ax.plot(h_op, iron[:, 5], "o", color="C1", ms=4, label="largest |B| per step")
ax.set_xscale("symlog", linthresh=100); ax.set_xlabel("H [A/m]"); ax.set_ylabel("B [T]"); ax.grid(alpha=0.3)
ax.legend(); plt.show()
''')
    nb.md("""
## The field at the end of the transient

The flux is guided by the frame; the colour shows |B| on the nodes, the arrows the in-plane B.
""")
    nb.code('''
grid.point_data["|B|"] = np.linalg.norm(grid.point_data["B"], axis=1)
surface = grid.extract_surface(algorithm=None)
outline = surface.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False)
regions = grid.threshold(value=(22, 30), scalars="TYPE").extract_surface(algorithm=None).extract_feature_edges(
    boundary_edges=True, feature_edges=False, manifold_edges=False)
arrows = surface.glyph(orient="B", scale=False, factor=0.035, tolerance=0.02)

pl = pv.Plotter(window_size=(700, 900))
pl.add_mesh(surface, scalars="|B|", cmap="plasma", scalar_bar_args={"title": "|B| [T]"})
pl.add_mesh(grid.threshold(value=(23, 24), scalars="TYPE").extract_surface(algorithm=None), color="#c87533")
pl.add_mesh(regions, color="white", line_width=2)
pl.add_mesh(arrows, color="black")
pl.view_xy(); pl.camera.zoom(1.3)
pl.show()
''')
    nb.md("""
## Try it yourself

* Halve the source voltage: the iron saturates less and the two curves get closer.
* Drive the coil by its current instead (`CURR` on the master node, `_VOLT` free): the EMF is then the only output of
  the circuit, and `VOLT = R i + EMF` gives the voltage the source must supply.
* Use the time-integrated element type (`"element_type": "AV2D_nodal_timint"` in the model data): `VOLT` and `EMF`
  become the time integrals of the voltages; `_VOLT` must then be a ramp `V t`.
""")
    return nb


# ======================================================================================================================
def coil_3d_sector():
    nb = nbbuild.Notebook()
    nb.md("""
# A 3D coil on a cyclic sector

This notebook follows `tests/coil_3D_avoltemf`: the coil and the passive ring of the 2D examples, now in 3D, modelled
on a **30 degree sector** about the y axis with cyclic symmetry. It shows:

* selections in a cylindrical coordinate system, and far-field / axis conditions on a nodal vector potential;
* an **element coordinate system** in which the winding direction of the coil is simply `T = (0, 1, 0)` (toroidal);
* **cyclic symmetry**: `periodic_boundaries()` pairs the nodes of the two sector faces, and the nodal frames of one
  face are rotated (`node_angles`) so that the same cylindrical components are coupled;
* the `symmetry` real constant: the loads refer to the whole coil although only 1/12 is modelled;
* the eddy currents of the passive ring after a current step, and 3D views of the full ring with pyvista.

See the user manual (*DOF Couplings and Periodic Boundaries*, *Coordinate Systems*, *Boundary Conditions*).
""")
    nb.code(SETUP.format(case="coil_3D_avoltemf", files='("coil_3D_avoltemf.vtu", "model_data.json")'))
    nb.md("""
## Model and boundary conditions

The coil (TYPE 23, `STRANDED_COIL`, 2000 turns) at r = 0.5 m, y = -0.5 m; the passive ring (TYPE 24,
`SOLID_CONDUCTOR`, with the electric potential V) at y = +0.5 m. In the Y-axis cylindrical system,
coordinate 0 is r, 1 is theta (degrees, from +x towards -z), 2 is the axial coordinate y.
""")
    nb.code('''
from pytec3d import YAxisCylindricalCoordinateSystem
cyl = YAxisCylindricalCoordinateSystem()
tol = 1e-5

solver = Solver()
solver.read_model("coil_3D_avoltemf.vtu")
solver.model.read_json("model_data.json")
mesh = solver.model.mesh
print(solver.model, "- coil constants:", solver.model.data["constants"]["1"])

# far field: A = 0 on y = +-1 and r = 1
selection.select_pos(mesh, 2, (1 - tol, 1 + tol), coordinate_system=cyl)
selection.select_pos(mesh, 2, (-1 - tol, -1 + tol), mode="a", coordinate_system=cyl)
selection.select_pos(mesh, 0, (1 - tol, 1 + tol), mode="a", coordinate_system=cyl)
for c in (0, 1, 2):
    selection.set_value_on_selection(mesh, "_A", 0.0, c)
# axis r = 0: the components normal to the axis vanish
selection.select_pos(mesh, 0, (-tol, tol), coordinate_system=cyl)
selection.set_value_on_selection(mesh, "_A", 0.0, 0)
selection.set_value_on_selection(mesh, "_A", 0.0, 2)
''')
    nb.md("""
## Coil terminal, element coordinate system, cyclic symmetry

The coil cells use the cylindrical system 1, so `T = (0, 1, 0)` of the real constants is the toroidal direction in
every cell. Its nodes share one `VOLT` and one `EMF` (CP id 1). The nodes of the face theta = 0 are paired with their
images at theta = 30 degrees in a fresh array `CP_1`, coupling A and V.
""")
    nb.code('''
import vtk
solver.model.set_coordinate_system(1, cyl)
selection.select(mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")
selection.nsle(mesh)
selection.set_value_on_selection(mesh, "CP_0", 1)
master = int(selection.selected_ids(mesh, "node")[0])

meshes.remove_array(mesh, "CP_1")
meshes.add_array(mesh, "CP_1", dtype=vtk.VTK_INT)
selection.select_pos(mesh, 1, (-tol, tol), coordinate_system=cyl)
pairs = solver.model.periodic_boundaries(cyl, ("A", "V"), dy=30.0, array_name="CP_1")
selection.select_pos(mesh, 1, (30 - tol, 30 + tol), coordinate_system=cyl)
selection.set_value_on_selection(mesh, "node_angles", 30.0, 2)
print(len(pairs), "periodic pairs")

# the potential of the passive ring floats: fix it on one of its nodes
selection.select(mesh, "MATERIAL", 24, location="elem"); selection.nsle(mesh)
selection.select(mesh, "ID", int(selection.selected_ids(mesh, "node")[0]))
selection.set_value_on_selection(mesh, "_V", 0.0)
''')
    nb.md("""
## A current step

The coil current jumps from 0 to 100 A at t = 0 (`CURR` on the master node; it is the current of the whole coil) and
is kept constant; backward Euler, 10 steps of 5 ms. The passive ring reacts with an eddy current that opposes the
change of flux and then decays with the L/R time of the ring. We estimate the ring current from the element current
density: the toroidal component of `Jt`, integrated over the sector and divided by the sector length (its mean over
the 30 degrees).
""")
    nb.code('''
with quietly():
    solver.load_model_into_solver()
solver.time_step["theta"] = 1.0
selection.select(mesh, "ID", master)
selection.set_value_on_selection(mesh, "CURR", 100.0)

centroids = meshes.array_values(mesh, "cell_centroids", "elem")
volume = meshes.array_values(mesh, "cell_volume", "elem")
ring_cells = meshes.array_values(mesh, "MATERIAL", "elem") == 24
local = cyl.to_local(centroids[ring_cells])
rows = []
for n in range(1, 11):
    solver.time_step["time"] = n * 5e-3
    solver.solve()
    jt = meshes.array_values(mesh, "Jt", "elem")[ring_cells]
    j_theta = cyl.vector_components_to_local(jt, centroids[ring_cells])[:, 1]
    i_ring = np.sum(j_theta * volume[ring_cells] / (local[:, 0] * np.deg2rad(30.0)))
    volt = meshes.array_values(mesh, "VOLT", "node")[master]
    emf = meshes.array_values(mesh, "EMF", "node")[master]
    rows.append((n * 5e-3, volt, emf, i_ring))
solver.write_output()
rows = np.array(rows)
print("  t [ms]   VOLT [V]   EMF [V]   I ring [A]")
for t, v, e, i in rows:
    print(f"  {t * 1e3:5.0f}   {v:8.2f}   {e:8.2f}   {i:9.2f}")
''')
    nb.md("""
With the current imposed, the terminal voltage is `VOLT = R i + EMF`. At the first step the current jumps from 0 to
100 A in 5 ms: the EMF `L di/dt` (about 80 kV, L of the order of 4 H) dominates. Afterwards the current is constant and
`VOLT` is `R i = 2000 V` plus the small voltage induced by the decaying eddy current of the ring.
""")
    nb.code('''
fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
ax[0].plot(rows[:, 0] * 1e3, rows[:, 3], "o-"); ax[0].set_ylabel("passive ring current [A]")
ax[1].plot(rows[:, 0] * 1e3, rows[:, 1], "o-", label="VOLT"); ax[1].plot(rows[:, 0] * 1e3, rows[:, 2], "s-", label="EMF")
ax[1].set_ylabel("[V]"); ax[1].legend()
for a in ax:
    a.set_xlabel("t [ms]"); a.grid(alpha=0.3)
plt.tight_layout(); plt.show()
''')
    nb.md("""
## The full ring in 3D

pyvista reads the VTK grid directly. The sector is copied 12 times, rotated about the y axis, to show the whole
machine: the coil (copper colour), the passive ring coloured by its toroidal eddy current density (opposite to the
coil current), and |B| on the face theta = 0.
""")
    nb.code('''
grid = pv.wrap(mesh).copy()
jt = grid.cell_data["Jt"]
grid.cell_data["J_theta"] = cyl.vector_components_to_local(jt, grid.cell_centers().points)[:, 1]
grid.point_data["|B|"] = np.linalg.norm(grid.point_data["B"], axis=1)
coil = grid.threshold(value=(23, 23), scalars="TYPE")
ring = grid.threshold(value=(24, 24), scalars="TYPE")
face = grid.slice(normal=(0, 0, 1), origin=(0, 0, -1e-6))          # the plane z = 0 (theta = 0)
limit = np.abs(ring["J_theta"]).max()

pl = pv.Plotter(window_size=(950, 700))
for k in range(12):
    pl.add_mesh(coil.rotate_y(30.0 * k, inplace=False), color="#c87533")
    pl.add_mesh(ring.rotate_y(30.0 * k, inplace=False), scalars="J_theta", cmap="coolwarm", clim=[-limit, limit],
                show_scalar_bar=(k == 0), scalar_bar_args={"title": "J_theta ring", "vertical": True, "position_x": 0.03, "position_y": 0.25, "height": 0.5, "width": 0.06})
pl.add_mesh(face, scalars="|B|", cmap="plasma", opacity=0.6, scalar_bar_args={"title": "|B| [T]", "vertical": True, "position_x": 0.88, "position_y": 0.25, "height": 0.5, "width": 0.06})
pl.camera_position = [(2.4, 1.6, 2.6), (0, 0, 0), (0, 1, 0)]
pl.show()
with quietly():
    solver.close()
''')
    nb.md("""
## Try it yourself

* Drive the coil by a voltage (`_VOLT` on the master node, `CURR = 0`) and plot the coil current
  `(VOLT - EMF) / R`.
* Remove the frame rotation (`node_angles`) and compare: the periodic coupling then ties the global x and z components
  of the two faces, which is wrong for a rotated sector.
* Run the edge-element version (`tests/coil_3D_avoltemf_edge`): no frame rotation is needed, a tree-cotree gauge is.
""")
    return nb


# ======================================================================================================================
def edge_elements():
    nb = nbbuild.Notebook()
    nb.md("""
# Edge elements and the tree-cotree gauge

This notebook follows `tests/two_wires_3D_edge_el_2order`: two rods along z crossing a box, a source rod with an
imposed current density and a **solid** rod whose ends are shorted (V = 0 on both end faces), discretized with
**Nedelec edge elements** (`AV3D_edge`). It shows:

* why edge elements need a quadratic mesh (their DOFs live on the midside nodes) and how their Dirichlet values are
  written (component 2 of `_A` on the midside nodes: the tangential A along each edge);
* the tree-cotree gauge: what `tree_gauge()` builds from the current selection, and how to look at the tree;
* a static step followed by transient steps, and the eddy current induced in the shorted rod.

See the user manual (*Element Library: AV3D_edge*, *Boundary Conditions: gauging*).
""")
    nb.code(SETUP.format(case="two_wires_3D_edge_el_2order", files='("two_wires_3D_edge_el_2order.vtu", "model_data.json")'))
    nb.code('''
from pytec3d import tree_gauge
solver = Solver()
solver.read_model("two_wires_3D_edge_el_2order.vtu")
solver.model.read_json("model_data.json")
mesh = solver.model.mesh
midside = meshes.array_values(mesh, "MidsideNodes", "node") == 1
print(solver.model, f"- {midside.sum()} midside nodes (one edge DOF each)")
for k, row in solver.model.data["element_types"].items():
    print("TYPE", k, row)
''')
    nb.md("""
## Tangential A = 0 on the six faces, then the gauge

On each face of the box we keep the midside nodes of the selection and fix component 2 of `_A`: this fixes the
circulation of A along every edge of the face, i.e. n x A = 0. Then the gauge: `tree_gauge()` builds a spanning tree of
the edges of the **selected cells** (all of them), grown from the **selected nodes** (the external boundary, where A
is already fixed), and fixes A = 0 on the tree edges. The ends of the solid rod are shorted with V = 0.
""")
    nb.code('''
for axis in (0, 1, 2):
    selection.select_pos(mesh, axis, (0.499, 0.50001))
    selection.select_pos(mesh, axis, (-0.501, -0.49999), mode="a")
    selection.select(mesh, "MidsideNodes", 1, mode="r")
    selection.set_value_on_selection(mesh, "_A", 0.0, 2)

for face in ("LOWER_FACE_ROD", "UPPER_FACE_ROD"):
    selection.select(mesh, face, 1)
    selection.set_value_on_selection(mesh, "_V", 0.0)

selection.select_all(mesh)
selection.nsel_ext(mesh)
info = tree_gauge.tree_gauge(mesh)
print({k: info[k] for k in ("num_edges", "num_tree", "num_cotree", "num_roots", "tree_identity_ok")})
''')
    nb.md("""
The tree is stored in the point array `TREE` (1 on the midside node of each tree edge). In a connected graph a
spanning tree has (number of nodes - number of roots) edges: the remaining (cotree) edges carry the solution.
""")
    nb.code('''
grid = pv.wrap(mesh).copy()
tree_points = pv.PolyData(grid.points[grid.point_data["TREE"] == 1])
rods = grid.threshold(value=(22, 30), scalars="TYPE")
pl = pv.Plotter(window_size=(800, 650))
pl.add_mesh(grid.outline(), color="black")
pl.add_mesh(tree_points, color="#2a7ab0", point_size=4, render_points_as_spheres=True, label="tree edges (midside nodes)")
pl.add_mesh(rods, color="#c87533", label="rods")
pl.add_legend(bcolor="white", size=(0.35, 0.12))
pl.camera_position = "iso"
pl.show()
''')
    nb.md("""
## Steps: static, rise, switch-off

* step 1, static: Js = 1 A/m^2 in the source rod (no eddy currents on a static step);
* step 2, transient (dt = 1 ms): Js = 100 A/m^2, a sudden rise: the shorted rod reacts;
* step 3, transient: Js = 0, the source is switched off.

We use theta = 1 (backward Euler) for these steps: they are sudden changes of the source, and theta < 1 would carry
part of the rate before each change into the next step (see the remark at the end of section 8 of the first tutorial).

The current through the shorted rod is estimated from its element current density, `J_z * volume / length` summed
over its cells (the mean over its length).
""")
    nb.code('''
with quietly():
    solver.load_model_into_solver()
cells = meshes.array_values(mesh, "MATERIAL", "elem")
volume = meshes.array_values(mesh, "cell_volume", "elem")
rod, source = cells == 24, cells == 23
length = 1.0                                            # the rods cross the box along z

def set_source(js):
    selection.select(mesh, "MATERIAL", 23, location="elem")
    selection.set_value_on_selection(mesh, "Js", js, 2, location="elem")

rows = []
for time, theta, transient, js in ((1e-3, 1.0, False, 1.0), (2e-3, 1.0, True, 100.0), (3e-3, 1.0, True, 0.0)):
    set_source(js)
    solver.time_step.update(time=time, theta=theta, transient_effects=transient)
    solver.solve()
    jz = meshes.array_values(mesh, "Jt", "elem")[:, 2]
    rows.append((time, js * np.sum(volume[source]) / length, np.sum(jz[rod] * volume[rod]) / length))
    print(f"t = {time:.0e} s ({'transient' if transient else 'static'}): source {rows[-1][1]:.3e} A, "
          f"shorted rod {rows[-1][2]: .3e} A, linear its {solver.steps[-1]['linear_iterations']}")
solver.write_output()
''')
    nb.md("""
The flux density on the mid-plane z = 0, and the current density in the rods (last step: the source is off, the eddy
current of the shorted rod decays):
""")
    nb.code('''
grid = pv.wrap(mesh).copy()
grid.point_data["|B|"] = np.linalg.norm(grid.point_data["B"], axis=1)
plane = grid.slice(normal="z", origin=(0, 0, 0))
rods = grid.threshold(value=(22, 30), scalars="TYPE")
# Jt: the imposed density Js in the source rod, the computed (eddy) current density in the shorted rod.
rods.cell_data["J_z"] = rods.cell_data["Jt"][:, 2]
arrows = plane.glyph(orient="B", scale="|B|", factor=0.1 / max(plane["|B|"].max(), 1e-30), tolerance=0.03)

pl = pv.Plotter(window_size=(850, 650))
pl.add_mesh(plane, scalars="|B|", cmap="plasma", opacity=0.9, scalar_bar_args={"title": "|B| [T]", "vertical": True, "position_x": 0.88, "position_y": 0.25, "height": 0.5, "width": 0.06})
pl.add_mesh(arrows, color="white")
pl.add_mesh(rods, scalars="J_z", cmap="coolwarm", scalar_bar_args={"title": "J_z [A/m^2]", "vertical": True, "position_x": 0.03, "position_y": 0.25, "height": 0.5, "width": 0.06})
pl.add_mesh(grid.outline(), color="black")
pl.camera_position = [(1.3, -1.6, 1.4), (0, 0, 0), (0, 0, 1)]
pl.show()
with quietly():
    solver.close()
''')
    nb.md("""
## Try it yourself

* Call `tree_gauge()` with no cell selected (`selection.unselect_all(mesh)` first): the tree is empty (`num_tree = 0`)
  and the curl-curl matrix of the air is singular. Add `"-mat_mumps_icntl_24": "1"` to
  `solver.solver_settings["petsc_options"]` to see the null pivots in the step report.
* Leave the rod ends free (no `_V`): the rod is open, and no net current flows through it, only local eddy currents.
""")
    return nb


# ======================================================================================================================
def loops_3d():
    nb = nbbuild.Notebook()
    nb.md("""
# Symmetry planes: a quarter model of two loops

This notebook follows `tests/3D_loop_cond_base`: a square source loop (imposed current density) below an identical
passive solid loop, modelled on one quarter of the space thanks to the symmetry planes x = 0 and y = 0. It shows:

* boundary conditions selected through named point arrays of the mesh (faces exported by the mesher);
* symmetry conditions on a nodal vector potential: on a plane crossed normally by the current, the tangential
  components of A vanish;
* a vector source current density given component by component on the sides and the corner of the loop;
* the current induced in the passive loop when the source is switched off, and a mirrored view of the whole system.

See the user manual (*Boundary Conditions*, *Loads*) and the examples page *Wires, Loops and Rings*.
""")
    nb.code(SETUP.format(case="3D_loop_cond_base", files='("3D_loop_cond_base.vtu", "model_data.json")'))
    nb.code('''
solver = Solver()
solver.read_model("3D_loop_cond_base.vtu")
solver.model.read_json("model_data.json")
mesh = solver.model.mesh
names = [mesh.GetPointData().GetArrayName(i) for i in range(mesh.GetPointData().GetNumberOfArrays())]
print("named point arrays:", [n for n in names if n.endswith(("_1", "_2"))])
''')
    nb.md("""
## Boundary conditions

* far field (z = +-0.125 m, x = 0.25 m, y = 0.25 m): A = 0;
* symmetry plane x = 0: the loop crosses it along x, so A (parallel to the current) is normal to the plane: A_y = A_z = 0;
* symmetry plane y = 0: A_x = A_z = 0;
* the passive loop is cut by the two planes: its potential V is fixed to 0 there (it is odd with respect to them).

As in the test macro, the conditions are written after the model is loaded: array values may change at any time.
""")
    nb.code('''
with quietly():
    solver.load_model_into_solver()

selection.select(mesh, "Z_BOUNDS_1", 1)
for face in ("Z_BOUNDS_2", "Y_BOUNDS_2", "X_BOUNDS_2"):
    selection.select(mesh, face, 1, mode="a")
for c in (0, 1, 2):
    selection.set_value_on_selection(mesh, "_A", 0.0, c)
selection.select(mesh, "X_BOUNDS_1", 1)
selection.set_value_on_selection(mesh, "_A", 0.0, 1); selection.set_value_on_selection(mesh, "_A", 0.0, 2)
selection.select(mesh, "Y_BOUNDS_1", 1)
selection.set_value_on_selection(mesh, "_A", 0.0, 0); selection.set_value_on_selection(mesh, "_A", 0.0, 2)
selection.select(mesh, "VOLT_0_1", 1); selection.select(mesh, "VOLT_0_2", 1, mode="a")
selection.set_value_on_selection(mesh, "_V", 0.0)
''')
    nb.md("""
## The source and the time history

The source loop carries J = 1e6 A/m^2 on cells of 1 cm x 1 cm, i.e. 100 A: +y on side 1, -x on side 2, and at 45
degrees on the corner cell. Step 1 is static (the initial state, no eddy currents); at t = 0+ the source is switched
off and the passive loop, a good conductor, takes over part of the flux: we follow it for 8 steps of 0.1 ms.
""")
    nb.code('''
b = 2 ** 0.5 / 2
def set_source(j):
    for group, components in (("COIL_NOCOND_1", {1: j}), ("COIL_NOCOND_2", {0: -j}), ("COIL_NOCOND_3", {0: -b * j, 1: b * j})):
        selection.select(mesh, group, 1, location="elem")
        for c, value in components.items():
            selection.set_value_on_selection(mesh, "Js", value, c, location="elem")

side = meshes.array_values(mesh, "COIL_COND_1", "elem") == 1          # side 1 of the passive loop, current along y
cell_area = 1e-4
rows = []
for n in range(0, 9):
    solver.time_step["time"] = 1e-4 * (n + 1)
    solver.time_step["transient_effects"] = n > 0
    set_source(1e6 if n == 0 else 0.0)
    solver.solve()
    jy = meshes.array_values(mesh, "Jt", "elem")[side, 1]
    rows.append((n * 1e-4, 100.0 if n == 0 else 0.0, jy.mean() * cell_area))
    if n == 1:
        solver.write_output()          # keep the step right after the switch-off for the 3D view
rows = np.array(rows)
for t, i_s, i_p in rows:
    print(f"t = {t * 1e3:.1f} ms: source {i_s:5.1f} A, passive loop {i_p:8.3f} A")
''')
    nb.code('''
fig, ax = plt.subplots(figsize=(6, 3.4))
ax.step(rows[:, 0] * 1e3, rows[:, 1], where="post", label="source loop")
ax.plot(rows[:, 0] * 1e3, rows[:, 2], "o-", label="passive loop")
ax.set_xlabel("t [ms]"); ax.set_ylabel("current [A]"); ax.grid(alpha=0.3); ax.legend(); plt.show()
''')
    nb.md("""
## The whole system, mirrored

The quarter model is reflected about x = 0 and y = 0 for the picture (the current of the reflected copies follows from
the symmetry). Shown right after the switch-off, read back from the output file: the induced current in the passive
loop (arrows: its direction, which keeps the flux of the switched-off source), the source loop in grey, and |B| on the
plane z = 0 between the two loops.
""")
    nb.code('''
from pytec3d import Model
out = Model.read_model(sorted(Path("output").glob("3D_loop_cond_base_out_*.vtu"))[0], reset_arrays=False).mesh
grid = pv.wrap(out).copy()
grid.cell_data["|Jt|"] = np.linalg.norm(grid.cell_data["Jt"], axis=1)
grid.point_data["|B|"] = np.linalg.norm(grid.point_data["B"], axis=1)
quarter = [grid]
quarter.append(grid.reflect((1, 0, 0), point=(0, 0, 0)))
quarter += [q.reflect((0, 1, 0), point=(0, 0, 0)) for q in list(quarter)]

pl = pv.Plotter(window_size=(850, 650))
limit = max(q.threshold(value=(24, 24), scalars="TYPE")["|Jt|"].max() for q in quarter)
for q in quarter:
    passive = q.threshold(value=(24, 24), scalars="TYPE")
    pl.add_mesh(q.threshold(value=(23, 23), scalars="TYPE"), color="lightgray")
    pl.add_mesh(passive, scalars="|Jt|", cmap="inferno", clim=[0, limit], show_scalar_bar=q is grid,
                scalar_bar_args={"title": "|Jt| passive", "vertical": True, "position_x": 0.03, "position_y": 0.25, "height": 0.5, "width": 0.06})
    pl.add_mesh(passive.cell_centers().glyph(orient="Jt", scale=False, factor=0.03), color="cyan")
    pl.add_mesh(q.slice(normal="z", origin=(0, 0, 0)), scalars="|B|", cmap="viridis", opacity=0.5,
                show_scalar_bar=q is grid, scalar_bar_args={"title": "|B| [T]", "vertical": True, "position_x": 0.88, "position_y": 0.25, "height": 0.5, "width": 0.06})
pl.camera_position = [(0.25, -0.3, 0.22), (0, 0, 0), (0, 0, 1)]
pl.show()
with quietly():
    solver.close()
''')
    return nb


if __name__ == "__main__":
    target = os.path.abspath(sys.argv[1])     # the notebooks change the working directory
    which = sys.argv[2:] or ["stranded_coil_iron", "coil_3d_sector", "edge_elements", "loops_3d"]
    for name in which:
        nbbuild.run(globals()[name](), f"{target}/tutorial_{name}.ipynb")

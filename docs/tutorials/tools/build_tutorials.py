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

CASE = Path("../../examples/{case}").resolve()
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

This notebook follows `examples/coil_2D_avoltemf_ferr`: a 1000-turn **stranded coil** fed by a voltage source, a passive
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
4. the coil terminal: `CP_0 = 1` on the coil nodes, coupling `["EMF", "VOLT"]` (`add_dof_couplings`);
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
    solver.model.add_dof_couplings({1: ["EMF", "VOLT"]})
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

This notebook follows `examples/coil_3D_avoltemf`: the coil and the passive ring of the 2D examples, now in 3D, modelled
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
solver.model.add_dof_couplings({1: ["VOLT", "EMF"]})          # the coil terminal
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
* Run the edge-element version (`examples/coil_3D_avoltemf_edge`): no frame rotation is needed, a tree-cotree gauge is.
""")
    return nb


# ======================================================================================================================
def edge_elements():
    nb = nbbuild.Notebook()
    nb.md(r"""
# Edge elements: pre- and post-processing of a 3D coil sector

This notebook takes the example `examples/coil_3D_avoltemf_edge` through a complete workflow around its macro:

1. **pre-processing**: look at the mesh, the regions, the model data and the edge DOFs before running anything;
2. **run** the example macro as it is, together with its nodal counterpart `examples/coil_3D_avoltemf`;
3. **post-processing**: read the outputs, compute the coil circuit quantities, the eddy current of the passive ring
   and the field, and compare the two formulations of the vector potential.

The model is the coil and passive ring of the 2D examples, revolved about the y axis and modelled on a 30 degree
cyclic sector. The lower ring is a 2000-turn stranded coil fed through its terminals (VOLT, EMF): first by a current
of 100 A, then by a voltage of 10 kV. The upper ring is a passive massive conductor (A and V).

With `AV3D_edge` the vector potential is discretized with **edge** (Nedelec) elements: its DOFs are the circulations
of A along the mesh edges, stored on the midside nodes of a quadratic mesh, while V, VOLT and EMF stay nodal. There is
no Coulomb gauge term: a tree-cotree decomposition of the edges fixes the gauge. The nodal model `AV3D_nodal` solves
the same problem with a nodal A and a penalty gauge. Theory: *Element library* (AV3D_edge, AV3D_nodal); user manual:
*Boundary conditions* (gauging edge elements, cyclic symmetry).
""")
    nb.code(SETUP.format(case="coil_3D_avoltemf_edge",
                         files='("coil_3D_avoltemf_edge.vtu", "model_data.json", "macro.py")')
            + "\nimport json, subprocess, sys, time\nfrom pytec3d import element_library\n")
    nb.md("""
## 1. Pre-processing: the mesh and the regions

The mesh is read without loading it into a solver: `meshes.read_mesh` gives the VTK grid with its arrays (TYPE,
MATERIAL, CONSTANT, the components of the model). The edge DOFs need a quadratic mesh (HEX20 / WEDGE15): their values
live on the midside nodes, flagged by the point array `MidsideNodes`.
""")
    nb.code('''
grid = pv.wrap(meshes.read_mesh("coil_3D_avoltemf_edge.vtu"))
print(f"{grid.n_points} points, {grid.n_cells} cells:", dict(zip(*np.unique(grid.celltypes, return_counts=True))),
      "(25 = HEX20, 26 = WEDGE15)")
print("midside nodes:", int(np.sum(grid.point_data["MidsideNodes"])) if "MidsideNodes" in grid.point_data else "-")
data = json.load(open("model_data.json"))
families = dict(enumerate(element_library.element_types())) if hasattr(element_library, "element_types") else {}
domains = dict(enumerate(element_library.domains()))
for tid, row in sorted(data["element_types"].items(), key=lambda kv: int(kv[0])):
    n = int(np.sum(grid.cell_data["TYPE"] == int(tid)))
    print(f"TYPE {tid}: element_type {row['element_type']} {families.get(row['element_type'], '')}, "
          f"{domains[row['option1']]:26s} {n:6d} cells")
for mid, mat in sorted(data["materials"].items(), key=lambda kv: int(kv[0])):
    print(f"MATERIAL {mid}: mur = {mat['mur'][0][0]:g}, rho = {mat['rho'][0][0]:g} Ohm m")
coil = data["constants"]["1"]
print("coil real constants:", {k: coil[k] for k in coil if k in ("cross_section", "turns", "resistance", "volume", "symmetry", "T")})
''')
    nb.code('''
types = grid.cell_data["TYPE"]
pl = pv.Plotter(window_size=(900, 520))
air = grid.threshold(value=(float(types.min()), float(types.min())), scalars="TYPE")
pl.add_mesh(air.extract_feature_edges(), color="lightgray")
for tid, color, label in ((int(t), c, l) for t, c, l in zip(sorted(set(types.tolist()))[1:], ("#c87533", "dimgray", "#5a8fbf"),
                                                          ("coil", "passive ring", "other"))):
    part = grid.threshold(value=(tid, tid), scalars="TYPE")
    if part.n_cells:
        pl.add_mesh(part, color=color, show_edges=True, edge_color="black", line_width=0.3, label=f"TYPE {tid}: {label}")
pl.add_legend(bcolor="white", size=(0.3, 0.12))
pl.camera_position = [(1.6, 0.9, 1.2), (0.35, 0.0, -0.1), (0, 1, 0)]
pl.show()
''')
    nb.md("""
The 30 degree sector: the coil (lower ring) and the passive ring (upper ring) in the air box. The macro sets the
boundary conditions on this grid: cyclic symmetry between the faces at 0 and 30 degrees (edge values are scalars: no
frame rotation), tangential A = 0 on the outer boundary (component 2 of `_A` on the midside nodes), and the
tree-cotree gauge (`tree_gauge.tree_gauge`), which writes A = 0 on the tree edges into `_A`.

## 2. Run the examples

The macros are run as they are, each in its own copy of the example folder. `coil_3D_avoltemf` is the same model with
nodal elements (HEX8 / WEDGE6, nodal frames rotated on the second cyclic face).
""")
    nb.code('''
NODAL = WORK / "coil_3D_avoltemf"
NODAL.mkdir()
for name in ("coil_3D_avoltemf.vtu", "model_data.json", "macro.py"):
    shutil.copy(CASE.parent / "coil_3D_avoltemf" / name, NODAL / name)
for folder in (WORK, NODAL):
    started = time.perf_counter()
    run = subprocess.run([sys.executable, "macro.py"], cwd=folder, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-3000:]
    print(f"{folder.name:24s} solved in {time.perf_counter() - started:5.1f} s:",
          sorted(p.name for p in (folder / "output").glob("*_out_*.vtu")))
''')
    nb.md("""
## 3. Post-processing

Each step wrote `output/<case>_out_<n>.vtu` (the solved model: point arrays A, V, VOLT, EMF, CURR, cell arrays B,
Jt, Js) and a line of `<case>_steps.csv`. We read both runs and compare:

* the **circuit** of the coil: VOLT and EMF on its nodes (all the coil nodes share one value through the coupling
  `CP_0 = 1`) and the coil current from the circuit equation, $I = (\\mathrm{VOLT} - \\mathrm{EMF}) / R$ (the
  array `CURR` holds the imposed current, a load: 0 on the voltage-driven step);
* the toroidal **eddy current** of the passive ring, $I = \\sum J_\\varphi\\,a$ over a meridian cut, with the cell area
  $a = V_{cell} / (r\\,\\Delta\\theta)$ (the sector has one cell in theta per element layer, 30 degrees);
* the **field** $|B|$ on the meridian plane in the middle of the sector.
""")
    nb.code('''
def results(folder, case, step):
    g = pv.wrap(meshes.read_mesh(str(folder / "output" / f"{case}_out_{step}.vtu"), reset_arrays=False))
    out = {"grid": g}
    coil_nodes = g.point_data["CP_0"] == 1
    for name in ("VOLT", "EMF"):
        values = g.point_data[name][coil_nodes]
        out[name] = float(values[np.argmax(np.abs(values))]) if values.size else 0.0
    out["I_coil"] = (out["VOLT"] - out["EMF"]) / data["constants"]["1"]["resistance"]
    c = g.cell_data["cell_centroids"]
    r = np.hypot(c[:, 0], c[:, 2])
    theta = np.arctan2(-c[:, 2], c[:, 0])
    j_phi = (g.cell_data["Jt"][:, 0] * c[:, 2] - g.cell_data["Jt"][:, 2] * c[:, 0]) / r   # e_phi = (z, 0, -x) / r
    passive = g.cell_data["TYPE"] == np.unique(g.cell_data["TYPE"])[2]
    out["I_passive"] = float(np.sum(j_phi[passive] * np.abs(g.cell_data["cell_volume"][passive]) / (r[passive] * np.radians(30.0))))
    return out

rows = []
for step in (1, 2):
    for label, folder, case in (("edge", WORK, "coil_3D_avoltemf_edge"), ("nodal", NODAL, "coil_3D_avoltemf")):
        res = results(folder, case, step)
        rows.append((step, label, res))
        print(f"step {step} {label:6s} VOLT = {res['VOLT']:10.2f} V   EMF = {res['EMF']:10.2f} V   "
              f"I_coil = {res['I_coil']:8.3f} A   I_passive = {res['I_passive']:9.2f} A")
''')
    nb.md("""
Step 1 imposes the coil current (100 A): VOLT and EMF follow from the circuit equation, VOLT = R I + EMF. Their
values are about 40 MV (the current rises in 10 microseconds), so their 2 kV difference carries the rounding of the
stored values: the nodal run shows 100.18 A instead of 100 A. Step 2 imposes the terminal voltage (10 kV): the
current follows. The passive ring carries the eddy current induced by the
change of the coil flux, opposite to the coil current. Both formulations give the same circuit quantities and the same
ring current, to the discretization error: the edge model has more A DOFs on a mesh of the same cells.
""")
    nb.code('''
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, (label, folder, case) in zip(axes, (("edge (AV3D_edge)", WORK, "coil_3D_avoltemf_edge"),
                                              ("nodal (AV3D_nodal)", NODAL, "coil_3D_avoltemf"))):
    g = pv.wrap(meshes.read_mesh(str(folder / "output" / f"{case}_out_2.vtu"), reset_arrays=False))
    a = np.radians(15.0)
    cut = g.slice(normal=(-np.sin(a), 0.0, -np.cos(a)), origin=(0.0, 0.0, 0.0))
    pts = cut.points
    r, y = np.hypot(pts[:, 0], pts[:, 2]), pts[:, 1]
    cut = cut.cell_data_to_point_data()
    b = np.linalg.norm(cut.point_data["B"], axis=1)
    tri = cut.triangulate()
    faces = tri.faces.reshape(-1, 4)[:, 1:]
    m = ax.tripcolor(r, y, faces, b, shading="gouraud", cmap="viridis")
    fig.colorbar(m, ax=ax, label="|B| [T]")
    ax.set_aspect("equal"); ax.set_xlim(0, 1.0); ax.set_ylim(-1.0, 1.0)
    ax.set_xlabel("r [m]"); ax.set_ylabel("y [m]"); ax.set_title(f"{label}, step 2")
fig.tight_layout()
''')
    nb.md("""
The field of the voltage-driven coil on the meridian plane at 15 degrees, for both formulations: the same pattern,
concentrated around the coil and coupled to the passive ring. Edge elements represent the tangential continuity of A
exactly and need no gauge penalty; nodal elements are cheaper per cell. To go further: change the voltage or the time
step in the macro and compare again; add iron (`examples/coil_3D_avoltemf_edge_ferr`).
""")
    return nb


if __name__ == "__main__":
    target = os.path.abspath(sys.argv[1])     # the notebooks change the working directory
    which = sys.argv[2:] or ["stranded_coil_iron", "coil_3d_sector", "edge_elements"]
    for name in which:
        nbbuild.run(globals()[name](), f"{target}/tutorial_{name}.ipynb")

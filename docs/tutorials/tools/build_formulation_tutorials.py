"""Build the tutorial notebooks of the scalar-potential formulations (examples core_3D_phi, ring_3D_tphi,
coil_core_3D_tphi, wire_2D_tphi, bulk_2D_tphi_axisym), with their outputs.

Usage, from docs/tutorials/tools, with the pytec3d environment (pyvista installed):
    python build_formulation_tutorials.py ..                     # all of them (about 10 minutes)
    python build_formulation_tutorials.py .. tphi_cuts           # one of them (scalar_potential, tphi_cuts, superconductors)
"""
import os
import sys

import nbbuild

SETUP = '''
import contextlib, io, os, shutil, tempfile
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

import pytec3d
from pytec3d import Solver, selection, tree_gauge
from pytec3d import mesh as meshes

pv.set_jupyter_backend("static")          # images in the notebook; "trame" gives interactive views
pytec3d.set_log_level("ERROR")            # the job log keeps everything (<model>.log in the working folder)
MU0 = 4e-7 * np.pi
QUIET = {{"-ksp_monitor": None, "-ksp_view": None, "-ksp_converged_reason": None, "-log_view": None,
         "-ksp_monitor_true_residual": None, "-ksp_monitor_solution": None}}   # no PETSc reports in the notebook

EXAMPLES = Path("../../examples").resolve()
WORK = Path(tempfile.mkdtemp(prefix="pytec3d_{name}_"))
for case in {cases}:
    folder = EXAMPLES / case
    for f in (case + ".vtu", "model_data.json"):
        shutil.copy(folder / f, WORK / f"{{case}}_{{f}}" if f == "model_data.json" else WORK / f)
os.chdir(WORK)
print("working in a temporary folder:", WORK.name)

def quietly():
    """Keep the report tables of the solver (and libMesh's Newton trace) out of the notebook: they are in the job log."""
    return contextlib.redirect_stdout(io.StringIO())

def new_solver(case, log):
    """A solver on the example's mesh and model data (model_data.json copied as <case>_model_data.json)."""
    solver = Solver(log_file=log)
    solver.read_model(case + ".vtu")
    solver.model.read_json(case + "_model_data.json")
    solver.solver_settings = {{"petsc_options": QUIET}}
    return solver

def select_surface(grid, conductor, others, midside):
    """Select the points shared by the cells of TYPE ``conductor`` and those of TYPEs ``others`` (the surface of the
    conductor): its midside points (midside=True: the edge DOFs) or its corner points."""
    selection.select(grid, "TYPE", conductor, location="elem"); selection.nsle(grid)
    selection.select(grid, "TYPE", list(others), location="elem"); selection.nsle(grid, mode="r")
    selection.select(grid, "MidsideNodes", int(midside), mode="r")

def select_inside(grid, conductor, others):
    """Select the corner points of the cells of TYPE ``conductor`` that no cell of TYPEs ``others`` has."""
    selection.select(grid, "TYPE", conductor, location="elem"); selection.nsle(grid)
    selection.select(grid, "TYPE", list(others), location="elem"); selection.nsle(grid, mode="u")
    selection.select(grid, "MidsideNodes", 0, mode="r")
'''


# ======================================================================================================================
def scalar_potential():
    nb = nbbuild.Notebook()
    nb.md(r"""
# Magnetostatics with the scalar potential: a coil on a closed iron frame (PHI3D_nodal)

This notebook follows `examples/core_3D_phi`. A coil is wound on one leg of a closed square iron frame; we solve the
static field with the element type **PHI3D_nodal**:

$$\mathbf{H}=\mathbf{H}_s-\nabla\phi,\qquad \nabla\times\mathbf{H}_s=\mathbf{J}_s,\qquad
\nabla\cdot\mu(\mathbf{H}_s-\nabla\phi)=0 ,$$

with the reduced scalar potential $\phi$ (`PHI`, nodal) and the source field $\mathbf{H}_s$ (`Hs`, edge elements)
both unknowns of **one** system. Why this model:

* an iron loop around a current is the case where the classical reduced potential (a Biot-Savart $\mathbf{H}_s$ given
  per cell) loses its accuracy: in the iron $\mathbf{H}\approx0$, so $\mathbf{H}_s$ and $\nabla\phi$ nearly cancel, and
  if they are not in the same discrete space the cancellation error dominates. Here they are: the gradients of
  first-order Lagrange functions are edge (Nedelec) functions;
* $\phi$ stays single valued: $\mathbf{H}_s$ carries the current, no cut is needed.

The notebook shows:

1. the model, the source current `Js` and the gauge of $\mathbf{H}_s$;
2. the magnetization curve of the frame with a B-H curve (Newton, static load steps), against a linear iron and a
   magnetic-circuit estimate;
3. the alternative of an $\mathbf{H}_s$ **given on the edges** (the circulations of a known field: only $\phi$ is
   solved);
4. the field in the saturated frame.

See the user manual (*Element library: PHI3D_nodal*, *Nonlinear laws*) and the theory (*Scalar potentials*).
""")
    nb.code(SETUP.format(name="core_3D_phi", cases='("core_3D_phi",)'))
    nb.md("""
## The model

The box |x|, |y| < 2 m, |z| < 1.5 m, 3072 HEX20 cells of 0.25 m (`examples/files/formulation_meshes.py`). The frame
(TYPE 2, `SOFT_MAGNETIC`) is 0.5 < max(|x|, |y|) < 1 m, |z| < 0.5 m: legs of 0.5 x 1 m section, a mean magnetic path
of about 6 m. The coil (TYPE 3, `STRANDED_CONDUCTOR`) is a 0.25 m thick winding around the leg 0.5 < x < 1 m, over
|y| < 0.5 m. All three TYPE rows are `PHI3D_nodal` with `option2` = 1: the type is first order only.
""")
    nb.code('''
import json
print(json.dumps(json.load(open("core_3D_phi_model_data.json"))["element_types"]))
grid = pv.wrap(meshes.read_mesh("core_3D_phi.vtu"))
print(grid.n_cells, "cells, cells per TYPE:", dict(zip(*np.unique(grid.cell_data["TYPE"], return_counts=True))))

p = pv.Plotter(window_size=(700, 500))
for t, color, label in ((2, "silver", "iron"), (3, "orange", "coil")):
    part = grid.threshold((t, t), scalars="TYPE")
    p.add_mesh(part, color=color, label=label)
    p.add_mesh(part.extract_feature_edges(), color="black")
p.add_mesh(grid.outline(), color="black")
p.add_legend(); p.view_isometric(); p.show()
''')
    nb.md(r"""
## The source current and the gauge of Hs

`Js` is a cell array of the `STRANDED_CONDUCTOR` cells. The current circulates about +y around the leg; per
ampere-turn its density is $J_0=1/(0.25\times1.0)$ A/m². The corner cells of the winding take half of each direction:
that cell-wise `Js` is not exactly divergence-free, and need not be. The weak form of $\mathbf{H}_s$,
$\int\nabla\times\mathbf{h}\cdot\nabla\times\mathbf{H}_s=\int\nabla\times\mathbf{h}\cdot\mathbf{J}_s$, only sees the
solenoidal part of `Js`; the element output `Jt` $=\nabla\times\mathbf{H}_s$ shows the current the source field really
carries.

$\mathbf{H}_s$ is a curl-curl unknown: it needs a gauge. A spanning **tree** of all the mesh edges, grown from the box,
has its edges fixed to zero (`tree_gauge(..., array_name="_Hs")`), and $\mathbf{n}\times\mathbf{H}_s=0$ on the box
(`_Hs` = 0 on the box edges). Edge values live on the midside nodes, component 2 of `_Hs`. $\phi$ needs one fixed value
(the box is a free boundary: a flux wall, $\mathbf{n}\cdot\mathbf{B}=0$).
""")
    nb.code('''
BH_H = np.concatenate(([0.0], np.geomspace(1.0, 1e8, 400)))
BH = np.column_stack((BH_H, MU0 * BH_H + 1.8 * BH_H / (BH_H + 300.0)))     # B = mu0 H + Bs H / (H + Hk)

def coil_current(centroids, coil):
    """Js per ampere-turn: about +y around the leg x in (0.5, 1), |z| < 0.5."""
    u, w = centroids[:, 0] - 0.75, centroids[:, 2]
    du, dw = np.abs(u) - 0.25, np.abs(w) - 0.5
    along_z = np.where(np.isclose(du, dw), 0.5, (du > dw).astype(float))
    j = np.zeros((len(centroids), 3))
    j[coil, 0] = np.sign(w[coil]) * (1.0 - along_z[coil])
    j[coil, 2] = -np.sign(u[coil]) * along_z[coil]
    return j / 0.25

def t0(p, ni=1.0):
    """A source field with curl t0 = the winding current (per ampere-turn): t0 = f(x, z) e_y for |y| < 0.5, f = ni / 1
    over the leg section, falling linearly to 0 across the 0.25 m of the winding (piecewise linear on the grid)."""
    ramp = lambda d, inner: np.clip((inner + 0.25 - d) / 0.25, 0.0, 1.0)
    f = ni * ramp(np.abs(p[:, 0] - 0.75), 0.25) * ramp(np.abs(p[:, 2]), 0.5) * (np.abs(p[:, 1]) < 0.5)
    return np.column_stack((np.zeros(len(p)), f, np.zeros(len(p))))

def run(ni_list, bh=True, hs_given=False, log="core.log"):
    """The frame at each NI of ni_list (static steps). bh: B-H curve, else mur = 1000. hs_given: _Hs = the edge
    circulations of t0 on every edge (no Js, no tree), else Js and the gauge. Returns the mean |B_y| of the opposite
    leg, the Newton iterations, the free DOFs and the grid of the last step."""
    solver = new_solver("core_3D_phi", log)
    grid = solver.model.mesh
    if bh:
        solver.model.set_bh_curve(2, BH)
    selection.select_all(grid); selection.nsel_ext(grid)
    box = selection.selected_ids(grid, "node")
    info = None
    if hs_given:
        mids, start, end = tree_gauge.edges(grid)
        circulation = tree_gauge.edge_circulations(grid, start, end, t0)
    else:
        info = tree_gauge.tree_gauge(grid, array_name="_Hs")       # selected cells: all; roots: the box points
        selection.set_value_on_selection(grid, "_Hs", 0.0, 2)      # n x Hs = 0 on the box
    selection.select(grid, "ID", int(box[0]))
    selection.set_value_on_selection(grid, "_PHI", 0.0)
    centroids, volumes = meshes.cell_centroids(grid), meshes.cell_sizes(grid)
    types = meshes.array_values(grid, "TYPE", "elem")
    coil = np.flatnonzero(types == 3)
    j_unit = coil_current(centroids, types == 3)[coil]
    leg = (types == 2) & (centroids[:, 0] < -0.5) & (np.abs(centroids[:, 1]) < 0.5)
    with quietly():
        solver.load_model_into_solver()
    solver.time_step["transient_effects"] = False
    b_leg, newton = [], []
    for step, ni in enumerate(ni_list, start=1):
        if hs_given:                                               # one value per edge (midside point)
            selection.select(grid, "ID", mids)
            selection.set_value_on_selection(grid, "_Hs", ni * circulation, 2)
        else:                                                      # one vector per coil cell
            selection.select(grid, "ID", coil, location="elem")
            selection.set_value_on_selection(grid, "Js", ni * j_unit, None, location="elem")
        solver.time_step["time"] = float(step)
        with quietly():
            solver.solve()
        b = meshes.array_values(grid, "B", "elem")
        b_leg.append(abs(volumes[leg] @ b[leg, 1]) / volumes[leg].sum())
        newton.append(solver.steps[-1]["newton_iterations"])
    out = {"B": np.array(b_leg), "newton": newton, "info": info, "grid": pv.wrap(grid).copy()}
    with quietly():
        solver.close()
    return out
''')
    nb.md("""
## The magnetization curve of the frame

Static steps from 250 A-turns to 50 kA-turns with the B-H curve B = mu0 H + 1.8 H / (H + 300) (initial mu_r about
4800, saturating at 1.8 T), each solved by Newton from the previous step, and the same with a linear iron of
mu_r = 1000. A magnetic-circuit estimate: H = NI / l with the mean path l = 6 m, B = B(H) of the curve (no leakage).
""")
    nb.code('''
ni = np.array([250.0, 500.0, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4])
nonlinear = run(ni, bh=True, log="nonlinear.log")
linear = run(ni[:4], bh=False, log="linear.log")
info = nonlinear["info"]
print(f"Hs gauge: {info['num_edges']} edges, {info['num_tree']} tree edges fixed, {info['num_roots']} root(s)")
for k, (n, b, it) in enumerate(zip(ni, nonlinear["B"], nonlinear["newton"])):
    print(f"NI = {n:7.0f} A-turns: B_leg = {b:.3f} T, Newton iterations {it}")

estimate = np.interp(ni / 6.0, BH[:, 0], BH[:, 1])
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.semilogx(ni, nonlinear["B"], "o-", label="PHI3D_nodal, B-H curve")
ax.semilogx(ni[:4], linear["B"], "s--", label="PHI3D_nodal, mu_r = 1000")
ax.semilogx(ni, estimate, "k:", label="magnetic circuit, l = 6 m")
ax.set_xlabel("NI (A-turns)"); ax.set_ylabel("mean |B_y| of the opposite leg (T)"); ax.grid(True, which="both", alpha=0.3)
ax.legend(); plt.show()
''')
    nb.md(r"""
The flux in the opposite leg follows the magnetic-circuit estimate closely: below the knee the steep part of the
curve, above 10 kA-turns the saturation (the curves flatten at about 1.7 T; the leakage through the window and the
corners make the difference). The linear iron (mu_r = 1000, below the initial 4800 of the curve) grows without limit. Newton needs only a few
iterations per step (the law `bh_curve_h`: input H, output $\mu$, its tangent $d\mathbf{B}/d\mathbf{H}$).

## Hs given on the edges

Instead of solving $\mathbf{H}_s$, it may be **given**: `_Hs` on every edge, the circulations of any field with the
right curl, for example the Biot-Savart field of the coils, here the piecewise-linear field `t0` of the winding
($\mathbf{t}_0=f\,\mathbf{e}_y$, $f$ = NI over the leg section, falling to 0 across the winding: its curl is the coil
current). Circulations, not point values: their discrete curl is then exactly the flux of the current through each
mesh face (Stokes), and the cancellation in the iron stays exact. All the $\mathbf{H}_s$ rows are then eliminated:
only $\phi$ is solved (no tree, no `Js`).
""")
    nb.code('''
given = run(ni, bh=True, hs_given=True, log="given.log")
for n, b1, b2 in zip(ni, nonlinear["B"], given["B"]):
    print(f"NI = {n:7.0f}: B_leg solved Hs {b1:.4f} T, given Hs {b2:.4f} T, difference {100 * (b2 / b1 - 1):+.2f} %")
''')
    nb.md(r"""
The two agree within about a percent: the same current, two source fields differing by a gradient, which $\phi$
absorbs (up to the discretization: the solved $\mathbf{H}_s$ carries the solenoidal part of the cell-wise `Js`, the
given one exactly the winding current). The benchmark `closed_core` of the `formulations` suite compares both against
the A formulation: 1.9 % and 3.6 % on this mesh size.

Which to use: solving $\mathbf{H}_s$ needs only `Js` (any shape, any number of coils) and adds one unknown per edge of
the whole mesh; giving it needs a field whose circulations you can compute (a Biot-Savart integral per edge), and
removes those unknowns.

## The field in the saturated frame

|B| on the mid-plane z = 0 at 50 kA-turns: the legs are saturated (1.6-1.7 T); the outer corners of the frame carry
less flux (it takes the inner, shorter path), and some flux leaks into the window and around the coil (orange: the
coil).
""")
    nb.code('''
g = nonlinear["grid"]
g.cell_data["|B|"] = np.linalg.norm(g.cell_data["B"], axis=1)
p = pv.Plotter(window_size=(700, 560))
p.add_mesh(g.slice(normal="z", origin=(0, 0, 0.01)), scalars="|B|", cmap="viridis", clim=(0, 2.0),
           scalar_bar_args={"title": "|B| (T), z = 0, NI = 50 kA"})
p.add_mesh(g.threshold((3, 3), scalars="TYPE").extract_feature_edges(), color="orange")
p.view_xy(); p.show()
''')
    nb.md("""
## Takeaways

* `PHI3D_nodal`: `Js` on `STRANDED_CONDUCTOR` cells, a tree gauge of `_Hs` rooted on the box plus `_Hs` = 0 on the
  box edges, one fixed `_PHI`; first order (`option2` = 1); static steps (`transient_effects` False).
* An iron loop around a current needs nothing special: no cut, no cancellation error.
* A B-H curve on the iron is the law `bh_curve_h` (`set_bh_curve` selects it for this type); Newton converges in a
  few iterations per load step.
* `_Hs` given on every edge (`tree_gauge.edges`, `tree_gauge.edge_circulations` of a known field) is the cheap
  variant when the source field is known.
""")
    return nb


# ======================================================================================================================
def tphi_cuts():
    nb = nbbuild.Notebook()
    nb.md(r"""
# T-phi: cuts and fixed-T coils (TPHI3D_edge)

This notebook follows `examples/ring_3D_tphi` and `examples/coil_core_3D_tphi`. In the T-phi formulation

$$\mathbf{H}=\mathbf{T}-\nabla\phi,\qquad \mathbf{J}=\nabla\times\mathbf{T},$$

$\phi$ (`PHI`, nodal) lives everywhere and $\mathbf{T}$ (`T`, edge elements) only in the conductors and in
`CUT_LAYER` cells. There is no source current and no circuit: currents are carried by $\mathbf{T}$. It shows:

1. **why a cut** is needed for a conductor with a net current, and how to build one (a ring);
2. the three **drives** of a cut: an induced current, an imposed current, an imposed voltage, and their consistency
   (the same resistance and time constant);
3. what happens **without the cut**;
4. a current-driven **coil as T fixed on its edges** around a conducting core.

See the user manual (*Element library: TPHI3D_edge* and *Cuts*) and the theory (*Scalar potentials*).
""")
    nb.code(SETUP.format(name="tphi", cases='("ring_3D_tphi", "coil_core_3D_tphi")'))
    nb.md(r"""
## Why a cut

Outside the conductors $\mathbf{H}=-\nabla\phi$ with a single-valued $\phi$: its circulation along any closed loop is
zero. A loop around a ring carrying the current I must give I (Ampère): the air must be made simply connected by a
**cut**, a surface spanning the hole on which $\mathbf{T}$ carries the jump. In practice: a one-cell layer of
`CUT_LAYER` cells on one side of the surface, with $\mathbf{T}=I\nabla c$, c going from 0 to 1 across the layer:

* the edges in the two faces of the layer (*planar*) have $\int\mathbf{T}\cdot d\mathbf{l}=0$;
* the edges *crossing* the layer have $\int\mathbf{T}\cdot d\mathbf{l}=\pm I$ (+ when the solver's orientation of the
  edge goes from the face c = 0 to the face c = 1).

The ring of the example: 0.5 < max(|x|, |y|) < 1 m, |z| < 0.25 m (TYPE 2), the layer max(|x|, |y|) < 0.5 m,
0 < z < 0.25 m (TYPE 3). `tree_gauge.edges(grid, cells)` lists the edges of the layer cells with their midside point
and their end points in the solver's orientation (from the lexicographically smaller end point, x then y then z):
here the crossing edges are the vertical ones, all oriented along +z.
""")
    nb.code('''
solver = new_solver("ring_3D_tphi", "look.log")
grid = solver.model.mesh
xyz = meshes.points(grid)
layer_cells = np.flatnonzero(meshes.array_values(grid, "TYPE", "elem") == 3)
mids, start, end = tree_gauge.edges(grid, layer_cells)
crossing = ~np.isclose(xyz[start, 2], xyz[end, 2])
print(f"{len(layer_cells)} layer cells, {len(mids)} edges: {crossing.sum()} crossing, {(~crossing).sum()} planar")
print("crossing edges oriented along +z:", bool(np.all(xyz[end[crossing], 2] > xyz[start[crossing], 2])))

g = pv.wrap(grid)
p = pv.Plotter(window_size=(700, 500))
for t, color, label, opacity in ((2, "peru", "ring", 0.5), (3, "lightblue", "cut layer", 1.0)):
    part = g.threshold((t, t), scalars="TYPE")
    p.add_mesh(part, color=color, opacity=opacity, label=label)
    p.add_mesh(part.extract_feature_edges(), color="black")
p.add_points(xyz[mids[crossing]], color="red", point_size=10, render_points_as_spheres=True, label="crossing edges")
p.add_legend(); p.camera_position = [(3.0, -3.5, 2.5), (0, 0, 0), (0, 0, 1)]; p.show()
with quietly():
    solver.close()
''')
    nb.md(r"""
## One function, three drives

`ring(drive)` sets the model up as the example does, then runs 30 backward Euler steps of 10 ms:

* gauges of T-phi: $\mathbf{n}\times\mathbf{T}=0$ on the ring surface (`_T` = 0 on its edges), `_PHI` = 0 at the
  nodes inside the ring (first-order $\phi$: its gradients are in the edge space of $\mathbf{T}$);
* the cut: `_T` = 0 on the planar edges; the crossing edges, those on the ring's inner surface included, are taken out
  of the $\mathbf{n}\times\mathbf{T}=0$ set; then
  * `"induced"`: the crossing edges coupled on **one DOF**, the ring current (CP id 1 on their midside nodes, ids 2
    and 3 on the corner nodes of the two faces: they give the solver the orientation of each crossing edge); a
    uniform field of 1 T along z switched on at t = 0, the total potential $\phi=-H_0z$ on the box;
  * `"current"`: `_T` = I0 = 10 kA on every crossing edge (all along +z), no coupling, no field;
  * `"voltage"`: the coupling of `"induced"` and `VLTG` = V0 on one crossing midside node: the voltage around the
    ring, no field;
  * `"none"`: the field of `"induced"`, but **no cut**: the layer is air (TYPE 1).

Every value is set as in the example macros: select the points (or cells), then `set_value_on_selection`, with one
value for all or one value per selected id in increasing id order (e.g. the potential of each box point). It returns
the ring current from the cell values of `Jt` across a section of the ring, and the Joule power $\int\rho J^2$.
""")
    nb.code('''
RHO, DT, STEPS = 1e-6, 0.01, 30

def ring(drive, b0=1.0, i0=1e4, v0=0.24):
    solver = new_solver("ring_3D_tphi", f"ring_{drive}.log")
    grid = solver.model.mesh
    xyz = meshes.points(grid)
    types = meshes.array_values(grid, "TYPE", "elem")
    if drive == "none":
        selection.select(grid, "TYPE", 3, location="elem")
        selection.set_value_on_selection(grid, "TYPE", 1, location="elem")      # the layer becomes air
    select_surface(grid, 2, [1, 3], midside=True)
    selection.set_value_on_selection(grid, "_T", 0.0, 2)           # n x T = 0 on the ring surface
    select_inside(grid, 2, [1, 3])
    selection.set_value_on_selection(grid, "_PHI", 0.0)            # PHI = 0 inside the ring
    if drive != "none":
        selection.select(grid, "TYPE", 3, location="elem")
        mids, start, end = tree_gauge.edges(grid, selection.selected_ids(grid, "elem"))
        cross = ~np.isclose(xyz[start, 2], xyz[end, 2])
        selection.select(grid, "ID", mids[~cross])
        selection.set_value_on_selection(grid, "_T", 0.0, 2)       # planar edges
        selection.select(grid, "ID", mids[cross])                  # crossing edges: imposed, or free and coupled
        selection.set_value_on_selection(grid, "_T", i0 if drive == "current" else np.nan, 2)
        if drive in ("induced", "voltage"):
            selection.set_value_on_selection(grid, "CP_0", 1)
            corners = np.unique(np.concatenate((start, end)))
            for face, z in ((2, 0.0), (3, 0.25)):                  # the faces c = 0 and c = 1 of the layer
                selection.select(grid, "ID", corners[np.isclose(xyz[corners, 2], z)])
                selection.set_value_on_selection(grid, "CP_0", face)
            solver.model.add_dof_couplings({1: ["T"], 2: ["T"], 3: ["T"]})
        if drive == "voltage":
            selection.select(grid, "ID", int(mids[cross][0]))
            selection.set_value_on_selection(grid, "VLTG", v0)
    selection.select_all(grid); selection.nsel_ext(grid)
    box = selection.selected_ids(grid, "node")
    field = b0 if drive in ("induced", "none") else 0.0
    selection.set_value_on_selection(grid, "_PHI", -field / MU0 * xyz[box, 2])   # one value per box point
    with quietly():
        solver.load_model_into_solver()
    solver.time_step["theta"] = 1.0
    c, vol = meshes.cell_centroids(grid), meshes.cell_sizes(grid)
    section = (types == 2) & (c[:, 0] > 0) & (c[:, 1] > 0) & (c[:, 1] < 0.25)
    rows = []
    for step in range(1, STEPS + 1):
        solver.time_step["time"] = step * DT
        with quietly():
            solver.solve()
        j = meshes.array_values(grid, "Jt", "elem")
        rows.append((step * DT, j[section, 1].sum() * 0.25 * 0.25, RHO * (np.einsum("ij,ij->i", j, j) @ vol)))
        if step == 1:
            first = pv.wrap(grid).copy()                          # the field of the first step
    out = {"t": np.array([r[0] for r in rows]), "I": np.array([r[1] for r in rows]),
           "P": np.array([r[2] for r in rows]), "grid": first}
    with quietly():
        solver.close()
    return out

runs = {d: ring(d) for d in ("induced", "current", "voltage", "none")}
''')
    nb.md(r"""
## Resistance and time constant, three ways

* `"current"`: the current is imposed; once the transient eddy currents inside the ring section have died out, the
  power gives the resistance, $R=P/I_0^2$;
* `"voltage"`: the current rises to $V_0/R$ with the time constant $\tau=L/R$; with backward Euler,
  $I_k-I_\infty\propto(1+\Delta t/\tau)^{-k}$;
* `"induced"`: the current jumps at the field step, then decays as $(1+\Delta t/\tau)^{-k}$.

The three must give one ring: the same R and the same $\tau$.
""")
    nb.code('''
def tau_of(d):
    """tau from the step ratio of a backward-Euler exponential, d_{k+1} / d_k = 1 / (1 + dt / tau), steps 15 to 25."""
    ratio = (d[25] / d[15]) ** (1 / 10)
    return DT / (1.0 / ratio - 1.0)

r_current = runs["current"]["P"][-1] / 1e4**2
i_inf = 0.24 / r_current
iv = runs["voltage"]["I"]
tau_v = tau_of(np.diff(iv))                   # the increments of the rising current: no need of its final value
tau_i = tau_of(runs["induced"]["I"])
i_end = iv[-1] + np.diff(iv)[-1] / (DT / tau_v)   # the rising current extrapolated to t -> infinity
print(f"R from the imposed current      : {r_current * 1e6:.3f} uOhm   (rho x mean perimeter / section: {1e-6 * 6 / 0.25 * 1e6:.3f})")
print(f"voltage drive, I extrapolated   : {i_end / 1e3:.2f} kA  (V0 / R = {i_inf / 1e3:.2f} kA)")
print(f"tau from the voltage drive      : {tau_v * 1e3:.1f} ms")
print(f"tau from the induced current    : {tau_i * 1e3:.1f} ms")
print(f"L = tau R                       : {tau_i * r_current * 1e6:.3f} uH")

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
ax[0].plot(runs["induced"]["t"], -runs["induced"]["I"] / 1e3, "o-", ms=3, label="induced, with cut (-I)")
ax[0].plot(runs["none"]["t"], -runs["none"]["I"] / 1e3, "s-", ms=3, label="induced, no cut (-I)")
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("ring current (kA)"); ax[0].legend(); ax[0].grid(alpha=0.3)
ax[1].plot(runs["voltage"]["t"], runs["voltage"]["I"] / 1e3, "o-", ms=3, label="voltage drive V0 = 0.24 V")
ax[1].plot(runs["current"]["t"], runs["current"]["I"] / 1e3, "s-", ms=3, label="current drive I0 = 10 kA")
ax[1].axhline(i_inf / 1e3, color="k", ls=":", label="V0 / R")
ax[1].set_xlabel("t (s)"); ax[1].legend(); ax[1].grid(alpha=0.3)
plt.show()
''')
    nb.md(r"""
The two time constants agree, and the voltage drive tends to $V_0/R$ with R from the current drive: the cut DOF is
the net current of the ring and the voltage load `VLTG` is the voltage around it (the discrete energy balance
$P=V_0I$ holds at steady state; the benchmark `ring_transport` checks it). The 2 % left is the estimate of R: the
power from the cell means of J, $\sum\rho\bar J^2V$, is a little below $\int\rho J^2$, so R is low and $V_0/R$ high. The section current of the current drive
is $I_0$ at every step: an imposed current through a cut is exact. R is a little below
$\rho\,\ell/S$ with the mean perimeter: the current takes the inner, shorter path at the corners.

**Without the cut**, the field step induces only local eddy loops inside the ring section, with no net current: the
section current is about zero, the ring does not shield its hole. The solver gives no warning: a model of a
multiply connected conductor without its cut is simply another physical problem (a ring cut open). Check the net
current wherever it matters.

The field of the first step (t = 10 ms), on the plane y = 0: with the cut the ring current shields the hole;
without it the field goes through the hole, only the ring section itself is shielded by its local eddy loops.
""")
    nb.code('''
p = pv.Plotter(shape=(1, 2), window_size=(1000, 420))
for k, (name, title) in enumerate((("induced", "with the cut"), ("none", "without the cut"))):
    solver_grid = runs[name]["grid"]
    s = solver_grid.slice(normal="y", origin=(0, 0.01, 0))
    s.cell_data["B_z"] = s.cell_data["B"][:, 2]
    p.subplot(0, k)
    p.add_mesh(s, scalars="B_z", cmap="coolwarm", clim=(0.0, 2.0), scalar_bar_args={"title": f"B_z (T), {title}"})
    p.add_mesh(solver_grid.threshold((2, 2), scalars="TYPE").slice(normal="y", origin=(0, 0.01, 0)).extract_feature_edges(), color="black")
    p.view_xz()
p.show()
''')
    nb.md(r"""
## A coil as T fixed on its edges

T-phi has no source current: a current-driven coil is $\mathbf{T}$ **fixed** on its edges. Take any $\mathbf{T}_0$
with $\nabla\times\mathbf{T}_0=\mathbf{J}$ (per ampere-turn) that is zero outside the winding and the hole it encloses,
and impose `_T` = $NI(t)\int_e\mathbf{T}_0\cdot d\mathbf{l}$ on every edge of those cells. For the square solenoid
a < max(|x|, |y|) < b, |z| < h:

$$\mathbf{T}_0=f\,\mathbf{e}_z,\qquad f=\begin{cases}J\,(b-\max(|x|,|y|)) & \text{winding}\\ J\,(b-a) & \text{hole column}
\\ 0 & \text{elsewhere}\end{cases}\qquad |z|<h .$$

The *hole column* is "the disk filling the hole": the region the winding encloses, over the height of the winding.
There $\mathbf{T}_0$ is uniform (no current), and $\mathbf{n}\times\mathbf{T}_0=0$ on the boundary of winding + hole.
A **conductor in the hole** (the core of the example) has its own $\mathbf{T}$, which cannot be fixed: $\mathbf{T}_0$ is
moved out of it by a gradient, $\mathbf{T}_0'=\mathbf{T}_0-\nabla g$, $g=J(b-a)\,z\,\chi$, $\chi$ = 1 on the nodes of the
core cells (exact on the edges). Below: the edge values of $\mathbf{T}_0'$ on the plane y = 0.
""")
    nb.code('''
A_, B_, H_ = 0.5, 1.0, 0.5
J_UNIT = 1.0 / ((B_ - A_) * 2 * H_)

def t0(p):
    r = np.clip(np.maximum(np.abs(p[:, 0]), np.abs(p[:, 1])), A_, B_)
    return np.column_stack((np.zeros(len(p)), np.zeros(len(p)), J_UNIT * (B_ - r) * (np.abs(p[:, 2]) < H_)))

def coil_core(core=True, ni_max=5e5, t_ramp=0.1, dt=0.02, steps=15):
    """The solenoid ramped to ni_max over t_ramp; core False: the core cells are part of the hole (CUT_LAYER)."""
    solver = new_solver("coil_core_3D_tphi", f"coil_core_{core}.log")
    grid = solver.model.mesh
    xyz = meshes.points(grid)
    if not core:                                                   # the core cells become part of the hole
        selection.select(grid, "TYPE", 3, location="elem")
        selection.set_value_on_selection(grid, "TYPE", 4, location="elem")
        selection.set_value_on_selection(grid, "MATERIAL", 4, location="elem")
    types = meshes.array_values(grid, "TYPE", "elem")
    selection.select(grid, "TYPE", [2, 4], location="elem")
    mids, start, end = tree_gauge.edges(grid, selection.selected_ids(grid, "elem"))
    values = tree_gauge.edge_circulations(grid, start, end, t0)
    if core:
        selection.select(grid, "TYPE", 3, location="elem"); selection.nsle(grid)
        g = J_UNIT * (B_ - A_) * xyz[:, 2] * selection.selection_mask(grid, "node")   # chi: the core points
        values -= g[end] - g[start]
        select_inside(grid, 3, [1, 2, 4])
        selection.set_value_on_selection(grid, "_PHI", 0.0)
    selection.select_all(grid); selection.nsel_ext(grid)
    selection.select(grid, "ID", int(selection.selected_ids(grid)[0]))
    selection.set_value_on_selection(grid, "_PHI", 0.0)
    with quietly():
        solver.load_model_into_solver()
    solver.time_step["theta"] = 1.0
    c, vol = meshes.cell_centroids(grid), meshes.cell_sizes(grid)
    t0_cells = t0(c)
    centre = int(np.argmin(np.linalg.norm(xyz, axis=1)))
    rows = []
    for step in range(1, steps + 1):
        t = step * dt
        ni = ni_max * min(t / t_ramp, 1.0)
        selection.select(grid, "ID", mids)
        selection.set_value_on_selection(grid, "_T", ni * values, 2)
        solver.time_step["time"] = t
        with quietly():
            solver.solve()
        b = meshes.array_values(grid, "B", "elem")
        j = meshes.array_values(grid, "Jt", "elem")
        rows.append((t, ni, meshes.array_values(grid, "B", "node")[centre, 2],
                     float(np.einsum("ij,ij->i", t0_cells, b) @ vol),
                     np.linalg.norm(j[types == 3], axis=1).max() if core else 0.0))
    out = {k: np.array([r[i] for r in rows]) for i, k in enumerate(("t", "NI", "Bc", "Psi", "Jcore"))}
    out.update(mids=mids, values=values, xyz=xyz, grid=pv.wrap(grid).copy())
    with quietly():
        solver.close()
    return out

with_core, no_core = coil_core(True), coil_core(False)
m, v, x = with_core["mids"], with_core["values"], with_core["xyz"]
plane = np.isclose(x[m, 1], 0.0)
fig, ax = plt.subplots(figsize=(6.5, 3.6))
sc = ax.scatter(x[m[plane], 0], x[m[plane], 2], c=v[plane], cmap="viridis", s=18)
ax.add_patch(plt.Rectangle((-0.25, -0.25), 0.5, 0.5, fill=False, ec="red", label="core"))
ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)"); ax.legend(loc="upper right")
fig.colorbar(sc, label="edge value of T0' per A-turn (A)"); ax.set_title("T0' on the edges of the plane y = 0")
plt.show()
''')
    nb.md(r"""
The vertical edges of the winding carry $J(b-r)\,\Delta z$, those of the hole column $J(b-a)\,\Delta z$ (0.25 per
ampere-turn here), the core edges zero, and the edges leaving the core carry the correction of the gradient (0.5 on
the vertical edges above and below the core, -0.25 on the horizontal ones at its corners): the flux of the hole is
routed around the core. The core keeps a free $\mathbf{T}$ (its eddy currents) with the usual gauges. The winding cells are
`SOLID_CONDUCTOR` (their $\mathbf{T}$ is fixed: no eddy current, their $\rho$ unused), the rest of the hole
`CUT_LAYER`.

The response to the ramp: without the core B at the centre follows the current; with it, the eddy currents of the
core delay the field there (the core's diffusion time $\mu_0\sigma c^2$ is about 0.08 s) and decay after the ramp.
The flux of a turn, $\Psi=\int\mathbf{T}_0\cdot\mathbf{B}\,dV$ per ampere-turn, gives the induced voltage of an
N-turn coil, $N\,d\Psi/dt$: T-phi has no voltage drive, the voltage comes afterwards from the flux.
""")
    nb.code('''
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
for out, name in ((no_core, "no core"), (with_core, "conducting core")):
    ax[0].plot(out["t"], out["Bc"], "o-", ms=3, label=name)
ax[0].plot(with_core["t"], with_core["Bc"][-1] * with_core["NI"] / with_core["NI"][-1], "k:", label="current (scaled)")
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("B_z at the centre (T)"); ax[0].legend(); ax[0].grid(alpha=0.3)
ax[1].plot(with_core["t"], with_core["Jcore"] / 1e3, "o-", ms=3, color="C1")
ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("largest |J| in the core (kA/m^2)"); ax[1].grid(alpha=0.3)
plt.show()
for out, name in ((no_core, "no core"), (with_core, "conducting core")):
    v = 1000 * np.diff(out["Psi"]) / np.diff(out["t"])
    print(f"{name:16s}: Psi = {out['Psi'][-1]:.4f} Wb per turn at 500 kA-turns; N dPsi/dt during the ramp "
          f"(N = 1000) {v[2]:.0f} V, after it {v[6]:.1f} V")
''')
    nb.md("""
## Takeaways

* A conductor carrying a net current needs a cut: a one-cell `CUT_LAYER` layer, `_T` = 0 on its planar edges, the
  crossing edges (those on the conductor surface included) coupled (CP ids 1, 2, 3; induced), imposed (`_T` = I) or
  coupled and loaded by `VLTG` (voltage). Without it the solver silently models the conductor cut open.
* `tree_gauge.edges()` gives the edges of any cells in the solver's orientation; `tree_gauge.edge_circulations()`
  the edge values of a known field.
* A current-driven coil is `_T` = NI(t) x (edge values of T0) on the winding and the hole column; a conductor in the
  hole is cleared by a gradient; the induced voltage is N dPsi/dt afterwards.
* Several conductors: one cut each, one free cut DOF per hole; series and parallel connections: user manual, *Cuts*.
""")
    return nb


# ======================================================================================================================
def superconductors():
    nb = nbbuild.Notebook()
    nb.md(r"""
# Superconductors: power laws, AC loss and trapped fields (TPHI2D_edge)

This notebook follows `examples/wire_2D_tphi` and `examples/bulk_2D_tphi_axisym`: the 2D H-formulation (T-phi on 2D
cells, `TPHI2D_edge`) with the E-J power law of superconductors,

$$\mathbf{E}=E_c\left(\frac{|\mathbf{J}|}{J_c}\right)^n\frac{\mathbf{J}}{|\mathbf{J}|}\ (+\rho_0\mathbf{J}),$$

a resistivity $\rho(J)=E(J)/J$ given as a table. It shows:

1. how to write the E-J table (and why it ends with a normal-state segment);
2. the AC transport loss of a round wire (planar 2D, a cut) against Norris's critical-state formula, its current
   fronts, and the convergence to the critical state as n grows;
3. the nonlinear settings that work at REBCO-level n (and the ones that fail);
4. a REBCO bulk (axisymmetric) magnetized by field cooling with a field-dependent $J_c(B)$ (`power_law_b`): the
   trapped field and its creep.

See the user manual (*Element library: TPHI2D_edge*, *Nonlinear laws*: *Superconductors: the nonlinear settings*).
The wire alone runs about 40 s per exponent.
""")
    nb.code(SETUP.format(name="sc", cases='("wire_2D_tphi", "bulk_2D_tphi_axisym")'))
    nb.md(r"""
## The E-J table

`set_power_law(mat, table)` takes rows $[J, E]$ from $J=0$, both increasing; the law `power_law` interpolates it
linearly and **continues along its last segment**. Three rules (measured on the `formulations` benchmarks):

* enough points where the curve bends: log-spaced from $10^{-4}J_c$, 0.6 % apart;
* a negligible floor $\rho_0$ (its diffusion time $\mu_0L^2/\rho_0$ much longer than the run): it is the resistivity
  at $J\to0$ (the first segment);
* **a normal-state end**: above $E_\text{max}=100E_c$ a segment of the normal resistivity (current sharing). A pure
  power law continued by a Newton iterate far above $J_c$ gives tangents of $10^{30}$, and the factorization fails.
""")
    nb.code('''
def power_law_table(jc, n, ec=1e-4, rho0=1e-14, e_max=1e-2, rho_n=1e-6):
    j = np.concatenate(([0.0], jc * np.geomspace(1e-4, 10.0, 2001)))
    e = rho0 * j + ec * (j / jc) ** n
    keep = e < e_max
    j_max = jc * (e_max / ec) ** (1.0 / n)
    return np.column_stack((np.concatenate((j[keep], [j_max, 10 * j_max])),
                            np.concatenate((e[keep], [e_max, e_max + rho_n * 9 * j_max]))))

fig, ax = plt.subplots(figsize=(6.5, 4))
for n in (10, 30, 100):
    t = power_law_table(1e4, n)
    ax.loglog(t[1:, 0] / 1e4, t[1:, 1], label=f"n = {n}")
ax.axhline(1e-4, color="k", ls=":", lw=1); ax.text(1.2e-4, 1.5e-4, "Ec")
ax.set_xlim(1e-4, 10); ax.set_ylim(1e-20, 1e2)
ax.set_xlabel("J / Jc"); ax.set_ylabel("E (V/m)"); ax.grid(True, which="major", alpha=0.3); ax.legend()
ax.set_title("power law + floor + normal-state end"); plt.show()
''')
    nb.md(r"""
## A round wire with an AC transport current

The wire (radius R = 1, scaled units: $J_c=10^4$ A/m², $I_c=J_c\pi R^2$) carries $I(t)=0.5\,I_c\sin(2\pi t/T)$,
T = 40 s (E about $E_c$). Its net current needs a **cut**: a strip of `CUT_LAYER` cells from the wire to the box along
the ray y = 0+, x > R. `_T` = 0 on the strip's edges along the ray and on its other side; `_T` = $\pm I(t)$ on its
crossing edges (the sign from `tree_gauge.edges`); the box carries the exact potential outside a round wire,
$\phi=-I\theta/2\pi$, with $\theta=2\pi$ on the ray: $\phi$ jumps by I across the strip (T carries the jump there).

`wire(n, settings)` runs 1.25 cycles (80 steps per cycle) and returns the loss per cycle, from the first peak to the
next, $Q=\oint\int\mathbf{E}\cdot\mathbf{J}\,dA\,dt$ (the `E` and `Jt` outputs, J/m), the Newton iterations, and the
current density along the x axis at two instants.
""")
    nb.code('''
JC, EC, F, PERIOD, PER_CYCLE = 1e4, 1e-4, 0.5, 40.0, 80
IC = JC * np.pi
SETTINGS = {"absolute_residual_tolerance": 1e-9, "relative_residual_tolerance": 1e-30, "max_iterations": 100,
            "require_residual_reduction": True}

def norris(f, ic=IC):
    return MU0 * ic**2 / np.pi * ((1 - f) * np.log(1 - f) + (2 - f) * f / 2)

def wire(n, settings=SETTINGS, steps=PER_CYCLE + PER_CYCLE // 4):
    solver = new_solver("wire_2D_tphi", f"wire_n{n:g}.log")
    solver.solver_settings = {"petsc_options": QUIET, "nonlinear": dict(settings, verbose=False)}
    grid = solver.model.mesh
    xyz = meshes.points(grid)
    types = meshes.array_values(grid, "TYPE", "elem")
    solver.model.set_power_law(2, power_law_table(JC, n))
    select_surface(grid, 2, [1, 3], midside=True)
    selection.set_value_on_selection(grid, "_T", 0.0, 2)           # n x T = 0 on the wire surface
    select_inside(grid, 2, [1, 3])
    selection.set_value_on_selection(grid, "_PHI", 0.0)            # PHI = 0 inside the wire
    on_ray = np.isclose(xyz[:, 1], 0.0) & (xyz[:, 0] > 0.0)
    selection.select(grid, "TYPE", 3, location="elem")
    mids, start, end = tree_gauge.edges(grid, selection.selected_ids(grid, "elem"))
    cross = on_ray[start] != on_ray[end]
    sign = np.where(on_ray[start[cross]], 1.0, -1.0)
    selection.select(grid, "ID", mids[~cross])
    selection.set_value_on_selection(grid, "_T", 0.0, 2)           # the strip edges along the cut
    selection.select_all(grid); selection.nsel_ext(grid)
    box = selection.selected_ids(grid, "node")
    theta = np.mod(np.arctan2(xyz[box, 1], xyz[box, 0]), 2 * np.pi)
    theta[on_ray[box]] = 2 * np.pi
    with quietly():
        solver.load_model_into_solver()
    solver.time_step["theta"] = 1.0
    c, vol = meshes.cell_centroids(grid), meshes.cell_sizes(grid)
    axis = (types == 2) & (c[:, 1] > 0) & (c[:, 1] < 0.07 * np.abs(c[:, 0]) + 0.03)   # the first row above y = 0
    order = np.argsort(c[axis, 0])
    dt = PERIOD / PER_CYCLE
    out = {"loss": [], "newton": [], "profiles": {}}
    try:
        for step in range(1, steps + 1):
            t = step * dt
            i = F * IC * np.sin(2 * np.pi * t / PERIOD)
            selection.select(grid, "ID", mids[cross])
            selection.set_value_on_selection(grid, "_T", sign * i, 2)
            selection.select(grid, "ID", box)
            selection.set_value_on_selection(grid, "_PHI", -i * theta / (2 * np.pi))
            solver.time_step["time"] = t
            with quietly():
                solver.solve()
            e, j = meshes.array_values(grid, "E", "elem"), meshes.array_values(grid, "Jt", "elem")
            out["loss"].append(float(np.einsum("ij,ij->i", e, j) @ vol) * dt)
            out["newton"].append(solver.steps[-1]["newton_iterations"])
            if step in (PER_CYCLE // 4, PER_CYCLE // 2):
                out["profiles"][step] = (c[axis, 0][order], j[axis, 2][order] / JC)
    finally:
        with quietly():
            solver.close()
    q = sum(out["loss"][PER_CYCLE // 4:PER_CYCLE // 4 + PER_CYCLE])
    out.update(Q=q, ratio=q / norris(F))
    return out

w30 = wire(30)
print(f"n = 30: Q = {w30['Q']:.3f} J/m per cycle, Norris {norris(F):.3f} J/m, ratio {w30['ratio']:.3f}; "
      f"Newton iterations per step: mean {np.mean(w30['newton']):.1f}, max {max(w30['newton'])}")
''')
    nb.md(r"""
The current fronts: at the first peak (t = T/4) the current flows in a shell at the surface, $J\approx J_c$ for
$R\sqrt{1-F}<r<R$ in the critical state (dotted); at the zero crossing (t = T/2) the current of the previous half
cycle stays in the inner part of the shell and a reversed layer has entered from the surface.
""")
    nb.code('''
fig, ax = plt.subplots(figsize=(7, 3.8))
for step, label in ((PER_CYCLE // 4, "t = T/4 (peak)"), (PER_CYCLE // 2, "t = T/2 (I = 0)")):
    x, j = w30["profiles"][step]
    ax.plot(x, j, "o-", ms=3, label=label)
r_front = np.sqrt(1 - F)
ax.plot([-1, -r_front, -r_front, r_front, r_front, 1], [1, 1, 0, 0, 1, 1], "k:", label="Bean, at the peak")
ax.set_xlabel("x / R (the cells just above y = 0)"); ax.set_ylabel("J_z / Jc"); ax.grid(alpha=0.3); ax.legend(); plt.show()
''')
    nb.md(r"""
## The exponent: from a resistor to the critical state

The power law tends to Bean's critical state as n grows. The ratio of the loss to Norris's formula goes from well
above 1 at low n towards 1; at n = 30-100 (REBCO, Nb-Ti) it is 1.2-1.1: compare with the critical state at the n of
the material, not with $n\to\infty$. The Newton iterations grow with n.
""")
    nb.code('''
sweep = {30: w30}
for n in (10, 100):
    sweep[n] = wire(n)
print(" n    Q / Q_Norris   Newton its (mean, max)")
for n in sorted(sweep):
    print(f"{n:3d}   {sweep[n]['ratio']:10.3f}     {np.mean(sweep[n]['newton']):5.1f}, {max(sweep[n]['newton'])}")
''')
    nb.md(r"""
## Nonlinear settings: what fails, and why

The defaults of `solver_settings["nonlinear"]` (relative residual tolerance $10^{-8}$, 50 iterations, no residual
reduction required) are made for iron. On the wire at n = 30:
""")
    nb.code('''
for name, settings in (("defaults", {}),
                       ("absolute tolerance, no residual reduction",
                        {"absolute_residual_tolerance": 1e-9, "relative_residual_tolerance": 1e-30, "max_iterations": 100})):
    try:
        wire(30, settings, steps=6)
        print(f"{name}: 6 steps solved")
    except pytec3d.errors.ConvergenceError as error:
        print(f"{name}: {str(error).split('. See')[0]}")
''')
    nb.md(r"""
Both stop at the third step: the transport current changes at the few cells next to the cut, the first iterate has
$J\gg J_c$ there, and the full Newton steps run away. The rules (user manual, *Superconductors: the nonlinear
settings*):

* an **absolute** residual tolerance scaled to the problem (a fraction of $E_c$ times a face area per T row; here
  $10^{-9}$), the relative one off ($10^{-30}$): a relative tolerance is relative to the first residual of the step,
  huge when the starting guess has $J>J_c$ (at n = 100 it accepted $E\approx10^{20}E_c$ on the slab benchmark); a loose
  absolute tolerance stops at once on a wrong state; check $\max|E|$ against $E_c$;
* **require the residual to decrease** for a transport current through a cut (as here), **but not** from a
  current-free conductor in a changing field (the bulk below: the first full step overshoots and must be taken);
* Brent and plain backtracking give the same results; Picard does not converge on the power law.

## A field-cooled REBCO bulk with Jc(B)

The bulk disk (radius 10 mm, thickness 10 mm, axisymmetric half section, `option3` = 1: (x, y) = (r, z), component 2
of `Jt` and `E` azimuthal) is cooled in 2 T, then the field is ramped to zero in 100 s and held. The critical current
depends on the field, the Kim model $J_c(B)=J_{c0}/(1+|B|/B_0)$, $J_{c0}=3\times10^8$ A/m², $B_0$ = 0.5 T, n = 25:
the law `power_law_b`, $E=E_0(J/s(|B|))$, two tables (`set_power_law(mat, table, jc_b=rows)`).

* no cut: the section touches the axis, the air of the half plane is simply connected;
* the **static step** first: the cooled state, B uniform and no current. In a static step the model sets T = 0 on
  every conductor edge and frees $\phi$ inside; the transient steps restore $\mathbf{n}\times\mathbf{T}=0$ on the
  surface and $\phi$ = 0 inside;
* the applied field on the box, $\phi=-H_a(t)z$.

We run the example macro as it is (`examples/bulk_2D_tphi_axisym/macro.py`) and read its output files.
""")
    nb.code('''
import csv, runpy
shutil.copy("bulk_2D_tphi_axisym_model_data.json", "model_data.json")    # the macro reads model_data.json
with quietly(), contextlib.redirect_stderr(io.StringIO()):
    namespace = runpy.run_path(str(EXAMPLES / "bulk_2D_tphi_axisym" / "macro.py"), run_name="__main__")
    namespace["pytec3d"].close()                                         # writes the steps file
steps = list(csv.DictReader(open("output/bulk_2D_tphi_axisym_steps.csv")))
print(len(steps), "steps; columns of the steps file:", ", ".join(steps[0]))
print("Newton iterations per step:", " ".join(s["newton_iterations"] for s in steps))
''')
    nb.md(r"""
The steps file (`output/<model>_steps.csv`) lists each step with its time, Newton iterations, convergence and output
file; from the output files we read B at the probe (1 mm above the centre of the top face) and the current density.
""")
    nb.code('''
probe_xyz = np.array([0.0, 6e-3, 0.0])
t, b_probe = [], []
for s in steps:
    g = pv.wrap(meshes.read_mesh(s["output"]))
    k = int(np.argmin(np.linalg.norm(g.points - probe_xyz, axis=1)))
    t.append(float(s["time"])); b_probe.append(g.point_data["B"][k, 1])
t, b_probe = np.array(t), np.array(b_probe)
print(f"trapped field at the probe: {b_probe[np.argmin(np.abs(t - 101.0))]:.3f} T at the end of the ramp, "
      f"{b_probe[-1]:.3f} T at t = {t[-1]:.0f} s")

ba = 2.0 * np.clip(1 - (t - 1.0) / 100.0, 0, 1)
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
ax[0].plot(t, ba, "k--", label="applied B")
ax[0].plot(t, b_probe, "o-", ms=3, label="B_z, probe 1 mm above the bulk")
ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("T"); ax[0].legend(); ax[0].grid(alpha=0.3)
hold = t > 101.5
ax[1].semilogx(t[hold] - 101.0, b_probe[hold], "o-", ms=3)
ax[1].set_xlabel("time after the ramp (s)"); ax[1].set_ylabel("trapped B_z at the probe (T)")
ax[1].set_title("flux creep"); ax[1].grid(True, which="both", alpha=0.3)
plt.show()

g = pv.wrap(meshes.read_mesh(steps[-1]["output"]))
bulk = g.threshold((2, 2), scalars="TYPE")
bulk.cell_data["J_theta / Jc0"] = bulk.cell_data["Jt"][:, 2] / 3e8
p = pv.Plotter(window_size=(700, 420))
p.add_mesh(bulk, scalars="J_theta / Jc0", cmap="coolwarm", clim=(-0.8, 0.8))
p.view_xy(); p.add_text("trapped current, r-z half section (end of the hold)", font_size=10); p.show()
''')
    nb.md(r"""
The bulk traps about 0.5 T at the probe; during the hold the currents relax (with a power law E > 0 wherever
$J\approx J_c$) and the trapped field decays roughly linearly in log t, as measured on real bulks. The current fills
the whole section, lowest where |B| is largest (the Kim model: the centre of the disk).

## Takeaways

* `TPHI2D_edge`: planar (per `thickness`) or axisymmetric (`option3` = 1); J and E out of the plane (component 2).
* Power laws as tables with a floor and a normal-state end; `power_law_b` for $J_c(B)$ (the Jacobian has no
  $\partial/\partial B$ term: more iterations, the same solution).
* High n needs an absolute residual tolerance; `require_residual_reduction` for transport currents through cuts, not
  for current-free starts in a changing field.
* A static step gives an initial magnetostatic state (field cooling); the transient steps restore the gauges.
* The power-law loss approaches the critical state as n grows (1.2 at n = 30, 1.1 at n = 100).
""")
    return nb


if __name__ == "__main__":
    target = os.path.abspath(sys.argv[1])     # the notebooks change the working directory
    which = sys.argv[2:] or ["scalar_potential", "tphi_cuts", "superconductors"]
    for name in which:
        nbbuild.run(globals()[name](), f"{target}/tutorial_{name}.ipynb")

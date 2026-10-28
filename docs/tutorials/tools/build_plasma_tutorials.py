"""Build the plasma tutorial notebooks (ITER equilibrium and evolution) from the test cases, with their outputs.

Usage, from docs/tutorials/tools, with the pytec3d environment (pyvista installed):
    python build_plasma_tutorials.py ..                    # both
    python build_plasma_tutorials.py .. iter_equilibrium   # one of them (iter_3d_halo: about 20 minutes)
"""
import os
import sys

import nbbuild
from build_tutorials import SETUP

# Plot helpers shared by the two notebooks: cell polygons and triangles of a 2D mesh for matplotlib.
PLOT_HELPERS = '''
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation

def polygons(grid, mask=None):
    """The (R, Z) polygons of the cells of a 2D grid (all of them, or those where mask is True)."""
    g = pv.wrap(grid)
    points, conn, offset = g.points[:, :2], g.cell_connectivity, g.cell_offsets
    ids = np.arange(g.n_cells) if mask is None else np.flatnonzero(mask)
    return [points[conn[offset[i]:offset[i + 1]]] for i in ids]

def fill(ax, grid, mask, values=None, **kwargs):
    """Draw cells: one colour (facecolor=...) or coloured by the cell values (cmap=..., clim=...)."""
    collection = PolyCollection(polygons(grid, mask), **kwargs)
    if values is not None:
        collection.set_array(np.asarray(values)[mask])
    ax.add_collection(collection)
    return collection

def triangulation(grid):
    """The grid as matplotlib triangles over its points (for contours of point arrays such as PSI)."""
    t = pv.wrap(grid).triangulate()
    return Triangulation(t.points[:, 0], t.points[:, 1], t.cells.reshape(-1, 4)[:, 1:])

def frame(ax, r=(0.0, 13.0), z=(-9.0, 9.0)):
    ax.set_xlim(*r); ax.set_ylim(*z); ax.set_aspect("equal"); ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
'''


# ======================================================================================================================
def iter_equilibrium():
    nb = nbbuild.Notebook()
    nb.md(r"""
# ITER free-boundary plasma equilibrium (2D axisymmetric)

This notebook computes the free-boundary equilibrium of the ITER 15 MA reference plasma, following
`examples/ITER_2D_fine2` (the refined linear mesh), and looks inside every step of the computation. It is the starting point of every plasma
evolution (see the next tutorial).

**The problem.** In an axisymmetric plasma the poloidal flux per radian $\psi(R,Z) = R A_\varphi$ obeys the
Grad-Shafranov equation

$$
\Delta^*\psi \;=\; R\frac{\partial}{\partial R}\Big(\frac{1}{R}\frac{\partial\psi}{\partial R}\Big)
  + \frac{\partial^2\psi}{\partial Z^2} \;=\; -\mu_0 R\,J_\varphi,
\qquad
J_\varphi = R\,p'(\psi) + \frac{F F'(\psi)}{\mu_0 R}\;\text{in the plasma},
$$

with $J_\varphi$ the imposed current density in the coils and zero in vacuum. It is *free boundary*: the plasma
is the region inside the last closed flux surface, which is itself part of the solution. pyTEC3D does not solve
the Grad-Shafranov equation as such: the field equation is the finite-element magnetostatic problem of the
element library ($\nabla\times\nu\nabla\times A = J$, which in axisymmetry is the same operator), and the
`Algorithm` class wraps it in a **Picard iteration** on the plasma current density:

1. from $\psi$, find the magnetic axis (maximum of $\psi$), the X point (saddle) and the limiter point (the wall
   point of largest $\psi$); the larger of $\psi_X$, $\psi_{lim}$ is the boundary flux $\psi_b$, and the cells
   inside it connected to the axis are the plasma;
2. set $J_\varphi = \lambda\,(1-\bar\psi^{\alpha_m})^{\alpha_n}\big(\beta\,R/R_0 + (1-\beta)\,R_0/R\big)$ in the
   plasma, with $\bar\psi = (\psi-\psi_{axis})/(\psi_b-\psi_{axis})$ and $\lambda$ scaling the total to $I_p$;
3. correct the currents of the control coils so that the poloidal field at the target axis vanishes (position
   control);
4. solve the field with the new currents; stop when $\psi$ no longer changes.

What we do: look at the machine model, set the coil currents and the algorithm parameters (and check where the
control parameters come from), run the iterations one by one, then analyse the converged equilibrium: flux map,
separatrix, current density, pressure and toroidal-field profiles, global quantities, and a 3D view. The run takes
about 10 s. Theory: *Theory - Plasma equilibrium and evolution*; options: *User manual - Plasma algorithm*.
""")
    nb.code(SETUP.format(case="ITER_2D_fine2", files='("ITER_2D_fine2.vtu", "model_data.json")')
            + "\nimport json, time\nfrom pytec3d import Algorithm, ConvergenceError, element_library\n" + PLOT_HELPERS)
    nb.md("""
## The machine

`ITER_2D_fine2.vtu` is the poloidal cross-section (x = R, y = Z) of the whole machine in a 100 m x 200 m half-plane of
air: the six central-solenoid modules CS3U ... CS3L, the six poloidal-field coils PF1 ... PF6, the plasma region
(where the plasma can be) and the **passive structures**: the double-wall vacuum vessel and the in-vessel supports,
massive steel conductors. Every component is a cell array with value 1 on its cells. The model data gives two TYPE
ids: 21 air (`NON_CONDUCTING`) and 23 for all the rest (`STRANDED_CONDUCTOR`, a region where the current density `Js`
is imposed). The equilibrium is a sequence of static solves, so no eddy current flows and the conductivities play no
role. In the evolution (next tutorials) the passive structures become **solid conductors** (`SOLID_CONDUCTOR`): the
massive conductors in which the plasma motion induces eddy currents.
""")
    nb.code('''
solver = Solver(log_file="ITER_2D_fine2.log")
algorithm = Algorithm(solver)
solver.read_model("ITER_2D_fine2.vtu")
solver.model.read_json("model_data.json")
mesh = solver.model.mesh
grid = pv.wrap(mesh)

components = ["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L", "PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
              "PLASMA_REGION", "PASSIVE_STRUCTURES"]
area = grid.cell_data["cell_volume"]            # in 2D the cell size is its area [m^2]
centroid = grid.cell_data["cell_centroids"]
print(f"{grid.n_points} points, {grid.n_cells} cells; TYPE ids", np.unique(grid.cell_data["TYPE"]))
print(f"{'component':>18} {'cells':>6} {'area [m2]':>10} {'R [m]':>6} {'Z [m]':>6}")
for name in components:
    on = grid.cell_data[name] == 1
    r, z = (centroid[on, :2] * area[on, None]).sum(axis=0) / area[on].sum()
    print(f"{name:>18} {on.sum():6d} {area[on].sum():10.4f} {r:6.2f} {z:6.2f}")
''')
    nb.code('''
fig, ax = plt.subplots(figsize=(6, 7.5))
fill(ax, mesh, grid.cell_data["PLASMA_REGION"] == 1, facecolor="#dbe9f6", edgecolor="none")
fill(ax, mesh, grid.cell_data["PASSIVE_STRUCTURES"] == 1, facecolor="#2c6e9e", edgecolor="#2c6e9e", linewidth=0.6)
for name in components[:6]:
    fill(ax, mesh, grid.cell_data[name] == 1, facecolor="#c05a2a", edgecolor="none")
for name in components[6:12]:
    fill(ax, mesh, grid.cell_data[name] == 1, facecolor="#e3a33b", edgecolor="none")
for name in components[:12]:
    on = grid.cell_data[name] == 1
    r, z = (centroid[on, :2] * area[on, None]).sum(axis=0) / area[on].sum()
    ax.annotate(name.replace("PFC_", "PF"), (r, z), xytext=(4, 4) if r > 2 else (-4, 0),
                textcoords="offset points", fontsize=8, ha="left" if r > 2 else "right")
ax.plot(6.52136102486433, 0.536537145137914, "k+", ms=10, label="target magnetic axis")
frame(ax, r=(0, 13.5)); ax.legend(loc="lower right", fontsize=8)
ax.plot([], [], color="#2c6e9e", lw=4, label="passive structures (solid conductors)")
ax.set_title("CS modules, PF coils, plasma region (light blue), passive structures (dark blue)", fontsize=9)
plt.show()
''')
    nb.md("""
The model data rows (`model_data.json`): `option1` is the domain of each TYPE, `option2 = 1` the vector potential
`A_z` as unknown, `option3 = 1` axisymmetry. Rows 24 and 25 (solid conductor, stranded coil) are not used by any
cell here: they serve the evolution cases, which share the same file layout.
""")
    nb.code('''
data = json.load(open("model_data.json"))
domains = dict(enumerate(element_library.domains()))
for tid, row in data["element_types"].items():
    print(f"TYPE {tid}: {domains.get(row['option1'], row['option1']):26s} {row}")
''')
    nb.md("""
## Boundary condition and analysis type

$A_z = 0$ on the external boundary: on the axis R = 0 (where $\\psi = R A_\\varphi = 0$ anyway) and far away
(the field of a current loop decays fast, so 100 m is far enough). Every Picard iteration is a magnetostatic solve.
The model is loaded into a solver session once: the arrays (`Js` in particular) can be changed between solves, and
the factorized matrix is reused at every iteration, since the step does not change.
""")
    nb.code('''
selection.select_all(mesh)
selection.nsel_ext(mesh)
selection.set_value_on_selection(mesh, "_A", 0.0, 2)

solver.time_step.update(time=1.0, transient_effects=False, theta=1.0)
with quietly():
    solver.load_model_into_solver()
''')
    nb.md(r"""
## Position control: computing the `b*_value` coefficients

At every iteration the algorithm measures the poloidal field $(B_R, B_Z)$ at the target position of the magnetic
axis and computes a *radial* and a *vertical* correction current, $I_{rad}$ and $I_{ver}$, that would cancel it:

$$
\begin{pmatrix} b_{rad,x} & b_{ver,x}\\ b_{rad,y} & b_{ver,y}\end{pmatrix}
\begin{pmatrix} I_{rad}\\ I_{ver}\end{pmatrix} = -\begin{pmatrix} B_R\\ B_Z\end{pmatrix}.
$$

The four coefficients are the field at the target of a **unit correction current**: a radial correction of 1 A is
+1 A in PF2 and -1 A in PF5 (the two outer coils, which push the plasma in and out), a vertical one +1 A in PF3 and in
PF4 (which push it up and down). They are properties of the coils and of the mesh only, so we compute them before the
equilibrium, with two magnetostatic solves of the model: 1 MA of radial (then vertical) correction current in the
coils and no other current, and the flux density read at the target point (the nodal field ``B``, interpolated at
the point).
""")
    nb.code('''
target = (6.52136102486433, 0.536537145137914, 0.0)
coil_area = {"PFC_1": 0.94375190, "PFC_2": 0.41453946, "PFC_3": 0.66413094, "PFC_4": 0.60871516,
             "PFC_5": 0.77496250, "PFC_6": 1.72659250, "CS": 1.55608}
unit_corrections = {"radial": {"PFC_2": +1.0, "PFC_5": -1.0}, "vertical": {"PFC_3": +1.0, "PFC_4": +1.0}}
probe = pv.PolyData(np.array([target]))
js = meshes.array_values(mesh, "Js", "elem")
b_unit = {}
for name, coils in unit_corrections.items():
    js[:] = 0.0
    for coil, sign in coils.items():
        js[grid.cell_data[coil] == 1, 2] = sign * 1e6 / coil_area[coil]      # 1 MA through the coil cross-section
    meshes.get_array(mesh, "Js", "elem").Modified()
    with quietly():
        solver.solve()
    sampled = probe.sample(pv.wrap(mesh))                                    # interpolate the nodal B at the target
    b_unit[name] = sampled.point_data["B"][0, :2] / 1e6                       # T per A of correction current

macro = {"radial": (-3.873415861e-08, 5.313857037e-09), "vertical": (1.822467692e-10, 1.133287835e-07)}
print(f"{'':10s} {'computed':>24s} {'test macro':>24s}")
for name, label in (("radial", "rad"), ("vertical", "ver")):
    for k, axis in enumerate("xy"):
        print(f"b{label}{axis}_value  {b_unit[name][k]:+24.4e} {macro[name][k]:+24.4e}")
''')
    nb.md("""
The computed coefficients agree with the values written in the test macro, which were obtained the same way. In the
test the control coils have the gains +-0.5 (`control_coils` below), while the coefficients are those of unit gains:
the algorithm computes the correction current that would cancel the field, and the coils receive half of it. The
gains 0.5 are an **under-relaxation of 0.5** of the position control, which damps the coupling between the moving
plasma and its correction from one iteration to the next.

Note also that the correction is recomputed at every iteration from the field of the current state and *replaces*
the previous one (the coils keep their scenario current as the base): the control is proportional, so a residual
field, and a small offset of the axis from the target, remain at convergence. We will measure them below.

## Coil currents and the initial guess

The coil currents of the scenario point (ampere-turns of each coil) become current densities `Js_z = I / a` on
the cells of each coil. The plasma current is not known yet: the first solve uses 15 MA spread uniformly over the
whole plasma region, a rough guess that only has to give a flux with an axis and an X point inside the chamber.
""")
    nb.code('''
coil_current = {"CS3U": -4.8204e06, "CS2U": -4.0480e06, "CS1U": -2.0804e07, "CS1L": -2.0804e07,
                "CS2L": -1.0297e07, "CS3L": 5.3012e06,
                "PFC_1": 5.7325e06, "PFC_2": -2.8834e06, "PFC_3": -5.8084e06, "PFC_4": -4.7762e06,
                "PFC_5": -7.8133e06, "PFC_6": 1.6941e07}
plasma_current, chamber_area = 15e6, 28.6027095

selection.select_all(mesh)
selection.set_value_on_selection(mesh, "Js", 0.0, location="elem")       # forget the response solves above
for coil, current in coil_current.items():
    selection.select(mesh, coil, 1, location="elem")
    selection.set_value_on_selection(mesh, "Js", current / coil_area["CS" if coil.startswith("CS") else coil], 2,
                                     location="elem")
selection.select(mesh, "PLASMA_REGION", 1, location="elem")
selection.set_value_on_selection(mesh, "Js", plasma_current / chamber_area, 2, location="elem")

js = meshes.array_values(mesh, "Js", "elem")
for name in ("COILS", "PLASMA_REGION"):
    on = grid.cell_data[name] == 1
    print(f"total current in {name}: {(js[on, 2] * area[on]).sum() / 1e6:+.2f} MA(-turns)")
''')
    nb.md(r"""
## Algorithm parameters

| parameter | value | meaning |
|---|---|---|
| `target_central_axis_x/y` | 6.521, 0.537 m | where the position control wants the magnetic axis |
| `R0` | 6.2 m | major radius of the current profile |
| `plasma_current` | 15 MA | $I_p$, imposed through $\lambda$ |
| `beta_p` | 0.682 | $\beta$ of the profile: the share of the pressure term $\beta R/R_0$ |
| `alpha_mp`, `alpha_np` | 1.0001, 0.83197 | shape of the profile $(1-\bar\psi^{\alpha_m})^{\alpha_n}$ |
| `AF0` | $5.3\times6.2$ T m | vacuum $F = R B_\varphi$: 5.3 T at 6.2 m (used for $F$ and the poloidal currents) |
| `control_coils`, `b*_value` | gains +-0.5, coefficients computed above | position control |
| `tolerance_value` | 1e-4 Wb/rad | convergence: $\max\lvert\psi_k-\psi_{k-1}\rvert$ |
| `max_picard_iterations` | 70 | the solve raises `ConvergenceError` beyond |
""")
    nb.code('''
algorithm.input_parameters.update(
    target_central_axis_x=target[0], target_central_axis_y=target[1], target_central_axis_z=0.0,
    R0=6.2, plasma_current=plasma_current, beta_p=0.68216, alpha_mp=1.0001, alpha_np=0.83197, AF0=5.3 * 6.2,
    control_coils=[{"component": "PFC_2", "radial_gain": 0.5, "vertical_gain": 0.0},
                   {"component": "PFC_3", "radial_gain": 0.0, "vertical_gain": 0.5},
                   {"component": "PFC_4", "radial_gain": 0.0, "vertical_gain": 0.5},
                   {"component": "PFC_5", "radial_gain": -0.5, "vertical_gain": 0.0}],
    bradx_value=b_unit["radial"][0], brady_value=b_unit["radial"][1],          # computed above
    bverx_value=b_unit["vertical"][0], bvery_value=b_unit["vertical"][1],
    tolerance_value=1e-4, max_picard_iterations=70)
''')
    nb.md("""
## The Picard iterations, one at a time

`algorithm.equilibrium_solve()` runs the whole loop and returns when $\\psi$ has converged, as the test macro does.
To watch the iterations we call it with `max_picard_iterations = 1`: each call does one iteration and raises
`ConvergenceError` until the last one; the state of the algorithm (previous flux, control-coil currents) is kept
between calls, so the sequence is the same as one call with 70 iterations. After each iteration we keep the flux of
the plasma region (`algorithm.plasma_region_mesh`, the mesh on which the algorithm works) and the status.
""")
    nb.code('''
algorithm.input_parameters["max_picard_iterations"] = 1
history, started = [], time.perf_counter()
for iteration in range(1, 71):
    try:
        with quietly():
            converged = algorithm.equilibrium_solve()
    except ConvergenceError:
        converged = False
    chamber = algorithm.plasma_region_mesh
    history.append({"psi": meshes.array_values(chamber, "PSI", "node").copy(), **algorithm.status})
    if converged:
        break
status = algorithm.status
print(f"converged at iteration {iteration} in {time.perf_counter() - started:.1f} s")
offset = np.subtract(status["magnetic_axis"][:2], target[:2])
print(f"axis - target = ({offset[0] * 100:+.1f}, {offset[1] * 100:+.1f}) cm, residual field at the target "
      f"(B_R, B_Z) = ({status['b_target_axis'][0] * 1e3:+.2f}, {status['b_target_axis'][1] * 1e3:+.2f}) mT")
for key in ("separatrix_configuration", "magnetic_axis", "x_point", "limiter_point", "psi_magnetic_axis",
            "psi_min", "plasma_element_count", "plasma_current_scaling_lambda",
            "radial_correction_current", "vertical_correction_current", "flux_norm"):
    value = status[key]
    print(f"  {key:32s}", tuple(round(v, 4) for v in value[:2]) if isinstance(value, tuple) else value)
''')
    nb.md("""
The status after every iteration is also appended to `output/ITER_2D_fine2_iteration_diagnostics.csv` in the working
folder (one row per iteration, the points split in x/y/z columns), which is what you would look at after a batch
run:
""")
    nb.code('''
diag = np.genfromtxt("output/ITER_2D_fine2_iteration_diagnostics.csv", delimiter=",", names=True, dtype=None,
                     encoding="utf-8")
it = np.arange(1, len(diag) + 1)
fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
ax = axes[0, 0]
ax.semilogy(it[1:], diag["flux_metric"][1:], "o-", ms=3); ax.axhline(1e-4, color="k", ls="--", lw=0.8)
ax.set_ylabel("max |PSI_k - PSI_k-1| [Wb/rad]"); ax.set_title("convergence", fontsize=10)
ax = axes[0, 1]
ax.plot(it[1:], diag["magnetic_axis_x"][1:] - target[0], "o-", ms=3, label="R_axis - R_target")
ax.plot(it[1:], diag["magnetic_axis_y"][1:] - target[1], "s-", ms=3, label="Z_axis - Z_target")
ax.set_ylabel("[m]"); ax.set_ylim(-0.5, 0.2); ax.legend(fontsize=8); ax.set_title("magnetic axis vs target", fontsize=10)
ax = axes[1, 0]
ax.plot(it, diag["radial_correction_current"] / 1e3, "o-", ms=3, label="radial")
ax.plot(it, diag["vertical_correction_current"] / 1e3, "s-", ms=3, label="vertical")
ax.set_yscale("symlog", linthresh=10); ax.set_ylabel("correction current [kA]"); ax.legend(fontsize=8); ax.set_title("position control", fontsize=10)
ax = axes[1, 1]
ax.plot(it, diag["plasma_element_count"], "o-", ms=3, color="C2"); ax.set_ylabel("plasma cells", color="C2")
ax2 = ax.twinx(); ax2.plot(it, diag["plasma_current_scaling_lambda"] / 1e6, "s-", ms=3, color="C3")
ax2.set_ylabel("lambda [MA/m^2]", color="C3"); ax.set_title("plasma region and profile scaling", fontsize=10)
for ax in axes[1]:
    ax.set_xlabel("Picard iteration")
for ax in axes.flat:
    ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()
''')
    nb.md("""
How to read it:

* the first iteration only solves the initial guess (no plasma region yet: zero cells, no correction); from the
  second on, the plasma takes its shape: too small at first, then too large and touching the wall (limited, in
  some of the early iterations, while the control pulls it back), then diverted again for good (picture below);
* the flux change falls roughly geometrically: the Picard iteration is a fixed-point iteration, linearly
  convergent; its rate is set by how strongly the plasma shape reacts to its own current;
* the correction currents swing by several MA in the first iterations, while the plasma forms, then settle to a
  steady correction of a few tens of kA: being proportional to the residual field, they do not go to zero, and the
  axis stops a few cm from the target (the magnetic axis is also a mesh node: it moves in steps of the node
  spacing);
* $\\lambda$ settles when the plasma region stops changing.

If an equilibrium does not converge, these curves are the first thing to look at: an axis that oscillates around
the target asks for smaller gains (more relaxation); a large offset from the target, for larger gains; a plasma region that keeps changing size, for a better initial
guess or coil currents closer to the scenario.

### The plasma boundary during the iterations
""")
    nb.code('''
chamber = algorithm.plasma_region_mesh
chamber_tri = triangulation(chamber)
fig, ax = plt.subplots(figsize=(5.5, 7))
fill(ax, mesh, grid.cell_data["PASSIVE_STRUCTURES"] == 1, facecolor="0.8", edgecolor="none")
shown = [2, 4, 6, 8, len(history)]
colors = plt.cm.viridis(np.linspace(0, 0.9, len(shown)))
for k, color in zip(shown, colors):
    h = history[k - 1]            # the status of iteration k: boundary found on the flux of iteration k - 1
    ax.tricontour(chamber_tri, history[k - 2]["psi"], levels=[h["psi_min"]], colors=[color], linewidths=1.2, linestyles="solid")
    ax.plot(*h["magnetic_axis"][:2], "o", color=color, ms=4)
    ax.plot([], [], color=color, label=f"iteration {k}: {h['separatrix_configuration']}, {h['plasma_element_count']} cells")
ax.plot(*target[:2], "k+", ms=10, label="target axis")
frame(ax, r=(3.5, 9.5), z=(-5.5, 5.5)); ax.legend(fontsize=8, loc="lower right")
ax.set_title("boundary flux surface psi = psi_b and axis", fontsize=10)
plt.show()
''')
    nb.md("""
## The converged equilibrium

### Flux surfaces, separatrix and plasma region

The flux of the whole model is in the point array `PSI` of the solver mesh ($\\psi = R A_z$, Wb/rad). The plasma
(`PLASMA = 1` cells of the plasma region) is bounded by the separatrix $\\psi = \\psi_b = \\psi_X$ since the
configuration is diverted (`X-point`): the X point is below the plasma, and its two legs reach the divertor. The
limiter point, the wall point of largest flux, has $\\psi_{lim} < \\psi_X$: the plasma does not touch the wall.
""")
    nb.code('''
psi = meshes.array_values(mesh, "PSI", "node")
psi_axis, psi_b = status["psi_magnetic_axis"], status["psi_min"]
plasma = meshes.array_values(chamber, "PLASMA", "elem") == 1
jphi = meshes.array_values(chamber, "Js", "elem")[:, 2]

fig, ax = plt.subplots(figsize=(6.5, 8))
fill(ax, mesh, grid.cell_data["PASSIVE_STRUCTURES"] == 1, facecolor="0.75", edgecolor="none")
for name in components[:12]:
    fill(ax, mesh, grid.cell_data[name] == 1, facecolor="#e3a33b", edgecolor="none")
pc = fill(ax, chamber, plasma, jphi / 1e6, cmap="inferno", edgecolor="none")
fig.colorbar(pc, ax=ax, shrink=0.6, label="J_phi [MA/m^2]")
levels = psi_b + (psi_axis - psi_b) * np.linspace(-1.5, 1.0, 26)
ax.tricontour(triangulation(mesh), psi, levels=np.sort(levels), colors="steelblue", linewidths=0.5,
              linestyles="solid")
ax.tricontour(triangulation(mesh), psi, levels=[psi_b], colors="cyan", linewidths=1.5, linestyles="solid")
ax.plot(*status["magnetic_axis"][:2], "wo", ms=5, mec="k", label="magnetic axis")
ax.plot(*status["x_point"][:2], "rx", ms=9, mew=2, label="X point")
ax.plot(*status["limiter_point"][:2], "gs", ms=6, label="limiter point (not touching)")
ax.plot(*target[:2], "k+", ms=10, label="target axis")
frame(ax, r=(2.5, 10.5), z=(-6.5, 6.5)); ax.legend(fontsize=8, loc="lower left")
ax.set_title(f"psi contours, separatrix (cyan), J_phi in the plasma", fontsize=10)
plt.show()
print(f"psi_axis = {psi_axis:.4f}, psi_b = psi_X = {psi_b:.5f}, psi_limiter = {status['last_psi_limiter']:.4f} Wb/rad")
''')
    nb.md(r"""
### Current density across the plasma

On the midplane through the axis, $J_\varphi$ is the product of the flux profile $(1-\bar\psi^{\alpha_m})^{\alpha_n}$,
which vanishes at the boundary, and of $\beta R/R_0 + (1-\beta)R_0/R$, which shifts the peak outwards for large
$\beta$ (the pressure term grows with R). The same profile also gives the poloidal current density
$\mathbf J_p = F'(\psi)\,\mathbf B_p$, in `Js_x`, `Js_y` (the current that flows along the flux surfaces in the
poloidal plane, linked to the toroidal field change $F(\psi)$ inside the plasma).
""")
    nb.code('''
cc = meshes.array_values(chamber, meshes.CENTROIDS_ARRAY, "elem")
js_chamber = meshes.array_values(chamber, "Js", "elem")
band = plasma & (np.abs(cc[:, 1] - status["magnetic_axis"][1]) < 0.12)
order = np.argsort(cc[band, 0])
fig, ax = plt.subplots(figsize=(7, 3.4))
ax.plot(cc[band, 0][order], js_chamber[band, 2][order] / 1e6, ".-", label="J_phi (toroidal)")
ax.plot(cc[band, 0][order], np.linalg.norm(js_chamber[band, :2], axis=1)[order] / 1e6 * 10, ".-",
        label="|J_poloidal| x 10")
ax.axvline(status["magnetic_axis"][0], color="k", lw=0.8, ls=":"); ax.text(status["magnetic_axis"][0], 0.1, " axis")
ax.set_xlabel("R [m]"); ax.set_ylabel("[MA/m^2]"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
ax.set_title("plasma current density on the midplane Z = Z_axis", fontsize=10); plt.show()
cs = meshes.array_values(chamber, meshes.SIZES_ARRAY, "elem")
print(f"I_p = sum J_phi a = {(js_chamber[plasma, 2] * cs[plasma]).sum() / 1e6:.3f} MA over {plasma.sum()} cells")
''')
    nb.md(r"""
### Pressure and toroidal-field function

Comparing the profile with $J_\varphi = R\,p' + FF'/(\mu_0 R)$ gives the two flux functions of the equilibrium:

$$
p'(\psi) = \frac{\lambda\beta}{R_0}\,g(\bar\psi), \qquad FF'(\psi) = \mu_0\lambda(1-\beta)R_0\,g(\bar\psi), \qquad
g = (1-\bar\psi^{\alpha_m})^{\alpha_n}.
$$

Integrating from the boundary ($p = 0$, $F = F_0$ = `AF0`), with $d\psi = (\psi_b - \psi_{axis})\,d\bar\psi$:

$$
p(\bar\psi) = \frac{\lambda\beta}{R_0}(\psi_{axis}-\psi_b)\int_{\bar\psi}^1 g\,ds, \qquad
F^2(\bar\psi) = F_0^2 + 2\mu_0\lambda(1-\beta)R_0(\psi_{axis}-\psi_b)\int_{\bar\psi}^1 g\,ds .
$$

The algorithm uses the same $F$ for the poloidal currents. From $p$ we get the stored thermal energy and the
poloidal beta $\beta_p = 2\mu_0\langle p\rangle / \bar B_p^2$, with $\bar B_p = \mu_0 I_p / L$ and $L$ the length
of the separatrix.
""")
    nb.code('''
p = algorithm.input_parameters
lam, beta, r0, am, an, f0 = (status["plasma_current_scaling_lambda"], p["beta_p"], p["R0"], p["alpha_mp"],
                             p["alpha_np"], p["AF0"])
mu0 = 4e-7 * np.pi
s = np.linspace(0.0, 1.0, 2001)
g = np.maximum(1.0 - s**am, 0.0) ** an
tail = np.concatenate(([0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(s))))
tail = tail[-1] - tail                                                       # int_s^1 g
pressure = lam * beta / r0 * (psi_axis - psi_b) * tail
F = np.sqrt(f0**2 + 2 * mu0 * lam * (1 - beta) * r0 * (psi_axis - psi_b) * tail)

fig, axes = plt.subplots(1, 2, figsize=(10, 3.4))
axes[0].plot(s, pressure / 1e6); axes[0].set_ylabel("p [MPa]")
axes[1].plot(s, F / r0, label="F / R0 = B_phi at R0"); axes[1].axhline(f0 / r0, color="k", ls=":", lw=0.8)
axes[1].set_ylabel("[T]"); axes[1].legend(fontsize=8)
for ax in axes:
    ax.set_xlabel("normalized flux psi_bar (0 axis, 1 boundary)"); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# global quantities over the plasma cells (volume element 2 pi R dA)
psi_bar = np.clip((meshes.array_values(chamber, "PSIe", "elem") - psi_axis) / (psi_b - psi_axis), 0, 1)
dv = 2 * np.pi * cc[:, 0] * cs
volume = dv[plasma].sum()
p_mean = (np.interp(psi_bar, s, pressure) * dv)[plasma].sum() / volume
scratch = plt.figure()                                  # the separatrix as a polyline, from a contour
contour = scratch.gca().tricontour(chamber_tri, history[-1]["psi"], levels=[psi_b + 1e-4 * (psi_axis - psi_b)])
length = max(np.linalg.norm(np.diff(seg, axis=0), axis=1).sum() for seg in contour.allsegs[0]); plt.close(scratch)
bp = mu0 * p["plasma_current"] / length
print(f"plasma volume {volume:.0f} m^3, cross-section {cs[plasma].sum():.1f} m^2, separatrix length {length:.2f} m")
print(f"p(axis) = {pressure[0] / 1e6:.2f} MPa, <p> = {p_mean / 1e3:.0f} kPa, W_th = 3/2 <p> V = {1.5 * p_mean * volume / 1e6:.0f} MJ")
print(f"mean poloidal field {bp:.3f} T -> beta_p = 2 mu0 <p> / Bp^2 = {2 * mu0 * p_mean / bp**2:.3f} (profile beta = {beta})")
print(f"F on axis {F[0]:.3f} T m vs F0 = {f0:.3f} T m: the plasma is {'para' if F[0] > f0 else 'dia'}magnetic")
''')
    nb.md("""
The numbers are those of the ITER reference plasma: about 830 m$^3$, a central pressure around 1 MPa, some 400 MJ
of thermal energy. Note that the profile parameter `beta_p` is *not* exactly the poloidal beta of the result:
it fixes the share of the $R/R_0$ term in the current density, while $\\beta_p$ also depends on the shape of the
plasma and on the profile exponents. To match a given $\\beta_p$ one adjusts the parameter (an outer iteration).

### The current density of the output

The solver writes the total current density of the cells in `Jt`. In a stranded conductor it is the imposed `Js`,
evaluated at every output from the current loads: after the Picard iterations `Jt` of the plasma cells holds the
converged plasma current density, like `Js`:
""")
    nb.code('''
jt = meshes.array_values(mesh, "Jt", "elem")[:, 2]
on = meshes.array_values(mesh, "PLASMA", "elem") == 1
js_solver = meshes.array_values(mesh, "Js", "elem")[:, 2]
print(f"over the PLASMA cells of the solver mesh: sum Js a = {(js_solver * area)[on].sum() / 1e6:.2f} MA, "
      f"sum Jt a = {(jt * area)[on].sum() / 1e6:.2f} MA")
''')
    nb.md("""
## A 3D view

The plasma and the coils revolved by 270 degrees about the vertical axis: the cut faces show $J_\\varphi$ in the
plasma and the flux surfaces.
""")
    nb.code('''
plasma_cells = pv.wrap(chamber).threshold(value=(1, 1), scalars="PLASMA").extract_surface(algorithm=None)
plasma_cells.cell_data["J_phi [MA/m^2]"] = plasma_cells.cell_data["Js"][:, 2] / 1e6
outline = plasma_cells.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                             non_manifold_edges=False)
shell = outline.extrude_rotate(resolution=90, angle=270, rotation_axis=(0, 1, 0), capping=False)
pl = pv.Plotter(window_size=(900, 700))
pl.add_mesh(shell, color="orchid", opacity=0.35, smooth_shading=True)
for cut in (plasma_cells, plasma_cells.rotate_y(270, point=(0, 0, 0))):
    pl.add_mesh(cut, scalars="J_phi [MA/m^2]", cmap="inferno",
                scalar_bar_args={"vertical": True, "position_x": 0.85, "position_y": 0.3})
    pl.add_mesh(cut.cell_data_to_point_data().contour(10, scalars="PSIe"), color="white", line_width=1)
for name in components[:12]:
    coil = grid.threshold(value=(1, 1), scalars=name).extract_surface(algorithm=None)
    pl.add_mesh(coil.extrude_rotate(resolution=90, angle=270, rotation_axis=(0, 1, 0), capping=True),
                color="#c87533", smooth_shading=True)
pl.camera_position = [(22, 14, 30), (0, 0, 0), (0, 1, 0)]
pl.show()
''')
    nb.md("""
## Save the result

`write_output()` writes the solver mesh with all the arrays (A, PSI, Js of the plasma and the coils, PLASMA, B...)
as `output/ITER_2D_fine2_out_<n>.vtu`, where `n` numbers the outputs of the session (the latest output plus one: the
Picard iterations and the two response solves do not count, only the outputs written). This file is what the evolution reads as its initial state (next tutorial).
""")
    nb.code('''
written = solver.write_output()
print(Path(written).name, "in", Path(written).parent.name)
with quietly():
    solver.close()
''')
    return nb


# ======================================================================================================================
def iter_evolution():
    nb = nbbuild.Notebook()
    nb.md(r"""
# ITER plasma evolution: a vertical displacement event (2D axisymmetric)

This notebook follows `examples/ITER_2D_evol_fine2`: a 15 MA equilibrium, computed on the refined `fine2` mesh as in the
previous tutorial, evolves for 0.63 s
while its poloidal beta drops from 0.68 to 0.23, at nearly constant current. Losing pressure, the plasma loses the
equilibrium imposed by the (fixed) coil currents and drifts vertically: a **vertical displacement event** (VDE).
The motion induces eddy currents in the vessel, which oppose it and set its speed; the plasma moves up, its X point
disappears, and it ends against the wall (a *limited* plasma).

**The coupled problem.** At each time step the field is a transient eddy-current solution (backward Euler):
$\nabla\times\nu\nabla\times A + \sigma\,\partial_t A = J_s$ in 2D, with $\sigma$ in the passive structures and
the plasma current density $J_s$ from the Grad-Shafranov profile, as in the equilibrium. `evolution_solve()` runs
the Picard iteration of the equilibrium **inside the time step**: each iteration recomputes the plasma current
from the last flux and re-solves *the same step* (same $t$, same previous state), until $\psi$ converges. The plasma
current and the eddy currents are thus found together, implicitly. Unlike the equilibrium there is no position
control: the plasma goes where the fields take it. The coil currents stay as they were (the model imposes their
current densities, as if ideal power supplies kept them constant).

The run takes a few minutes (the equilibrium, then 42 time steps of 15 ms). Theory: *Theory - Plasma
equilibrium and evolution*.
""")
    nb.code(SETUP.format(case="ITER_2D_evol_fine2", files='("ITER_2D_evol_fine2.vtu", "model_data.json")')
            + "\nimport json, subprocess, sys, time\nfrom pytec3d import Algorithm, Model, element_library\n" + PLOT_HELPERS)
    nb.md("""
## The initial equilibrium

The evolution starts from the converged equilibrium. We run the test macro `examples/ITER_2D_fine2/macro.py` itself, in a
sub-folder of the working folder (about 10 s), and take its output: `output/ITER_2D_fine2_out_<n>.vtu`, the most recent
of that series, as the evolution macro does. (The previous tutorial explains this computation.)
""")
    nb.code('''
EQ = WORK / "ITER_2D_fine2"
EQ.mkdir()
for name in ("ITER_2D_fine2.vtu", "model_data.json", "macro.py"):
    shutil.copy(CASE.parent / "ITER_2D_fine2" / name, EQ / name)
started = time.perf_counter()
run = subprocess.run([sys.executable, "macro.py"], cwd=EQ, capture_output=True, text=True)
assert run.returncode == 0, run.stderr[-3000:]
eq_file = max(EQ.glob("output/ITER_2D_fine2_out_*.vtu"), key=os.path.getmtime)
print(f"equilibrium computed in {time.perf_counter() - started:.0f} s:", eq_file.name)
equilibrium = Model.read_model(str(eq_file), reset_arrays=False)
''')
    nb.md("""
## The evolution model

`ITER_2D_evol_fine2.vtu` is the same mesh as the equilibrium (same cells, same component arrays), with a third TYPE: 24,
`SOLID_CONDUCTOR`, on the **passive structures**: the vacuum vessel and its supports are **solid conductors**,
massive steel (rho = 8e-7 Ohm m) in which the eddy currents are computed. In 2D a solid conductor has no electric
unknown: its current density is $-\\sigma\\,\\partial_t A$, the current of a toroidally closed shell. The coils and the
plasma region stay `STRANDED_CONDUCTOR` (imposed `Js`).
""")
    nb.code('''
solver = Solver(log_file="ITER_2D_evol_fine2.log")
algorithm = Algorithm(solver)
solver.read_model("ITER_2D_evol_fine2.vtu")
solver.model.read_json("model_data.json")
mesh = solver.model.mesh
grid = pv.wrap(mesh)
data = json.load(open("model_data.json"))
domains = dict(enumerate(element_library.domains()))
types = grid.cell_data["TYPE"]
for tid in np.unique(types):
    row, rho = data["element_types"][str(tid)], data["materials"][str(tid)]["rho"][0][0]
    print(f"TYPE {tid}: {domains[row['option1']]:26s} {np.sum(types == tid):5d} cells, rho = {rho:g} Ohm m")
structures = types == 24
area = grid.cell_data["cell_volume"]
print("passive structures = PASSIVE_STRUCTURES cells:", np.array_equal(structures, grid.cell_data["PASSIVE_STRUCTURES"] == 1))
''')
    nb.md(r"""
How fast can the vessel react? A rough time constant of a thin toroidal shell is $\tau \approx \mu_0\sigma\,d\,b$
(thickness $d$, minor radius $b$). It sets the time step: the vertical motion of an unstable plasma grows on the
time scale of the vessel currents, here a few hundred milliseconds, so steps of 15 ms resolve it, with backward
Euler (theta = 1) for stability. (On this refined mesh 30 ms steps are too coarse: around t = 0.33 s, just before
the plasma touches the wall, the Picard iteration of a 30 ms step stalls; 15 ms steps converge in at most 28
iterations.)
""")
    nb.code('''
d, b = 0.06, 2.5              # ITER vessel: about 60 mm of steel per shell, 2-2.5 m from the plasma
print(f"tau ~ mu0 sigma d b = {4e-7 * np.pi / 8e-7 * d * b * 1e3:.0f} ms")
''')
    nb.md("""
## Initial state

1. `algorithm.plasma_region_mesh = equilibrium.mesh`: the algorithm works on a copy of the PLASMA_REGION cells of
   the equilibrium, with its flux and plasma region, which is where it looks for the plasma at the first iteration.
2. $A_z = 0$ on the external boundary, as in the equilibrium.
3. The current densities of the equilibrium (coils and plasma) become the `Js` of the evolution model: the `Js` of
   the equilibrium output (its `Jt` holds the same values on the stranded conductors).
""")
    nb.code('''
algorithm.plasma_region_mesh = equilibrium.mesh
with quietly():
    solver.load_model_into_solver()
selection.select_all(mesh)
selection.nsel_ext(mesh)
selection.set_value_on_selection(mesh, "_A", 0.0, 2)

meshes.array_values(mesh, "Js", "elem")[:] = meshes.array_values(equilibrium.mesh, "Js", "elem")
meshes.get_array(mesh, "Js", "elem").Modified()
js = meshes.array_values(mesh, "Js", "elem")[:, 2]
print(f"plasma current {(js * area)[grid.cell_data['PLASMA_REGION'] == 1].sum() / 1e6:.2f} MA, "
      f"coils {(js * area)[grid.cell_data['COILS'] == 1].sum() / 1e6:.2f} MA-turns")

algorithm.input_parameters.update(R0=6.2, alpha_mp=1.0001, alpha_np=0.83197, AF0=5.3 * 6.2,
                                  tolerance_value=1e-4, max_picard_iterations=70)
''')
    nb.md("""
## The scenario

`examples/files` holds the waveforms of $I_p(t)$ and $\\beta_p(t)$, sampled every 0.2 ms up to 0.636 s. The test
macro uses every fifth sample (1 ms steps, 637 steps); here we take one every 75 (15 ms), 42 steps. At $t = 0$
the plasma loses part of its energy ($\\beta_p$ 0.68 -> 0.58 within 0.4 ms, $I_p$ -110 kA), then $\\beta_p$ keeps
decaying while $I_p$ recovers slowly to 15.2 MA. The last samples are the **thermal quench** at 0.6356 s
($\\beta_p \\to 0$, $I_p$ starting to fall; `plasma_current.txt` continues with the current quench beyond the last
time): we stop at 0.63 s, just before it.
""")
    nb.code('''
files = CASE.parent / "files"
t_all = np.loadtxt(files / "time_steps.txt")
ip_all = np.loadtxt(files / "plasma_current.txt")[:len(t_all)]
bp_all = np.loadtxt(files / "beta_p.txt")
pick = np.concatenate(([0], np.arange(75, len(t_all), 75)))
times, ip, bp = t_all[pick], ip_all[pick], bp_all[pick]
fig, ax = plt.subplots(figsize=(7, 3.2))
ax.plot(t_all, ip_all / 1e6, color="C0", lw=0.8); ax.plot(times, ip / 1e6, "o", color="C0", ms=4)
ax.set_xlabel("t [s]"); ax.set_ylabel("I_p [MA]", color="C0")
ax2 = ax.twinx(); ax2.plot(t_all, bp_all, color="C3", lw=0.8); ax2.plot(times, bp, "s", color="C3", ms=4)
ax2.set_ylabel("beta_p", color="C3"); ax.grid(alpha=0.3); ax.set_title("waveforms and the time steps used", fontsize=10)
ax.annotate("thermal quench", (0.636, 14.5), xytext=(0.43, 14.6), arrowprops={"arrowstyle": "->"}, fontsize=8)
plt.show()
print(len(times), "steps, dt =", np.unique(np.round(np.diff(times), 6)), "s")
''')
    nb.md("""
## The time loop

* Step 0 is **static**, at a tiny time: it re-solves the equilibrium field with the copied currents (no eddy
  currents) and gives the initial state of the transient. The evolution solve of the first call only initializes
  the flux, with one solve.
* Then every step is transient, theta = 1, with the waveform values of $I_p$ and $\\beta_p$ of its end time. The
  matrix is factorized once at the first transient step and reused (constant dt).

After every step we record the status (axis, X point, configuration, Picard iterations), the current induced in
the structures ($\\sum J_t\\,a$ over TYPE 24: their `Jt` is all eddy current, they have no `Js`), the current
centroid of the plasma, and the flux of the plasma region for the pictures.
""")
    nb.code('''
solver.time_step["theta"] = 1.0
chamber = algorithm.plasma_region_mesh
c_area = meshes.array_values(chamber, meshes.SIZES_ARRAY, "elem")
c_centroid = meshes.array_values(chamber, meshes.CENTROIDS_ARRAY, "elem")
records, started = [], time.perf_counter()
for k, t in enumerate(times):
    solver.time_step["transient_effects"] = k > 0
    solver.time_step["time"] = t if k > 0 else 1e-10
    algorithm.input_parameters.update(plasma_current=ip[k], beta_p=bp[k])
    with quietly():
        algorithm.evolution_solve()
    st = dict(algorithm.status)
    psi_c = meshes.array_values(chamber, "PSI", "node")
    if k == 0:      # the first call only solves: no status yet, the axis is where the flux is largest
        st.update(magnetic_axis=tuple(meshes.points(chamber)[np.argmax(psi_c)]), separatrix_configuration="initial")
    jt = meshes.array_values(mesh, "Jt", "elem")
    jp = meshes.array_values(chamber, "Js", "elem")[:, 2] * c_area
    records.append({"t": t, "iterations": st["picard_iteration"], "configuration": st["separatrix_configuration"],
                    "axis": st["magnetic_axis"][:2], "x_point": st["x_point"][:2],
                    "limiter": st["limiter_point"][:2], "psi_b": st["psi_min"],
                    "centroid": (jp @ c_centroid[:, :2]) / jp.sum() if jp.sum() else (np.nan, np.nan),
                    "vessel": (jt[structures, 2] * area[structures]).sum(),
                    "vessel_j": jt[structures, 2].copy(),
                    "psi": psi_c.copy()})
    r = records[-1]
    print(f"step {k:2d}  t = {r['t']:.2f} s  {r['iterations']:2d} Picard it.  axis ({r['axis'][0]:.2f}, {r['axis'][1]:+.2f}) m"
          f"  {r['configuration']:8s}  {st['plasma_element_count']:4d} plasma cells"
          f"  vessel {r['vessel'] / 1e6:+.2f} MA   [{time.perf_counter() - started:.0f} s]")
''')
    nb.md("""
## Results

### Trajectory and vessel current
""")
    nb.code('''
t = np.array([r["t"] for r in records])
axis = np.array([r["axis"] for r in records]); centroid = np.array([r["centroid"] for r in records])
vessel = np.array([r["vessel"] for r in records]) / 1e6
limited = np.array([r["configuration"] == "limiter" for r in records])
switch = t[np.argmax(limited)] if limited.any() else None

fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
ax = axes[0]
ax.plot(t, axis[:, 1], "o-", ms=3, label="magnetic axis Z"); ax.plot(t, centroid[:, 1], "s-", ms=3, label="current centroid Z")
ax.set_ylabel("Z [m]"); ax.legend(fontsize=8)
ax = axes[1]
ax.plot(t, axis[:, 0], "o-", ms=3, label="magnetic axis R"); ax.plot(t, centroid[:, 0], "s-", ms=3, label="current centroid R")
ax.set_ylabel("R [m]"); ax.legend(fontsize=8)
ax = axes[2]
ax.plot(t, vessel, "o-", ms=3, color="C2"); ax.set_ylabel("net vessel current [MA]")
for ax in axes:
    ax.set_xlabel("t [s]"); ax.grid(alpha=0.3)
    if switch is not None:
        ax.axvline(switch, color="k", ls=":", lw=0.8)
axes[0].set_title("vertical position", fontsize=10); axes[1].set_title("radial position", fontsize=10)
axes[2].set_title("current induced in the structures", fontsize=10)
plt.tight_layout(); plt.show()
if switch is not None:
    print(f"the plasma becomes limited at t = {switch:.2f} s (dotted line)")
''')
    nb.md("""
What happens:

* **first step**: the loss of pressure reduces the outward force on the plasma ring (the coil currents, set for
  $\\beta_p$ = 0.68, now push too much inwards): the plasma moves inwards by about 9 cm in the first 40 ms, and
  keeps drifting inwards as $\\beta_p$ decays (the magnetic axis by 0.5 m in all). The plasma current falls by about
  0.1 MA, and the vessel picks up a co-current of about 0.1 MA that keeps the flux it encloses;
* **vertical drift**: an elongated plasma in the field of fixed coils is vertically unstable. The plasma rises
  slowly at first (15 cm in 0.3 s) and then faster and faster (1.1 m in the last 0.3 s): an exponential growth whose
  rate is set by the vessel, since the plasma can only move as fast as the eddy currents that oppose the motion
  decay. Upper structures get a counter-current, lower ones a co-current (see the snapshots);
* **limited plasma**: from 0.36 s the X point is no longer the boundary: the last closed surface touches the wall
  (first the inner wall, then the top), and the plasma shrinks as it is scraped off against it.

The magnetic axis (the node of largest $\\psi$, which moves in steps of the mesh) and the current centroid move
together; the centroid is what the magnetic diagnostics of a machine measure.

### Picard iterations per step
""")
    nb.code('''
fig, ax = plt.subplots(figsize=(7, 2.8))
ax.bar(t, [r["iterations"] for r in records], width=0.02, color=np.where(limited, "C1", "C0"))
ax.set_xlabel("t [s]"); ax.set_ylabel("Picard iterations"); ax.grid(alpha=0.3, axis="y")
ax.set_title("iterations to converge each step (orange: limited)", fontsize=10); plt.show()
''')
    nb.md("""
Each Picard iteration is one solve with the factorized matrix (a forward-backward substitution), so the cost of
a step is proportional to its iterations. Steps converge slowly when the plasma moves a lot during the step or
changes configuration: smaller time steps reduce both the motion per step and the iterations.

### Snapshots

The separatrix (boundary flux surface) at a few times, over the eddy current density in the structures at the
same time (red: co-current, blue: counter-current, with the sign of the plasma current). The `PLASMA` cell array of
the chamber marks the plasma at each step: the cells inside the boundary of the current, the X-point surface of a
diverted plasma or the limiter surface of a limited one.
""")
    nb.code('''
chamber_tri = triangulation(chamber)
shown = [1, 6, 12, len(records) - 1]
fig, axes = plt.subplots(1, len(shown), figsize=(13, 5.5), sharey=True)
jmax = max(np.abs(records[k]["vessel_j"]).max() for k in shown) / 1e6
sign = np.sign(ip[0])
for ax, k in zip(axes, shown):
    r = records[k]
    pc = fill(ax, mesh, structures, np.zeros(grid.n_cells), cmap="coolwarm", edgecolor="none")
    values = np.zeros(grid.n_cells); values[structures] = sign * r["vessel_j"] / 1e6
    pc.set_array(values[structures]); pc.set_clim(-jmax / 3, jmax / 3)
    ax.tricontour(chamber_tri, r["psi"], levels=[r["psi_b"]], colors="k", linewidths=1.3, linestyles="solid")
    ax.plot(*r["axis"], "ko", ms=4)
    if r["configuration"] == "X-point":
        ax.plot(*r["x_point"], "rx", ms=8, mew=2)
    else:
        ax.plot(*r["limiter"], "gs", ms=6)
    frame(ax, r=(3.5, 9.5), z=(-5.5, 6.0)); ax.set_title(f"t = {r['t']:.2f} s, {r['configuration']}", fontsize=10)
    if ax is not axes[0]:
        ax.set_ylabel("")
fig.colorbar(pc, ax=axes, shrink=0.7, label="eddy current density x sign(Ip) [MA/m^2]")
plt.show()
''')
    nb.md("""
## A 3D view of the last step

The vessel structures revolved around the machine axis (240 degrees), coloured by the induced current density,
and the plasma (its current-carrying cells) at the end of the event, pressed against the top of the chamber.
""")
    nb.code('''
def revolved(surface, angle=240, capping=False):
    """The surface of revolution of the boundary of a 2D region (its point data carried around)."""
    outline = surface.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                            non_manifold_edges=False)
    return outline.extrude_rotate(resolution=80, angle=angle, rotation_axis=(0, 1, 0), capping=capping)

last = grid.copy()
last.cell_data["J vessel [MA/m^2]"] = sign * meshes.array_values(mesh, "Jt", "elem")[:, 2] / 1e6
wall = last.threshold(value=(24, 24), scalars="TYPE").extract_surface(algorithm=None).cell_data_to_point_data()
current = pv.wrap(chamber).copy()
current.cell_data["J_phi"] = meshes.array_values(chamber, "Js", "elem")[:, 2]
plasma_now = current.extract_cells(np.flatnonzero(current.cell_data["J_phi"] != 0)).extract_surface(algorithm=None)
pl = pv.Plotter(window_size=(900, 700))
pl.add_mesh(revolved(wall), scalars="J vessel [MA/m^2]", cmap="coolwarm", clim=[-jmax / 3, jmax / 3],
            smooth_shading=True, scalar_bar_args={"vertical": True, "position_x": 0.85, "position_y": 0.3})
pl.add_mesh(revolved(plasma_now, capping=True), color="orchid", smooth_shading=True)
pl.camera_position = [(22, 12, 30), (0, 0.5, 0), (0, 1, 0)]
pl.show()
with quietly():
    solver.close()
''')
    nb.md("""
## Going further

* **The full test**: `examples/ITER_2D_evol_fine2/macro.py` runs the same scenario with 1 ms steps and writes every step to
  `./output` (open the series in ParaView; see *User manual - ParaView*).
* **Halo currents**: when the plasma is limited, `halo_delta1`, `halo_delta2`, `t0_halo` define a halo region outside
  the last closed surface, whose poloidal currents close through the wall (3D models, `ITER_3D_evol*`).
* **3D structures**: `ITER_3D_evol*` couple the same 2D plasma with a 3D mesh of the vessel; the algorithm projects
  the plasma current around the torus and averages the 3D flux back to the 2D chamber.
* **Coils with power supplies**: modelling the coils as `STRANDED_COIL` with their circuit equation (terminal
  voltages) lets their currents react to the plasma motion as well (`ITER_2D_evol_coils`).
""")
    return nb


# ======================================================================================================================
def iter_3d_halo():
    nb = nbbuild.Notebook()
    nb.md(r"""
# ITER 3D: vertical displacement with halo currents and voltage-driven coils

This notebook follows `examples/ITER_3D_evol_coils_edge`: the vertical displacement event of the 2D tutorial, now with the
field solved on a **3D** model of the machine and three new physical effects:

* **3D solid conductors**: the vessel and the passive structures are `SOLID_CONDUCTOR` cells of the 3D nodal A-V
  formulation (`AV3D_edge`, edge elements for $A$, nodal $V$), with the electric potential $V$ (charge conservation). Besides the toroidal eddy
  currents they can carry **poloidal** currents, which the 2D model cannot represent. Here the vector potential uses
  **edge** (Nedelec) elements (`AV3D_edge`): its DOFs are circulations along the mesh edges, and a tree-cotree gauge
  fixes it;
* **halo currents**: when the plasma is limited, a thin layer of open field lines outside the last closed flux
  surface (the halo) carries part of the plasma current; its poloidal current flows *into the wall* and back,
  closing through the vessel. The algorithm injects it as nodal currents (`AMPS`) on the wall nodes;
* **voltage-driven coils**: the central solenoid modules and the PF coils are `STRANDED_COIL` cells with their circuit
  equation $V = R\,i + d\Phi/dt$, fed by constant voltages: their currents react to the plasma motion, instead of
  being frozen as in the 2D tutorial.

**How the 3D model is coupled with the plasma.** The plasma stays axisymmetric and lives on the 2D plasma-region
mesh of the equilibrium. The field is solved on a thin toroidal **sector** of 2 degrees (four cells along the torus,
cyclic symmetry between its faces) obtained by revolving the 2D mesh. At every Picard iteration the 2D plasma current
density (toroidal and poloidal) and the halo currents are *projected* onto every toroidal station of the sector; after
the solve, the toroidal vector potential is *averaged* over the stations to give the flux $\psi$ seen by the plasma.

The run takes about 40 minutes (the equilibrium, then 22 steps of 30 ms, each a Picard iteration of 3D solves on a
quadratic mesh of 464 000 nodes).
Read the 2D tutorials first; theory: *Theory - Plasma equilibrium and evolution* (3D coupling, halo currents),
*Circuits* (stranded coils).
""")
    nb.code(SETUP.format(case="ITER_3D_evol_coils_edge", files='("ITER_3D_evol_coils_edge.vtu", "model_data.json")')
            + "\nimport json, subprocess, sys, time, vtk\nfrom pytec3d import Algorithm, element_library\n"
            + "from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem\n"
            + PLOT_HELPERS)
    nb.md("""
## The initial equilibrium

As in the 2D evolution, we run the equilibrium macro `examples/ITER_2D_carma/macro.py` (the quadratic mesh) in a sub-folder and take
its output. Its mesh becomes the plasma-region mesh of the algorithm: the 2D mesh on which the plasma is computed.
""")
    nb.code('''
EQ = WORK / "ITER_2D_carma"
EQ.mkdir()
for name in ("ITER_2D_carma.vtu", "model_data.json", "macro.py"):
    shutil.copy(CASE.parent / "ITER_2D_carma" / name, EQ / name)
started = time.perf_counter()
run = subprocess.run([sys.executable, "macro.py"], cwd=EQ, capture_output=True, text=True)
assert run.returncode == 0, run.stderr[-3000:]
eq_file = max(EQ.glob("output/ITER_2D_carma_out_*.vtu"), key=os.path.getmtime)
print(f"equilibrium computed in {time.perf_counter() - started:.0f} s:", eq_file.name)
equilibrium = meshes.read_mesh(str(eq_file), reset_arrays=False)
''')
    nb.md("""
## The 3D sector

`ITER_3D_evol_coils_edge.vtu` is the quadratic ITER cross-section of `ITER_2D_carma` revolved by 2 degrees about the
vertical (y) axis, with four cells along the torus (five node stations, every 0.5 degree, and the midside nodes in
between), inside a 100 m sphere of air. The angle is measured
from the x axis towards -z. Cells: HEX20, with WEDGE15 on the axis: a quadratic mesh, because the edge DOFs of
`AV3D_edge` live on its midside nodes.
""")
    nb.code('''
solver = Solver(log_file="ITER_3D_evol_coils_edge.log")
algorithm = Algorithm(solver)
solver.read_model("ITER_3D_evol_coils_edge.vtu")
algorithm.plasma_region_mesh = equilibrium
solver.model.read_json("model_data.json")
mesh = solver.model.mesh
grid = pv.wrap(mesh)
centroid = grid.cell_data["cell_centroids"]
theta = np.arctan2(-centroid[:, 2], centroid[:, 0])          # toroidal angle of the cell centres [rad]
radius = np.hypot(centroid[:, 0], centroid[:, 2])
angles = np.degrees(np.arctan2(-grid.points[:, 2], grid.points[:, 0]))
print(f"{grid.n_points} points, {grid.n_cells} cells:", dict(zip(*np.unique(grid.celltypes, return_counts=True))),
      "(25 = HEX20, 26 = WEDGE15)")
print("node stations [deg]:", np.unique(np.round(angles[np.hypot(grid.points[:, 0], grid.points[:, 2]) > 0.1], 3)))

data = json.load(open("model_data.json"))
domains = dict(enumerate(element_library.domains()))
for tid, row in data["element_types"].items():
    print(f"TYPE {tid}: AV3D_edge {domains[row['option1']]:26s} {np.sum(grid.cell_data['TYPE'] == int(tid)):6d} cells")
''')
    nb.code('''
machine = grid.threshold(value=(32.5, 35.5), scalars="TYPE")         # all but the air
pl = pv.Plotter(window_size=(900, 600))
colors = {33: "#dbe9f6", 34: "dimgray", 35: "#c87533"}
for tid, color in colors.items():
    part = machine.threshold(value=(tid, tid), scalars="TYPE")
    if part.n_cells:
        pl.add_mesh(part, color=color, show_edges=True, edge_color="gray", line_width=0.3,
                    label={33: "plasma region", 34: "vessel and structures",
                           35: "coils (stranded coils)"}[tid])
pl.add_legend(bcolor="white", size=(0.32, 0.15))
pl.camera_position = [(16, 4, -14), (6, 0, -0.1), (0, 1, 0)]
pl.show()
''')
    nb.md("""
The picture shows the 2 degree slice of the machine, seen from outside: thin, but a true 3D mesh. The coil cells
already have TYPE 35 in the file; the macro sets it again on the `COILS` cells, with their element coordinate
system, in the next section.

## Voltage-driven coils

Each CS module and PF coil becomes a `STRANDED_COIL` (TYPE 35): a winding of N turns with a current $i$, uniform
over its cross-section, along the direction `T` of the real constants. Its circuit equation links the terminal
voltage `VOLT`, the resistive drop and the induced voltage `EMF` (the time derivative of the flux linked with the
turns):

$$\\text{VOLT} = R\\,i + \\text{EMF}.$$

Points to note for a 3D sector (real constants of each coil):

* `T = (0, 1, 0)` is the **toroidal** direction in the cylindrical element coordinate system 1 (y axis), which the
  coil cells must use (`COORDINATE_SYSTEM = 1`); the plasma region keeps the global system (the algorithm writes
  global components in its `Js`);
* `volume` is the volume of the **whole** coil and `symmetry = 180` (360 / 2 degrees): the model holds 1/180 of each
  coil, and the circuit equation is that of the whole coil (R and N of the whole coil, the EMF of the sector times
  180);
* all the nodes of a coil share one `VOLT` and one `EMF`: a CP id per coil in `CP_0` (1-6, 11-16), and the couplings
  set with `add_dof_couplings` (the model data has none).

The source voltage of each coil is the one that keeps its equilibrium current in steady state, `VOLT = R i`, with
`i = S Js / N` from the current density of the coil in the equilibrium output (its `Jt` holds the same values: the
table shows both). For the four control coils this includes the correction currents of the position control.
""")
    nb.code('''
coil_names = ["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L"] + [f"PFC_{n}" for n in range(1, 7)]
coil_ids = dict(zip(coil_names, list(range(11, 17)) + list(range(1, 7))))   # CP id = real constants row

selection.select(mesh, "COILS", 1, location="elem")
selection.set_value_on_selection(mesh, "TYPE", 35, location="elem")
selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")
solver.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
selection.select(mesh, "PLASMA_REGION", 1, location="elem")
selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 0, location="elem")

meshes.remove_array(mesh, "CP_0")
meshes.add_array(mesh, "CP_0", dtype=vtk.VTK_INT)
solver.model.add_dof_couplings({cid: ["VOLT", "EMF"] for cid in coil_ids.values()})
eq_js, eq_jt = meshes.array_values(equilibrium, "Js", "elem"), meshes.array_values(equilibrium, "Jt", "elem")
coils = {}
print(f"{'coil':>6} {'N':>5} {'R [Ohm]':>8} {'i [kA]':>8} {'i from Jt':>9} {'VOLT [V]':>9}")
for name, cid in coil_ids.items():
    selection.select(mesh, name, 1, location="elem")
    selection.nsle(mesh)
    selection.set_value_on_selection(mesh, "CP_0", cid)
    master = int(selection.selected_ids(mesh, "node")[0])          # the lowest node id: the master of the coupling
    turns, resistance, section_area = (solver.model.real_constant(cid, key)
                                       for key in ("turns", "resistance", "cross_section"))
    on = meshes.array_values(equilibrium, name, "elem") == 1
    current = section_area * eq_js[on, 2].mean() / turns            # A per turn
    selection.select(mesh, "ID", master)
    selection.set_value_on_selection(mesh, "_VOLT", resistance * current)
    coils[name] = {"master": master, "turns": turns, "resistance": resistance, "i0": current}
    print(f"{name:>6} {turns:5.0f} {resistance:8.4f} {current / 1e3:8.3f} "
          f"{section_area * eq_jt[on, 2].mean() / turns / 1e3:9.3f} {resistance * current:9.2f}")
''')
    nb.md("""
## Cyclic symmetry, boundary conditions and gauge

* **Cyclic symmetry**: every node of the face at 0 degrees is coupled with its image at 2 degrees, for `A` and `V`
  (`periodic_boundaries`, in a second CP array, `CP_1`, since a node of a coil face already has a `CP_0` id). With
  edge elements the coupled values are scalars (circulations along the edges), so no nodal frame rotation is needed.
* **Far boundary**: the tangential $A$ vanishes on the sphere of radius 100 m: the edge DOFs, stored as component 2 of
  `_A` on the midside nodes, are fixed to zero there.
* **Gauge**: the edge discretisation of $A$ has no Coulomb gauge term; a tree-cotree decomposition of the mesh edges
  fixes $A = 0$ on the edges of a spanning tree (`tree_gauge`), built consistently with the periodic pairs.
* **Potential of the structures**: $V$ is defined up to a constant in each solid conductor not connected to a
  source: it is fixed to zero on one node of the vessel (node 16170).
""")
    nb.code('''
from pytec3d import tree_gauge
eps, sector, divisions = 0.1, 2.0, 4
cylindrical = YAxisCylindricalCoordinateSystem()
meshes.remove_array(mesh, "CP_1")
meshes.add_array(mesh, "CP_1", dtype=vtk.VTK_INT)
selection.select_pos(mesh, 1, (-eps, eps), coordinate_system=cylindrical)
pairs = solver.model.periodic_boundaries(cylindrical, dy=sector, array_name="CP_1", variables=("A", "V"))
print(len(pairs), "periodic node pairs")

selection.select_all(mesh)
selection.select_pos(mesh, 0, (90.0, 110.0), coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(mesh, "_A", 0.0, 2)                 # tangential A = 0 (edge DOFs, midside nodes)
tree_gauge.tree_gauge(mesh, periodic_point_coupling="CP_1")
selection.select(mesh, "ID", 16170)
selection.set_value_on_selection(mesh, "_V", 0.0)

with quietly():
    solver.load_model_into_solver()
''')
    nb.md(r"""
## Algorithm parameters: 3D coupling and halo

Besides the profile parameters of the 2D evolution:

* `extrusion_angle = 2`, `n_toroidal_divisions = 4`: the sector, to build the toroidal maps between the 2D chamber
  and the 3D stations;
* the **halo**: from `t0_halo` = 0.41 s on, when the plasma is limited, the current-carrying region extends beyond the
  last closed surface $\psi_b$ to

  $$\psi_2 = \psi_b - \Delta_1\,(1 + (t - t_0)\,\Delta_2)\,(\psi_{axis} - \psi_b),$$

  with `halo_delta1` $=\Delta_1 = 0.04$ (a halo 4 % of the core flux wide) and `halo_delta2` $= \Delta_2 = 0$ (constant
  width). The same profile applies over core and halo, with $\bar\psi$ normalized on $\psi_2$; along the wall wetted
  by the halo, the poloidal current that crosses the wall between two nodes follows from the change of $F(\psi)$
  between them, and is injected as `AMPS` (nodal current of the potential $V$ of the vessel).
""")
    nb.code('''
algorithm.input_parameters.update(R0=6.2, alpha_mp=1.0001, alpha_np=0.83197, AF0=5.3 * 6.2,
                                  tolerance_value=1e-4, max_picard_iterations=70,
                                  extrusion_angle=sector, n_toroidal_divisions=divisions,
                                  halo_delta1=0.04, halo_delta2=0.0, t0_halo=0.41)
files = CASE.parent / "files"
t_all = np.loadtxt(files / "time_steps.txt")
ip_all = np.loadtxt(files / "plasma_current.txt")[:len(t_all)]
bp_all = np.loadtxt(files / "beta_p.txt")
pick = np.concatenate(([0], np.arange(150, len(t_all), 150)))       # 30 ms steps
times, ip, bp = t_all[pick], ip_all[pick], bp_all[pick]
print(len(times), "steps up to", times[-1], "s")
''')
    nb.md("""
## Post-processing tools

The results are 3D vectors in global components. For pictures and integrals we use:

* `section(grid)`: the cut of the 3D grid by the meridian plane at 0.25 degrees (through the first layer of cells),
  flattened to (R, Z), with the vector cell arrays in **cylindrical** components (R, Z, toroidal) - the toroidal
  sign as in the 2D models, so the plasma current is positive;
* integrals over the first layer of cells (0 - 0.5 degree), times 720 for the whole torus: the toroidal current of a
  region is $\\sum J_\\varphi\\,a$ with the cell area $a = V_{cell} / (R\\,\\Delta\\theta)$;
* the **vertical force** on the structures, $F_Z = \\int (\\mathbf J\\times\\mathbf B)_Z\\,dV$. The model has no
  toroidal-field coils: its $\\mathbf B$ has only the poloidal field (and the small toroidal field of the plasma
  poloidal currents). The force of the poloidal (halo) currents in the vessel comes from the vacuum toroidal field
  $B_\\varphi = F_0/R$ (`AF0` = 5.3 T x 6.2 m), which we add here: $F_Z^{halo} = -\\int J_R\\,B_\\varphi\\,dV$.
""")
    nb.code('''
def cylindrical_components(vectors, angle):
    """(R, Z, toroidal) components of global vectors at toroidal angles `angle` (toroidal axis = local z at 0 deg)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.c_[vectors[:, 0] * c - vectors[:, 2] * s, vectors[:, 1], vectors[:, 0] * s + vectors[:, 2] * c]

def section(grid, angle_deg=0.25):
    """The meridian cut at angle_deg, flattened to (R, Z), with Jt, Js, B in cylindrical components (*_cyl)."""
    a = np.radians(angle_deg)
    cut = pv.wrap(grid).slice(normal=(-np.sin(a), 0.0, -np.cos(a)), origin=(0.0, 0.0, 0.0))
    points = cut.points.copy()
    cut.points = np.c_[np.hypot(points[:, 0], points[:, 2]), points[:, 1], np.zeros(len(points))]
    for name in ("Jt", "Js", "B"):
        cut.cell_data[name + "_cyl"] = cylindrical_components(cut.cell_data[name], a)
    return cut.cast_to_unstructured_grid()

dtheta = sector / divisions                            # angle of a cell layer [deg]
layer = theta < np.radians(dtheta)                     # first layer of cells
area = grid.cell_data["cell_volume"] / (radius * np.radians(dtheta))
volume = grid.cell_data["cell_volume"]
types = grid.cell_data["TYPE"]
structure = layer & (types == 34)
b_phi_vacuum = algorithm.input_parameters["AF0"] / radius

def integrals():
    """Toroidal currents of the plasma region, structures and coils, and vertical forces on the structures (torus)."""
    jt = cylindrical_components(meshes.array_values(mesh, "Jt", "elem"), theta)
    js = cylindrical_components(meshes.array_values(mesh, "Js", "elem"), theta)
    j_global, b_global = meshes.array_values(mesh, "Jt", "elem"), meshes.array_values(mesh, "B", "elem")
    f_model = np.cross(j_global, b_global)[:, 1] * volume
    return {"I_plasma": (js[:, 2] * area)[layer & (types == 33)].sum(),
            "I_vessel": (jt[:, 2] * area)[structure].sum(),
            "Fz_model": f_model[structure].sum() * 360 / dtheta,
            "Fz_halo": (-jt[:, 0] * b_phi_vacuum * volume)[structure].sum() * 360 / dtheta,
            "jpol_max": np.hypot(jt[:, 0], jt[:, 1])[structure].max()}
''')
    nb.md("""
## The time loop

As in 2D: a static step 0 (the coils at their source voltages, `VOLT = R i`, give back the equilibrium currents),
then transient backward-Euler steps of 30 ms, with an evolution solve (Picard iterations of 3D solves) per step. After
each step we record the status, the coil currents $i = (\\text{VOLT} - \\text{EMF})/R$ at the master node of each
coil, the integrals above, the halo injection of the 2D chamber (`AMPS` of `algorithm.plasma_region_mesh`), the
flux and current of the chamber, and the meridian section.
""")
    nb.code('''
solver.time_step["theta"] = 1.0
chamber = algorithm.plasma_region_mesh
records, started = [], time.perf_counter()
for k, t in enumerate(times):
    solver.time_step["transient_effects"] = k > 0
    solver.time_step["time"] = t if k > 0 else 1e-10
    algorithm.input_parameters.update(plasma_current=ip[k], beta_p=bp[k])
    step_start = time.perf_counter()
    with quietly():
        algorithm.evolution_solve()
    st = dict(algorithm.status)
    psi_c = meshes.array_values(chamber, "PSI", "node")
    if k == 0:
        st.update(magnetic_axis=tuple(meshes.points(chamber)[np.argmax(psi_c)]), separatrix_configuration="initial")
    volt, emf = meshes.array_values(mesh, "VOLT", "node"), meshes.array_values(mesh, "EMF", "node")
    amps = meshes.array_values(chamber, "AMPS", "node")
    records.append({"t": t, "status": st, "seconds": time.perf_counter() - step_start,
                    "coil_i": {n: (volt[c["master"]] - emf[c["master"]]) / c["resistance"] for n, c in coils.items()},
                    "psi": psi_c.copy(), "psie": meshes.array_values(chamber, "PSIe", "elem").copy(),
                    "js": meshes.array_values(chamber, "Js", "elem").copy(),
                    "amps": amps.copy(), "halo": amps[amps > 0].sum() * 360 / sector,
                    "section": section(mesh), **integrals()})
    r = records[-1]
    print(f"step {k:2d} t = {t:.2f} s {st['picard_iteration']:2d} it. axis ({st['magnetic_axis'][0]:.2f}, "
          f"{st['magnetic_axis'][1]:+.2f}) {st['separatrix_configuration']:8s} vessel {r['I_vessel'] / 1e6:+.3f} MA"
          f"  halo {r['halo'] / 1e3:5.1f} kA  [{time.perf_counter() - started:4.0f} s]")
''')
    nb.md("""
## Results

### Overview

The axis trajectory, the toroidal current of the structures, the change of the coil currents (ampere-turns), the
vertical force on the structures (eddy and halo terms), the Picard iterations per step, and the halo current
entering the wall (whole torus) with the largest poloidal current density in the structures. Dotted line: first
limited step; shaded: halo active.
""")
    nb.code('''
t = np.array([r["t"] for r in records])
axis = np.array([r["status"]["magnetic_axis"][:2] for r in records])
config = [r["status"]["separatrix_configuration"] for r in records]
limited = np.array([c == "limiter" for c in config])
halo_on = np.array([np.abs(r["amps"]).sum() > 0 for r in records])
fig, axes = plt.subplots(2, 3, figsize=(15, 6.5), sharex=True)
ax = axes[0, 0]
ax.plot(t, axis[:, 1], "o-", ms=3, label="Z axis"); ax.plot(t, axis[:, 0] - axis[0, 0], "s-", ms=3, label="R axis - R(0)")
ax.set_ylabel("[m]"); ax.legend(fontsize=8); ax.set_title("magnetic axis", fontsize=10)
ax = axes[0, 1]
ax.plot(t, [r["I_vessel"] / 1e6 for r in records], "o-", ms=3, color="C2")
ax.set_ylabel("[MA]"); ax.set_title("toroidal current in the structures", fontsize=10)
ax = axes[1, 0]
for name in coil_names:
    di = np.array([r["coil_i"][name] - coils[name]["i0"] for r in records]) * coils[name]["turns"] / 1e3
    ax.plot(t, di, lw=1.2, label=name.replace("PFC_", "PF"), ls="-" if name.startswith("PF") else "--")
ax.set_ylabel("change of N i [kA-turns]"); ax.set_title("coil currents", fontsize=10); ax.legend(fontsize=7, ncol=3)
ax = axes[1, 1]
ax.plot(t, [r["Fz_model"] / 1e6 for r in records], "o-", ms=3, label="J x B (model field)")
ax.plot(t, [r["Fz_halo"] / 1e6 for r in records], "s-", ms=3, label="J_pol x B_phi (vacuum toroidal field)")
ax.set_ylabel("[MN]"); ax.set_title("vertical force on the structures", fontsize=10); ax.legend(fontsize=8)
ax = axes[0, 2]
ax.bar(t, [r["status"]["picard_iteration"] for r in records], width=0.02, color=np.where(limited, "C1", "C0"))
ax.set_ylabel("Picard iterations"); ax.set_title("iterations per step (orange: limited)", fontsize=10)
ax = axes[1, 2]
ax.plot(t, [r["halo"] / 1e3 for r in records], "o-", ms=3, color="C2"); ax.set_ylabel("halo current [kA]", color="C2")
ax2 = ax.twinx(); ax2.plot(t, [r["jpol_max"] / 1e6 for r in records], "s-", ms=3, color="C4")
ax2.set_ylabel("max |J_pol| structures [MA/m^2]", color="C4"); ax.set_title("halo and poloidal currents", fontsize=10)
for ax in axes.flat:
    ax.grid(alpha=0.3)
    if limited.any():
        ax.axvline(t[np.argmax(limited)], color="k", ls=":", lw=0.8)
    if halo_on.any():
        ax.axvspan(t[np.argmax(halo_on)], t[-1], color="orange", alpha=0.1)
for ax in axes[1]:
    ax.set_xlabel("t [s]")
plt.tight_layout(); plt.show()
''')
    nb.md("""
How to read it (and how it compares with the 2D tutorial):

* **Motion.** The same vertical drift, a little slower: the plasma becomes limited at 0.39 s instead of 0.36 s and the
  axis ends at Z = 1.6 m instead of 1.8 m. The structures have the same resistivity, but the coils now react, the
  structures can carry poloidal currents, and the plasma sees the toroidal average of the 3D flux.
* **Vessel current and forces.** The toroidal current induced in the structures grows to about 0.22 MA. Its
  interaction with the poloidal field gives a vertical force on the structures that grows with the displacement, to
  about 27 MN at the end: the order of magnitude of the design loads of the ITER vessel for VDEs.
* **Coils.** With constant source voltages and a tiny resistance, each coil keeps (almost) the flux it links: as the
  plasma moves, its currents change by hundreds of kA-turns, up to 0.6 MA-turns in CS2U (about 15 % of its current).
  The upper CS modules (CS3U, CS2U), which see the plasma come closer, lose current; the coils it moves away from
  gain current. In the 2D tutorial the coil currents were frozen: here they react, through their circuit
  equations, and their response feeds back on the plasma motion. A real machine adds the power supplies and the
  vertical-stability controller, which change the voltages.
* **Poloidal currents in the structures.** Even before the halo the structures carry poloidal currents of about
  0.1 MA/m$^2$. They are induced by the change of the *toroidal* flux of the plasma: its poloidal currents (the
  paramagnetic $F(\\psi) > F_0$ of the equilibrium) decrease with $\\beta_p$. They peak (0.2 MA/m$^2$) right after
  the fast drop of $\\beta_p$ at $t = 0$, and decay with the time constant of the vessel. A 2D model, with only
  toroidal currents, cannot have them.
* **Halo.** From 0.42 s the halo is active and its poloidal current is injected into the wall: some 13 - 28 kA for
  the whole torus, a small fraction of $I_p$. With `halo_delta1 = 0.04` the halo is thin, and the profile
  $(1 - \\bar\\psi^{\\alpha_m})^{\\alpha_n}$ normalized on $\\psi_2$ is close to zero there, so little current flows in
  it. The halo force ($J_{pol}\\times B_\\varphi$ with the 5.3 T vacuum field) is negligible here.
  Wider halos (larger `halo_delta1`, growth with `halo_delta2`) or broader edge profiles give the halo fractions of
  tens of per cent seen in VDEs, and forces comparable with the eddy ones.
* **Cost.** Each Picard iteration is a 3D solve with the factorized matrix plus the projections and averages
  between the 2D chamber and the 3D stations: about 5 - 6 s here, against 0.15 s in 2D.
""")
    nb.md("""
### Core, halo and wall currents

Snapshots of the meridian section: the plasma core (inside $\\psi_b$, coloured by $J_\\varphi$) and the halo
(between $\\psi_b$ and $\\psi_2$, green), the toroidal current density in the structures, and the halo injection
points on the wall (`AMPS`: red where the current enters the wall, blue where it leaves, marker area proportional to
the current).
""")
    nb.code('''
chamber_tri = triangulation(chamber)
chamber_points = meshes.points(chamber)[:, :2]
shown = [11, 14, len(records) - 1]
fig, axes = plt.subplots(1, len(shown), figsize=(13.5, 5.6), sharey=True)
jmax = max(np.abs(records[k]["section"].cell_data["Jt_cyl"][records[k]["section"].cell_data["TYPE"] == 34, 2]).max()
           for k in shown) / 1e6
for ax, k in zip(axes, shown):
    r, st = records[k], records[k]["status"]
    sec = r["section"]
    walls = sec.cell_data["TYPE"] == 34
    pc = fill(ax, sec, walls, sec.cell_data["Jt_cyl"][:, 2] / 1e6, cmap="coolwarm", edgecolor="none")
    pc.set_clim(-jmax / 2, jmax / 2)
    carrying = r["js"][:, 2] != 0
    core = carrying & (r["psie"] >= st["psi_min"])
    halo = carrying & ~core
    fill(ax, chamber, core, r["js"][:, 2] / 1e6, cmap="inferno", edgecolor="none")
    if halo.any():
        fill(ax, chamber, halo, facecolor="limegreen", edgecolor="none")
    ax.tricontour(chamber_tri, r["psi"], levels=[st["psi_min"]], colors="k", linewidths=1.0, linestyles="solid")
    if st["psi_2_min"] < st["psi_min"]:
        ax.tricontour(chamber_tri, r["psi"], levels=[st["psi_2_min"]], colors="green", linewidths=0.8,
                      linestyles="solid")
    hit = r["amps"] != 0
    if hit.any():
        size = 120 * np.abs(r["amps"][hit]) / np.abs(r["amps"]).max()
        ax.scatter(*chamber_points[hit].T, s=size, c=np.where(r["amps"][hit] > 0, "red", "blue"), zorder=5)
    frame(ax, r=(3.5, 9.5), z=(-5.0, 5.5))
    ax.set_title(f"t = {r['t']:.2f} s, {st['separatrix_configuration']}" + (", halo" if halo.any() else ""),
                 fontsize=10)
    if ax is not axes[0]:
        ax.set_ylabel("")
fig.colorbar(pc, ax=axes, shrink=0.7, label="toroidal current density in the structures [MA/m^2]")
plt.show()
''')
    nb.md("""
### The halo current through the wall

A zoom on the contact region at the last step, with the **poloidal** current density in the structures (the R and Z
components of `Jt`, colour and arrows). The halo current enters the wall above the contact point (red) and leaves it
below (blue), flowing along the vessel shell in between: a current path that only the 3D A-V model, with the
potential $V$, can carry. Here it is small against the poloidal currents induced by the change of toroidal flux: 13 -
28 kA spread over a wall section of about 1.5 m$^2$ (the shell thickness times $2\\pi R$) is some 0.02 MA/m$^2$,
against 0.1 MA/m$^2$ of induced current.
""")
    nb.code('''
r, st = records[-1], records[-1]["status"]
sec = r["section"]
walls = sec.cell_data["TYPE"] == 34
j_pol = sec.cell_data["Jt_cyl"][:, :2]
centres = sec.cell_centers().points[:, :2]
contact = np.array(st["limiter_point"][:2])
window = walls & (np.abs(centres[:, 0] - contact[0]) < 0.6) & (np.abs(centres[:, 1] - contact[1]) < 1.3)
window[np.flatnonzero(window)[1::2]] = False                     # every second cell: readable arrows
fig, ax = plt.subplots(figsize=(5, 7))
pc = fill(ax, sec, walls, np.linalg.norm(j_pol, axis=1) / 1e6, cmap="viridis", edgecolor="none")
fig.colorbar(pc, ax=ax, shrink=0.8, label="|J poloidal| in the structures [MA/m^2]")
carrying = r["js"][:, 2] != 0
core = carrying & (r["psie"] >= st["psi_min"])
fill(ax, chamber, core, facecolor="#f3c7e6", edgecolor="none")
fill(ax, chamber, carrying & ~core, facecolor="limegreen", edgecolor="none")
ax.tricontour(chamber_tri, r["psi"], levels=np.sort(st["psi_min"] + np.linspace(-0.15, 0.6, 12)
              * (st["psi_magnetic_axis"] - st["psi_min"])), colors="0.5", linewidths=0.5, linestyles="solid")
ax.quiver(*centres[window].T, *j_pol[window].T, color="k", angles="xy", scale_units="xy",
          scale=np.linalg.norm(j_pol[window], axis=1).max() / 0.2, width=0.004, headwidth=4)
hit = r["amps"] != 0
ax.scatter(*chamber_points[hit].T, s=40 * np.abs(r["amps"][hit]) / np.abs(r["amps"]).max() + 4,
           c=np.where(r["amps"][hit] > 0, "red", "blue"), zorder=5)
ax.plot(*contact, "s", color="orange", ms=8, label="limiter (contact) point")
frame(ax, r=(contact[0] - 0.45, contact[0] + 0.85), z=(contact[1] - 1.2, contact[1] + 1.2))
ax.legend(fontsize=8, loc="lower left"); ax.set_title(f"t = {r['t']:.2f} s: halo (green) and wall currents", fontsize=10)
plt.show()
''')
    nb.md("""
## A 3D view

The last step: the vessel revolved by 240 degrees, coloured by the poloidal current density; the last closed flux
surface (violet) and the halo boundary $\\psi_2$ (green, translucent) as surfaces of revolution; the arrows on the 3D sector itself (at 0 - 2
degrees, the cut face on the right) show the direction of the poloidal current density in the structure cells
near the contact point.
""")
    nb.code('''
def revolved(surface, angle=240, capping=False):
    outline = surface.extract_feature_edges(boundary_edges=True, feature_edges=False, manifold_edges=False,
                                            non_manifold_edges=False)
    return outline.extrude_rotate(resolution=80, angle=angle, rotation_axis=(0, 1, 0), capping=capping)

sec = records[-1]["section"].copy()
sec.cell_data["|J pol| [MA/m^2]"] = np.linalg.norm(sec.cell_data["Jt_cyl"][:, :2], axis=1) / 1e6
wall = sec.extract_cells(np.flatnonzero(sec.cell_data["TYPE"] == 34)).extract_surface(algorithm=None)
wall = wall.cell_data_to_point_data()
plasma_2d = pv.wrap(chamber).copy()
core_cells = np.flatnonzero(core)
halo_cells = np.flatnonzero(carrying & ~core)
pl = pv.Plotter(window_size=(950, 720))
pl.add_mesh(revolved(wall), scalars="|J pol| [MA/m^2]", cmap="viridis", smooth_shading=True,
            scalar_bar_args={"vertical": True, "position_x": 0.86, "position_y": 0.3})
plasma_2d.point_data["PSI"] = records[-1]["psi"]
def flux_surface(level):
    """The flux surface psi = level revolved by 240 degrees (from the contour line of the 2D chamber)."""
    line = plasma_2d.contour([level], scalars="PSI")
    return line.extrude_rotate(resolution=80, angle=240, rotation_axis=(0, 1, 0), capping=False)
pl.add_mesh(flux_surface(st["psi_min"]), color="orchid", smooth_shading=True)
if halo_cells.size:
    pl.add_mesh(flux_surface(st["psi_2_min"]), color="limegreen", opacity=0.35, smooth_shading=True)
structure_cells = grid.extract_cells(np.flatnonzero((types == 34)
                                     & (np.hypot(radius - contact[0], centroid[:, 1] - contact[1]) < 2.0)))
centres_3d = structure_cells.cell_centers()
near = np.flatnonzero((types == 34) & (np.hypot(radius - contact[0], centroid[:, 1] - contact[1]) < 2.0))
j3 = structure_cells.cell_data["Jt"]
phi_hat = np.c_[np.sin(theta[near]), np.zeros(near.size), np.cos(theta[near])]
centres_3d["J_pol"] = j3 - np.sum(j3 * phi_hat, axis=1)[:, None] * phi_hat        # the poloidal part of Jt
pl.add_mesh(centres_3d.glyph(orient="J_pol", scale=False, factor=0.25), color="white")
pl.camera_position = [(19, 13, 24), (2, 1.5, 0), (0, 1, 0)]
pl.show()
with quietly():
    solver.close()
''')
    nb.md("""
## Notes and going further

* **The full test** (`examples/ITER_3D_evol_coils_edge/macro.py`) runs the same model with 1 ms steps (637 steps), a long
  run.
* **Halo injection on the stations.** The algorithm computes the halo current of the sector (scaled by
  `extrusion_angle / 360`) and shares it among the five node stations of the wall (every 0.5 degree) with the
  trapezoid weights 1/4, 1/2, 1/4: stations 0 and 2 are one periodic pair, whose loads add up, so the sector receives
  its share once (*Theory - Halo currents*).
* **Other 3D cases**: `ITER_3D_evol` (coils with fixed currents), `ITER_3D_evol_edge` and `ITER_3D_evol_coils_edge`
  (edge elements and the tree-cotree gauge), `ITER_3D_evol_fine*` (finer meshes).
* **Non-axisymmetric structures**: this sector is axisymmetric. A sector with ports or blanket modules meshed in 3D
  carries fully 3D eddy and halo current paths, coupled with the same 2D plasma (*Theory - 3D coupling*).
""")
    return nb


if __name__ == "__main__":
    target = os.path.abspath(sys.argv[1])     # the notebooks change the working directory
    which = sys.argv[2:] or ["iter_equilibrium", "iter_evolution", "iter_3d_halo"]
    for name in which:
        nbbuild.run(globals()[name](), f"{target}/tutorial_{name}.ipynb")

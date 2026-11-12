# =====================================================================================================================
# wire_2D_tphi - 2D H-formulation (TPHI2D_edge): AC transport loss of a superconducting wire, power law E-J
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The cross-section of a long round superconducting wire carrying the AC transport current I(t) = F Ic sin(2 pi t / T)
# (F = 0.5, Ic = Jc pi R^2), solved with TPHI2D_edge (planar): T in the plane (edge, NEDELEC_ONE on the QUAD8 cells),
# PHI nodal, J = curl T and E out of the plane (component 2 of Jt and E). The loss per cycle is compared with Norris's
# critical-state formula for a round wire.
#   * the power law E = Ec (J / Jc)^n + rho0 J as a TABLE (set_power_law: law power_law, rho(J) = E(J) / J), ending
#     with a normal-state segment above E = 1e2 Ec (current sharing): the curve goes on along its last segment, and a
#     pure power law extrapolated by a Newton iterate reaches tangents of 1e30 where the factorization fails;
#   * the nonlinear settings for high n: an ABSOLUTE residual tolerance scaled to the problem (a T row is about
#     Ec curl(w) A, here 1e-9) with the relative one off (a relative tolerance is relative to the first residual of the
#     step, huge when the starting guess has J > Jc), and require_residual_reduction (a transport current changing
#     at the cut at every step);
#   * the CUT of a 2D conductor with a net current: a one-cell strip of CUT_LAYER cells from the wire to the box along
#     the ray y = 0+, x > R. _T = 0 on its edges along the ray and on its other side, _T = +-I on its crossing edges
#     (the sign from the solver's orientation of each edge: tree_gauge.edges), the crossing edge on the wire surface
#     included (it belongs to the cut, not to n x T = 0);
#   * the box carries the exact potential outside the wire, PHI = -I theta / (2 pi); PHI jumps by I across the cut,
#     but it is one nodal field: the jump is across the strip (T carries it there), so the nodes on the ray take
#     theta = 2 pi and those just above it theta about 0;
#   * the loss from the E and Jt outputs: Q = int over a cycle of int E . J dA dt (J/m).
#
# Geometry and mesh (wire_2D_tphi.vtu, built by ../files/formulation_meshes.py): a scaled model, R = 1 m in the square
#   [-2, 2]^2 m (O-grid, 1008 QUAD8 cells). Wire TYPE / MATERIAL 2; cut strip TYPE 3 (MATERIAL 1); air TYPE 1.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 1: TPHI2D_edge / NON_CONDUCTING  / 1 / 0   air
#   * 2: TPHI2D_edge / SOLID_CONDUCTOR / 1 / 0   the superconductor (rho replaced by the power law)
#   * 3: TPHI2D_edge / CUT_LAYER       / 1 / 0   the cut
#   Materials: 1 air; 2 rho (replaced by the law). Real constants 1: thickness = 1 m (results per metre).
#   Nonlinear laws: power_law on MATERIAL 2 (set by the macro): Jc = 1e4 A/m^2, Ec = 1e-4 V/m, n = 30, rho0 = 1e-14.
#
# Boundary conditions: PHI = -I(t) theta / (2 pi) on the box; n x T = 0 on the wire surface, PHI = 0 inside the wire.
# Time history: backward Euler, period 40 s (E about Ec), 80 steps per cycle, 1.25 cycles (the loss from the first
#   peak to the next).
# Results: output/wire_2D_tphi_out_<n>.vtu, wire_2D_tphi_steps.csv, wire_2D_tphi.log; printed: the net current against
#   the imposed one, the Newton iterations, the loss per cycle against Norris's (ratio 1.22: the power law at n = 30
#   against the critical state, n -> infinity; 1.11 at n = 100).
#
# How to run:  cd examples/wire_2D_tphi && python macro.py          (about 40 s)
# Documentation: "Element library - TPHI2D_edge" and "Cuts", "User manual - Nonlinear laws (Superconductors: the
#   nonlinear settings)".
# =====================================================================================================================

import os

import numpy as np

from pytec3d import Solver, mesh, selection, tree_gauge

MU0 = 4e-7 * np.pi
RADIUS, JC, EC, N, RHO0 = 1.0, 1e4, 1e-4, 30.0, 1e-14
F, PERIOD, PER_CYCLE = 0.5, 40.0, 80
IC = JC * np.pi * RADIUS**2
DT = PERIOD / PER_CYCLE


def current(t):
    return F * IC * np.sin(2.0 * np.pi * t / PERIOD)


# 1. Solver, mesh, model data and the nonlinear settings. The boundary conditions are set by selections: select points
#    or cells, then set_value_on_selection() writes one value, or one value per selected id (in increasing id order).
pytec3d = Solver()
filename = "wire_2D_tphi"
pytec3d.read_model(os.path.join(os.getcwd(), filename + ".vtu"))
pytec3d.model.read_json("model_data.json")
grid = pytec3d.model.mesh
xyz = mesh.points(grid)
pytec3d.solver_settings = {"nonlinear": {"absolute_residual_tolerance": 1e-9, "relative_residual_tolerance": 1e-30,
                                         "max_iterations": 100, "require_residual_reduction": True}}

# 2. The power law as a table [J, E]: J = 0, then log-spaced from 1e-4 Jc (0.6 % apart: the piecewise-linear curve
#    follows the steep power law), up to E = 1e2 Ec, then a normal-state segment of resistivity 1e-6 Ohm m.
j = np.concatenate(([0.0], JC * np.geomspace(1e-4, 10.0, 2001)))
e = RHO0 * j + EC * (j / JC) ** N
keep = e < 1e2 * EC
j_max = JC * 1e2 ** (1.0 / N)
table = np.column_stack((np.concatenate((j[keep], [j_max, 10 * j_max])),
                         np.concatenate((e[keep], [1e2 * EC, 1e2 * EC + 1e-6 * 9 * j_max]))))
pytec3d.model.set_power_law(2, table)

# 3. Gauges of the wire: n x T = 0 on its surface edges (the midside points of the wire cells that also belong to the
#    other cells), PHI = 0 at its inner corner points.
selection.select(grid, "TYPE", 2, location="elem")
selection.nsle(grid)
selection.select(grid, "TYPE", [1, 3], location="elem")
selection.nsle(grid, mode="r")
selection.select(grid, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(grid, "_T", 0.0, 2)
selection.select(grid, "TYPE", 2, location="elem")
selection.nsle(grid)
selection.select(grid, "TYPE", [1, 3], location="elem")
selection.nsle(grid, mode="u")
selection.select(grid, "MidsideNodes", 0, mode="r")
selection.set_value_on_selection(grid, "_PHI", 0.0)

# 4. The cut: the edges of the strip (TYPE 3 cells), with their end points in the solver's orientation; crossing edges
#    have exactly one end on the ray y = 0, x > 0, the others _T = 0 (all the strip edges first, then the crossing
#    ones get the current in the loop: the crossing edge on the wire surface included).
on_ray = np.isclose(xyz[:, 1], 0.0) & (xyz[:, 0] > 0.0)
selection.select(grid, "TYPE", 3, location="elem")
mids, start, end = tree_gauge.edges(grid, selection.selected_ids(grid, "elem"))
selection.select(grid, "ID", mids)
selection.set_value_on_selection(grid, "_T", 0.0, 2)
cross = on_ray[start] != on_ray[end]
crossing = mids[cross]
sign = np.where(on_ray[start[cross]], 1.0, -1.0)          # +I from the ray side (c = 0) to the other side (c = 1)

# 5. The box: theta in (0, 2 pi], 2 pi on the ray.
selection.select_all(grid)
selection.nsel_ext(grid)
box = selection.selected_ids(grid, "node")
theta = np.mod(np.arctan2(xyz[box, 1], xyz[box, 0]), 2.0 * np.pi)
theta[on_ray[box]] = 2.0 * np.pi

# 6. Load the model; backward Euler; at each step the current on the crossing edges and the potential on the box
#    (one value per selected point).
pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 1.0
volumes = mesh.cell_sizes(grid)
wire = mesh.array_values(grid, "TYPE", "elem") == 2
loss = []
steps = PER_CYCLE + PER_CYCLE // 4
for step in range(1, steps + 1):
    t = step * DT
    selection.select(grid, "ID", crossing)
    selection.set_value_on_selection(grid, "_T", sign * current(t), 2)
    selection.select(grid, "ID", box)
    selection.set_value_on_selection(grid, "_PHI", -current(t) * theta / (2.0 * np.pi))
    pytec3d.time_step["time"] = t
    pytec3d.solve()
    pytec3d.write_output()
    e_cells, j_cells = mesh.array_values(grid, "E", "elem"), mesh.array_values(grid, "Jt", "elem")
    loss.append(float(np.einsum("ij,ij->i", e_cells, j_cells) @ volumes) * DT)
    if step % 10 == 0:
        net = j_cells[wire, 2] @ volumes[wire]
        print(f"t = {t:5.1f} s: I = {net / 1e3:7.3f} kA (imposed {current(t) / 1e3:7.3f}), "
              f"Newton iterations {pytec3d.steps[-1]['newton_iterations']}")

q = sum(loss[PER_CYCLE // 4:PER_CYCLE // 4 + PER_CYCLE])
norris = MU0 * IC**2 / np.pi * ((1 - F) * np.log(1 - F) + (2 - F) * F / 2)
print(f"loss per cycle {q:.4e} J/m, Norris (critical state) {norris:.4e} J/m, ratio {q / norris:.3f}")

pytec3d.save_model()

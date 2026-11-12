# =====================================================================================================================
# bulk_2D_tphi_axisym - axisymmetric H-formulation (TPHI2D_edge): field-cooled magnetization of a REBCO bulk,
#                       Jc(B) (law power_law_b), trapped field and flux creep
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A superconducting bulk disk (radius 10 mm, thickness 10 mm) is cooled in an axial field of 2 T, then the field is
# ramped down to zero: the bulk traps the field. TPHI2D_edge with option3 = 1 (axisymmetric about y: (x, y) = (r, z),
# the integrals times 2 pi r, component 2 of Jt and E the azimuthal one).
#   * no cut: the section of the bulk touches the axis, so the air of the half plane is simply connected (a ring
#     section away from the axis would need one);
#   * a critical current depending on the field, the Kim model Jc(B) = Jc0 / (1 + |B| / B0): the law power_law_b,
#     E = E0(J / s(|B|)), the power law E0 at B = 0 (a table) and s = Jc(B) / Jc0 (a second table), |B| = mu0 |H|
#     (set_power_law(..., jc_b=...));
#   * a STATIC step first, the field-cooled state: H = -grad(PHI) everywhere, B uniform, no current. In a static step
#     the model sets T = 0 on every conductor edge and frees PHI inside the conductor; the transient steps restore the
#     gauges (n x T = 0 on the surface, PHI = 0 inside);
#   * the applied field on the box, PHI = -Ha(t) z (the box is 6 radii away: the field of the bulk there is small);
#   * nonlinear settings for a power law: an absolute residual tolerance (relative off), and here NO
#     require_residual_reduction: from the current-free state the first full Newton step overshoots (the tangent at
#     J = 0 is the tiny floor rho0: the residual grows 1e4-1e5 times), and a line search that must reduce the
#     residual stalls on that first step (an "Unrecoverable failure to converge"); the full steps converge, in 15-18
#     iterations at the start of the ramp, 3-9 later (the Jacobian of power_law_b has no d/dB term);
#   * flux creep after the ramp: with a power law the currents relax (E > 0 at J near Jc), the trapped field decays
#     slowly, logarithmically in time.
#
# Geometry and mesh (bulk_2D_tphi_axisym.vtu, built by ../files/formulation_meshes.py): the half plane 0 < r < 60 mm,
#   |z| < 60 mm, 2016 QUAD8 cells; the bulk TYPE / MATERIAL 2, r < 10 mm, |z| < 5 mm, 0.5 mm cells; air TYPE 1.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 1: TPHI2D_edge / NON_CONDUCTING  / 1 / 1   air
#   * 2: TPHI2D_edge / SOLID_CONDUCTOR / 1 / 1   the bulk (rho replaced by the law)
#   Materials: 1 air; 2 rho (replaced by the law).
#   Nonlinear laws: power_law_b on MATERIAL 2 (set by the macro): Jc0 = 3e8 A/m^2, B0 = 0.5 T, Ec = 1e-4 V/m,
#     n = 25, floor rho0 = 1e-16 Ohm m (its diffusion time mu0 R^2 / rho0 is far longer than the run).
#
# Time history: step 1 static at t = 1 s, Ha = 2 T / mu0; then backward Euler, Ha ramped to 0 from t = 1 s to 101 s
#   (2 s steps), then held for 200 s (steps growing from 2 to 40 s).
# Results: output/bulk_2D_tphi_axisym_out_<n>.vtu, bulk_2D_tphi_axisym_steps.csv, bulk_2D_tphi_axisym.log; printed:
#   the applied field, B_z 1 mm above the centre of the top face (a Hall probe), the largest |J| / Jc0, Newton
#   iterations.
#
#   Trapped field at the probe: 0.55 T at the end of the ramp, 0.48 T 160 s later (creep).
#
# How to run:  cd examples/bulk_2D_tphi_axisym && python macro.py          (about 30 s)
# Documentation: "Element library - TPHI2D_edge", "User manual - Nonlinear laws (power_law_b; Superconductors: the
#   nonlinear settings)".
# =====================================================================================================================

import os

import numpy as np

from pytec3d import Solver, mesh, selection

MU0 = 4e-7 * np.pi
JC0, B0, EC, N, RHO0 = 3e8, 0.5, 1e-4, 25.0, 1e-16
BA, T_START, T_RAMP = 2.0, 1.0, 100.0      # the static step at T_START (times are > 0), then the ramp down

# 1. Solver, mesh, model data, nonlinear settings. The boundary conditions are set by selections: select points or
#    cells, then set_value_on_selection() writes one value, or one value per selected id (in increasing id order).
pytec3d = Solver()
filename = "bulk_2D_tphi_axisym"
pytec3d.read_model(os.path.join(os.getcwd(), filename + ".vtu"))
pytec3d.model.read_json("model_data.json")
grid = pytec3d.model.mesh
xyz = mesh.points(grid)
pytec3d.solver_settings = {"nonlinear": {"absolute_residual_tolerance": 1e-10, "relative_residual_tolerance": 1e-30,
                                         "max_iterations": 100}}

# 2. The law power_law_b: E0(J) at B = 0 as a table [J, E] (as in wire_2D_tphi: log-spaced, ending with a
#    normal-state segment), and s(|B|) = 1 / (1 + |B| / B0) as a table [|B|, s].
j = np.concatenate(([0.0], JC0 * np.geomspace(1e-4, 10.0, 2001)))
e = RHO0 * j + EC * (j / JC0) ** N
keep = e < 1e2 * EC
j_max = JC0 * 1e2 ** (1.0 / N)
table = np.column_stack((np.concatenate((j[keep], [j_max, 10 * j_max])),
                         np.concatenate((e[keep], [1e2 * EC, 1e2 * EC + 1e-6 * 9 * j_max]))))
b = np.concatenate(([0.0], np.geomspace(1e-4, 100.0, 1000)))
pytec3d.model.set_power_law(2, table, jc_b=np.column_stack((b, 1.0 / (1.0 + b / B0))))

# 3. The point sets of the gauges, as id lists (selected_ids), reselected at each step: the midside points of all the
#    bulk edges, those of its surface edges (points of the bulk cells also in air cells), its inner corner points.
selection.select(grid, "TYPE", 2, location="elem")
selection.nsle(grid)
selection.select(grid, "MidsideNodes", 1, mode="r")
bulk_edges = selection.selected_ids(grid)
selection.select(grid, "TYPE", 1, location="elem")
selection.nsle(grid, mode="r")
surface_edges = selection.selected_ids(grid)
selection.select(grid, "TYPE", 2, location="elem")
selection.nsle(grid)
selection.select(grid, "TYPE", 1, location="elem")
selection.nsle(grid, mode="u")
selection.select(grid, "MidsideNodes", 0, mode="r")
inner = selection.selected_ids(grid)


def gauges(static):
    """Static step: T = 0 on every edge of the bulk, PHI free inside. Transient: T free inside, n x T = 0 on the
    surface, PHI = 0 inside."""
    selection.select(grid, "ID", bulk_edges)
    selection.set_value_on_selection(grid, "_T", 0.0 if static else np.nan, 2)
    selection.select(grid, "ID", surface_edges)
    selection.set_value_on_selection(grid, "_T", 0.0, 2)
    selection.select(grid, "ID", inner)
    selection.set_value_on_selection(grid, "_PHI", np.nan if static else 0.0)


# 4. The box: the outer boundary r = 60 mm, z = +-60 mm, without the axis r = 0 (a natural boundary of PHI).
selection.select_all(grid)
selection.nsel_ext(grid)
selection.select_pos(grid, 0, (1e-9, 1.0), mode="r")
box = selection.selected_ids(grid)

# 5. The Hall probe 1 mm above the centre of the top face: the nearest node.
probe = int(np.argmin(np.linalg.norm(xyz - [0.0, 6e-3, 0.0], axis=1)))
bulk = mesh.array_values(grid, "TYPE", "elem") == 2

pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 1.0


def solve(t, ba, static=False):
    gauges(static)
    selection.select(grid, "ID", box)
    selection.set_value_on_selection(grid, "_PHI", -ba / MU0 * xyz[box, 1])
    pytec3d.time_step["time"] = t
    pytec3d.time_step["transient_effects"] = not static
    pytec3d.solve()
    pytec3d.write_output()
    jt = mesh.array_values(grid, "Jt", "elem")[bulk, 2]
    print(f"t = {t:6.1f} s: Ba = {ba:5.2f} T, B_z(probe) = {mesh.array_values(grid, 'B', 'node')[probe, 1]:6.3f} T, "
          f"max |J| / Jc0 = {np.abs(jt).max() / JC0:5.3f}, Newton iterations {pytec3d.steps[-1]['newton_iterations']}")


# 6. Field cooling: the static step at Ba, then the ramp down and the hold.
solve(T_START, BA, static=True)
times = np.concatenate((np.arange(2.0, T_RAMP + 1e-9, 2.0), T_RAMP + np.cumsum(np.geomspace(2.0, 40.0, 12))))
for t in times:
    solve(T_START + float(t), BA * max(1.0 - t / T_RAMP, 0.0))

pytec3d.save_model()

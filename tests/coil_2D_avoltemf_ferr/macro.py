# =====================================================================================================================
# coil_2D_avoltemf_ferr - coil_2D_avoltemf with a saturating iron frame (nonlinear B-H curve, Newton-Raphson)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The stranded-coil analysis of coil_2D_avoltemf (read that macro first) with the iron frame of coil_2D_avolt_ferr
# around the two rings. The coil inductance now depends on the saturation of the iron, and every step is solved
# with Newton's method. It also uses backward Euler (theta = 1), the most damped choice of the theta method.
#
# Geometry and mesh (coil_2D_avoltemf_ferr.vtu): axisymmetric about y, 800 QUAD8 cells. Coil TYPE / MATERIAL 23 at
#   r = 0.5 m, z = -0.5 m; passive ring TYPE / MATERIAL 24 at r = 0.5 m, z = +0.5 m. Iron frame TYPE / MATERIAL 26 set
#   by this macro: centroids inside 0.25 < r < 0.75, -0.75 < z < 0.75 and outside 0.4 < r < 0.6, -0.6 < z < 0.6.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV2D_nodal / NON_CONDUCTING  / 2 / 1   air
#   * 23: AV2D_nodal / STRANDED_COIL   / 2 / 1   coil: A_z + VOLT + EMF
#   * 26: AV2D_nodal / HARD_MAGNETIC   / 2 / 1   iron (M = 0 by default: soft iron), B-H curve below
#   * 24: AV2D_nodal / SOLID_CONDUCTOR / 2 / 1   passive ring (eddy currents only)
#   Materials: 21 air; 23 coil (mur = 1); 26 iron (mur = 1, replaced by the curve); 24 rho = 1e-8 Ohm m.
#   Real constants 1: thickness = 1, cross_section = 0.04 m^2, turns = 1000, radius = 0.5 m, resistance = 50 Ohm,
#     symmetry = 1, T = (0, 0, 1).
#   DOF couplings: 1 -> ["EMF", "VOLT"] (coil).
#   Nonlinear laws: bh_curve on MATERIAL 26, added by this macro.
#
# Boundary conditions:  A_z = 0 on the external boundary; step 2: coil VOLT = 10 kV.
# Loads and time history (theta = 1.0, dt = 1e-5 s): step 1 CURR = 100 A; step 2 CURR = 0, VOLT = 10 kV.
# Results: output/coil_2D_avoltemf_ferr_out_<n>.vtu, coil_2D_avoltemf_ferr_steps.csv, coil_2D_avoltemf_ferr.log
#
# How to run:  cd tests/coil_2D_avoltemf_ferr && python macro.py
# Documentation: "User manual - Nonlinear laws", "Theory - stranded coils", "Theory - Newton-Raphson".
# =====================================================================================================================


from pytec3d import mesh, selection, Solver
import os
import vtk
import vtk.util.numpy_support as vtknp

# 1. Solver and model files.
pytec3d = Solver()

filename = "coil_2D_avoltemf_ferr"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# Make sure the COORDINATE_SYSTEM cell array exists (no-op when present).
mesh.add_array(pytec3d.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

# 2. Iron frame: cells of the outer rectangle minus the inner window get TYPE and MATERIAL 26.
selection.select_all(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.4, 0.6), location="elem")                 # inner window, r
selection.select_pos(pytec3d.model.mesh, 1, (-0.6, 0.6), location="elem", mode="r")      # inner window, z
selection.inv(pytec3d.model.mesh)                                                         # all but the window
selection.select_pos(pytec3d.model.mesh, 0, (0.25, 0.75), location="elem", mode="r")     # within the outer rectangle
selection.select_pos(pytec3d.model.mesh, 1, (-0.75, 0.75), location="elem", mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "MATERIAL", 26, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "TYPE", 26, location="elem")

# 3. B-H curve of the iron, rows [H (A/m), B (T)]: selects the bh_curve law on MATERIAL 26.
pytec3d.model.set_bh_curve(26, [
    [0.0, 0.0],
    [100.0, 0.371974867],
    [200.0, 0.743949735],
    [255.0, 0.935935912],
    [299.086072, 1.102225479],
    [373.2104478, 1.286546042],
    [447.5127567, 1.415107144],
    [521.2510599, 1.487784841],
    [595.2477801, 1.552655102],
    [743.2890135, 1.640740625],
    [1113.801547, 1.762982443],
    [1484.740333, 1.837739212],
    [1855.209656, 1.878623862],
    [2226.874464, 1.921533805],
    [3661.029294, 1.998438582],
    [4396.08343, 2.034920312],
    [5135.739786, 2.06316393],
    [5870.655683, 2.091731846],
    [7345.089878, 2.14026725],
    [8819.744418, 2.18452467],
    [10289.25037, 2.21995201],
    [11749.12601, 2.260008846],
    [13201.51859, 2.297295089],
    [14644.04232, 2.326161358],
    [16091.24444, 2.357816875],
    [17528.87526, 2.386338652],
    [18959.20013, 2.408834349],
    [20396.79471, 2.434974089],
])

# 4. Boundary condition: A_z = 0 on the external boundary.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 5. Coil terminals: one VOLT and one EMF for the winding (CP id 1); master node = lowest node id.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

# 6. Passive ring: a SOLID_CONDUCTOR (massive, no electric DOF in 2D): a shorted ring, nothing to set.

# 7. Load the model; backward Euler (theta = 1).
pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 1.0

# 8. Time loop (Newton iterations at each step): current driven, then voltage driven.
selection.select_all(pytec3d.model.mesh)
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

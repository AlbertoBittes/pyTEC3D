# =====================================================================================================================
# coil_2D_avolt_ferr - coil_2D_avolt with a saturating iron yoke (nonlinear B-H curve, Newton-Raphson)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The analysis of coil_2D_avolt (read that macro first) with a rectangular iron frame around the two rings, modelled
# with a nonlinear B-H curve. A model with a nonlinear law is solved with Newton's method at every time step. It shows:
#   * how to change the TYPE and MATERIAL of a region of the mesh from a geometric selection (select_pos, inv);
#   * how to select the bh_curve law on a material with Model.set_bh_curve() (rows [H, B]);
#   * the large current (30 kA) that drives the iron into saturation; the log reports, for each step, the Newton
#     iterations and the largest |B| reached in the iron, with a warning if it goes beyond the last point of the curve.
#
# Geometry and mesh (coil_2D_avolt_ferr.vtu)
# ------------------------------------------
#   As coil_2D_avolt: 2D axisymmetric about the y axis (x = r, y = z), 0 <= r <= 1 m, -1 <= z <= 1 m, 800 QUAD8 cells.
#   * Lower ring (TYPE / MATERIAL 23) at r = 0.5 m, z = -0.5 m; upper ring (TYPE / MATERIAL 24) at r = 0.5, z = +0.5.
#   * Iron frame (TYPE / MATERIAL 26, set by this macro): the cells whose centroid is inside
#     0.25 < r < 0.75, -0.75 < z < 0.75 but outside 0.4 < r < 0.6, -0.6 < z < 0.6 (a 0.15 m thick rectangular frame
#     enclosing both rings).
#
# Model data (model_data.json)
# ----------------------------
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV2D_nodal / NON_CONDUCTING            / 2 / 1   air
#   * 23: AV2D_nodal / CIRCUIT_COUPLED_CONDUCTOR / 2 / 1   active ring (A_z + terminal V)
#   * 24: AV2D_nodal / SOLID_CONDUCTOR           / 2 / 1   passive shorted ring (eddy currents)
#   * 26: AV2D_nodal / HARD_MAGNETIC             / 2 / 1   iron: magnetic only, no current; its magnetization M is
#                                                          not given, so it is zero and the domain behaves as soft iron
#   Materials: 21 air (mur = 1), 23 rho = 1e-7, 24 rho = 1e-8, 26 mur = 1 (replaced by the B-H curve; rho unused).
#   Real constants: 1: thickness = 1 (full 360 degree ring).   DOF couplings: 1 -> ["V"].
#   Nonlinear laws: added by this macro: bh_curve on MATERIAL 26. The law replaces the reluctivity nu = 1/(mu0 mur)
#   by the secant nu(|B|) = H(|B|)/|B| in the residual and by the differential reluctivity in the Jacobian.
#
# Boundary conditions:  A_z = 0 on the external boundary; step 2: V = 0.232874 V on the terminal of the lower ring.
# Loads and time history (theta = 0.5, dt = 1e-5 s):
#   * step 1: total current AMPS = 30 kA in the lower ring;
#   * step 2: AMPS = 0, terminal voltage V = 0.232874 V imposed.
# Results: output/coil_2D_avolt_ferr_out_<n>.vtu (A, V, B, H, Jt), coil_2D_avolt_ferr_steps.csv, .log
#
# How to run:  cd tests/coil_2D_avolt_ferr && python macro.py
# Documentation: "User manual - Nonlinear laws", "Theory - Newton-Raphson", "Examples" (2D axisymmetric coils).
# =====================================================================================================================


# 0. Imports: mesh gives low-level access to the arrays of the VTK grid.
from pytec3d import mesh, selection, Solver
import os
import vtk.util.numpy_support as vtknp
import vtk

# 1. Solver and model files.
pytec3d = Solver()

filename = "coil_2D_avolt_ferr"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# Make sure the cell array COORDINATE_SYSTEM exists (add_array does nothing when it is already there). 0 = the global
# Cartesian system for every cell.
mesh.add_array(pytec3d.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

# ---------------------------------------------------------------------------------------------------------------------
# 2. Build the iron frame from geometric selections of cells (by centroid, location="elem")
#    select_pos(mesh, axis, (low, high)) selects by coordinate; mode "r" reselects within the current selection and
#    inv() inverts it. The sequence: inner window -> invert -> keep the outer rectangle = the frame.
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.4, 0.6), location="elem")                 # inner window, r
selection.select_pos(pytec3d.model.mesh, 1, (-0.6, 0.6), location="elem", mode="r")      # inner window, z
selection.inv(pytec3d.model.mesh)                                                         # everything but the window
selection.select_pos(pytec3d.model.mesh, 0, (0.25, 0.75), location="elem", mode="r")     # within the outer rectangle
selection.select_pos(pytec3d.model.mesh, 1, (-0.75, 0.75), location="elem", mode="r")
# The frame cells take the element type 26 (HARD_MAGNETIC) and the material 26 of the model data.
selection.set_value_on_selection(pytec3d.model.mesh, "MATERIAL", 26, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "TYPE", 26, location="elem")

# ---------------------------------------------------------------------------------------------------------------------
# 3. Nonlinear material: B-H curve of MATERIAL 26
#    set_bh_curve(material_id, rows) takes rows [H (A/m), B (T)], B strictly increasing, and writes
#    model.data["nonlinear_laws"]["26"] = {"type": "bh_curve", "B": [...], "H": [...]}. The curve is piecewise linear
#    and is extended beyond its last point along its last segment: give it up to deep saturation.
# ---------------------------------------------------------------------------------------------------------------------
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


# 4. Boundary condition: A_z = 0 on the external boundary (axis and far boundary).
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 5. Terminal of the circuit-coupled conductor: one V DOF for the whole ring (CP id 1, dof_couplings["1"] = ["V"]).
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

# 6. Load the model: validation, solver session; the model data (with the B-H curve) is now fixed.
pytec3d.load_model_into_solver()

# 7. Theta method, theta = 0.5.
selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 0.5

# ---------------------------------------------------------------------------------------------------------------------
# 8. Time loop: Newton iterations at every step (solver_settings["nonlinear"] holds the tolerances and the line search)
# ---------------------------------------------------------------------------------------------------------------------
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 1e-5
    if i == 1:
        # 30 kA through the lower ring: saturates the iron.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "AMPS", 30000)
    else:
        # Voltage driven: AMPS removed, terminal voltage imposed.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "AMPS", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0.232874)
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

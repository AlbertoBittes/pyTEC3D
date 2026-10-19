# =====================================================================================================================
# coil_3D_avoltemf_ferr - coil_3D_avoltemf with a saturating iron frame (nonlinear B-H curve), nodal A-V elements
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The 3D cyclic-sector stranded-coil model of coil_3D_avoltemf (read that macro first) with a toroidal iron frame
# around the two rings, the 3D counterpart of coil_2D_avoltemf_ferr. It shows:
#   * a geometric selection of cells in CYLINDRICAL coordinates (r and axial coordinate) to define the iron region;
#   * a B-H curve (bh_curve law) on a 3D nodal element type: the law replaces both the reluctivity nu and the
#     coefficient nu_gauge of the Coulomb gauge term, and the model is solved with Newton's method;
#   * very long time steps (dt = 1e5 s, far beyond the L/R time constants): each step is practically the DC steady
#     state of the coil at the imposed current or voltage.
#
# Geometry and mesh (coil_3D_avoltemf_ferr.vtu): 30 degree sector about y, 0 <= r <= 1 m, -1 <= y <= 1 m,
#   3800 HEX8 + 200 WEDGE6 cells. Coil TYPE / MATERIAL 23 at r = 0.5 m, y = -0.5 m; passive ring TYPE / MATERIAL 24 at
#   r = 0.5 m, y = +0.5 m. Iron frame TYPE / MATERIAL 26 set by this macro: centroids with 0.25 < r < 0.75 and
#   -0.75 < y < 0.75, outside 0.4 < r < 0.6, -0.6 < y < 0.6.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV3D_nodal / NON_CONDUCTING            / 1 / 0   air
#   * 23: AV3D_nodal / STRANDED_COIL             / 1 / 0   coil: A + VOLT + EMF
#   * 24: AV3D_nodal / SOLID_CONDUCTOR           / 1 / 0   passive ring: A + V
#   * 26: AV3D_nodal / HARD_MAGNETIC             / 1 / 0   iron (M = 0 by default: soft iron) with the B-H curve
#   Materials: 21 air; 23 coil; 26 iron (mur = 1, replaced by the curve); 24 rho = 1.09e-7 Ohm m.
#   Real constants 1: cross_section = 0.04 m^2, turns = 2000, resistance = 20 Ohm, volume = 0.12566 m^3 (whole coil),
#     symmetry = 12, T = (0, 1, 0) in the cylindrical element coordinate system 1 (toroidal).
#   DOF couplings: 1 -> ["VOLT", "EMF"]; periodic pairs (A, V) added by the macro in CP_1.
#   Nonlinear laws: bh_curve on MATERIAL 26 (added by the macro).
#
# Boundary conditions:  far field A = 0 (y = +-1, r = 1); axis A_x = A_z = 0; cyclic symmetry theta = 0 <-> 30 deg
#   for A and V (nodal frames of the 30 degree face rotated); V = 0 on node 2220 (on the passive ring).
# Loads and time history (theta = 1, dt = 1e5 s): step 1 coil current CURR = 15 A; step 2 CURR = 0, VOLT = 10 kV.
# Results: output/coil_3D_avoltemf_ferr_out_<n>.vtu, coil_3D_avoltemf_ferr_steps.csv, coil_3D_avoltemf_ferr.log
#
# How to run:  cd tests/coil_3D_avoltemf_ferr && python macro.py
# Documentation: "User manual - Nonlinear laws", "Theory - Newton-Raphson", "Element library" (AV3D_nodal).
# =====================================================================================================================


from pytec3d import mesh, selection, Solver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import vtk
import numpy as np

# 1. Solver and model files.
pytec3d = Solver()

filename = "coil_3D_avoltemf_ferr"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# Far boundary, coil size, sector angle (degrees), selection tolerance.
size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

# Empty point array for the periodic couplings (CP_0 holds the coil terminals).
mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

# ---------------------------------------------------------------------------------------------------------------------
# 2. Iron frame, selected by cell centroid in the Y-axis cylindrical system (index 0 = r, index 2 = axial y):
#    inner window -> invert -> keep the outer rectangle; the frame cells get TYPE and MATERIAL 26.
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.4, 0.6), location="elem", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 2, (-0.6, 0.6), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.inv(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.25, 0.75), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 2, (-0.75, 0.75), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
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

# 4. Far field: A = 0 on y = +1, y = -1 and r = 1.
selection.select_pos(
    pytec3d.model.mesh,
    2,
    (size - tol, size + tol),
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.select_pos(
    pytec3d.model.mesh,
    2,
    (-size - tol, -size + tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.select_pos(
    pytec3d.model.mesh,
    0,
    (size - tol, size + tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)


# 5. Axis r = 0: A_x = A_z = 0.
selection.select_pos(
    pytec3d.model.mesh, 0, (-tol, tol), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)


# 6. Coupling for stranded coils: element coordinate system 1 (cylindrical) for T, and the terminal coupling CP_0 = 1.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]


# 7. Cyclic symmetry: couple theta = 0 with theta = 30 degrees for A and V, and rotate the nodal frames of the 30 degree
#    face by 30 degrees about y.
selection.select_pos(
    pytec3d.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A", "V"),
        )

selection.select_pos(
    pytec3d.model.mesh,
    1,
    (extr_angle - tol, extr_angle + tol),
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", extr_angle, 2)

# 8. Ground the potential of the passive ring (node 2220 is on the ring).
selection.select(pytec3d.model.mesh, "ID", 2220)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 9. Load the model; backward Euler.
pytec3d.load_model_into_solver()

selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 1.0

# 10. Time loop with dt = 1e5 s (quasi-static steps, Newton iterations each): 15 A, then 10 kV.
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 1e5
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 15)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 11. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

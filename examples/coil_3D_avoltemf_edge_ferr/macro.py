# =====================================================================================================================
# coil_3D_avoltemf_edge_ferr - edge-element coil sector with a saturating iron frame (B-H curve, Newton-Raphson)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The edge-element model of coil_3D_avoltemf_edge (read that macro first) with a toroidal iron frame. It combines the tree-cotree gauge of the Nedelec field with a nonlinear B-H law, solved with
# Newton's method at every step. The iron region is defined BEFORE the model data is read: the TYPE and MATERIAL cell
# arrays can be edited at any time before the model is loaded into the solver.
#
# Geometry and mesh (coil_3D_avoltemf_edge_ferr.vtu): 30 degree sector about y, 0 <= r <= 1 m, -1 <= y <= 1 m,
#   3800 HEX20 + 200 WEDGE15 cells. Coil TYPE / MATERIAL 23 at r = 0.5 m, y = -0.5 m; passive ring TYPE / MATERIAL 24
#   at r = 0.5 m, y = +0.5 m; iron frame TYPE / MATERIAL 26 (set by the macro): 0.25 < r < 0.75, -0.75 < y < 0.75
#   minus 0.4 < r < 0.6, -0.6 < y < 0.6.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV3D_edge / NON_CONDUCTING            / 2 / 0   air
#   * 23: AV3D_edge / STRANDED_COIL             / 2 / 0   coil: edge A + VOLT + EMF
#   * 24: AV3D_edge / SOLID_CONDUCTOR           / 2 / 0   passive ring: edge A + V
#   * 26: AV3D_edge / HARD_MAGNETIC             / 2 / 0   iron (M = 0: soft iron) with the B-H curve
#   Materials: 21 air; 23 coil; 26 iron (mur = 1, replaced by the curve); 24 rho = 1.09e-7 Ohm m.
#   Real constants 1 (set on every cell): cross_section = 0.04 m^2, turns = 2000, resistance = 20 Ohm,
#     volume = 0.12566 m^3 (whole coil), symmetry = 12, T = (0, 1, 0) in the cylindrical element system 1.
#   DOF couplings (set by the macro, add_dof_couplings): 1 -> ["VOLT", "EMF"]; periodic pairs (A, V) added in CP_1.  Nonlinear laws: bh_curve on MATERIAL 26.
#
# Boundary conditions:  tangential A = 0 on r = 1 and y = +-1; tree gauge; cyclic symmetry theta = 0 <-> 30 deg
#   (A, V); V = 0 on the master node of the passive ring.
# Loads and time history (theta = 1, dt = 1e-5 s): step 1 CURR = 15 A (VOLT free); step 2 CURR = 0, VOLT = 10 kV.
# Results: output/coil_3D_avoltemf_edge_ferr_out_<n>.vtu, coil_3D_avoltemf_edge_ferr_steps.csv, .log
#
# How to run:  cd examples/coil_3D_avoltemf_edge_ferr && python macro.py
# Documentation: "Element library" (AV3D_edge), "User manual - Nonlinear laws, gauging edge elements".
# =====================================================================================================================


from pytec3d import mesh, selection, Solver, tree_gauge
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and mesh.
pytec3d = Solver()

filename = "coil_3D_avoltemf_edge_ferr"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

# 2. Iron frame, by cell centroid in the Y-axis cylindrical system (index 0 = r, 2 = axial y): inner window,
#    inverted, restricted to the outer rectangle; the frame cells get TYPE and MATERIAL 26.
selection.select_all(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.4, 0.6), location="elem", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 2, (-0.6, 0.6), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.inv(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.25, 0.75), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 2, (-0.75, 0.75), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "MATERIAL", 26, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "TYPE", 26, location="elem")


# 3. Every cell uses the real constants row 1, then read the model data.
selection.select_all(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CONSTANT", 1, location="elem")
pytec3d.model.read_json("model_data.json")
# DOF couplings (set here, not in the model data): the nodes with CP_0 = 1 share one VOLT and one EMF value.
pytec3d.model.add_dof_couplings({1: ["VOLT", "EMF"]})

# Far boundary, coil size, sector angle (degrees), selection tolerance.
size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

# 4. B-H curve of the iron, rows [H (A/m), B (T)]: selects the bh_curve law on MATERIAL 26.
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


# 5. Empty point array for the periodic couplings.
mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

# 6. Element coordinate system of the coil (cylindrical about y): T = (0, 1, 0) is toroidal.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")


# 7. Coil terminals: VOLT and EMF coupled (CP_0 = 1); master = lowest node id.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

# 8. Cyclic symmetry: nodes of theta = 0 paired with their images at theta = 30 degrees, coupling A and V.
selection.select_pos(
    pytec3d.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A", "V"),
        )

# 9. Far field: tangential A = 0 (component 2 of _A, midside nodes) on r = 1 and y = +-1; then the tree-cotree
#    gauge over the edges of all the (selected) cells, grown from the selected far-field nodes.
selection.select_all(pytec3d.model.mesh)
selection.select_pos(
    pytec3d.model.mesh,
    0,
    (size-tol, size+tol),
    mode="s",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.select_pos(
    pytec3d.model.mesh,
    2,
    (size-tol, size+tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.select_pos(
    pytec3d.model.mesh,
    2,
    (-size-tol, -size+tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)
tree_gauge.tree_gauge(
    pytec3d.model.mesh,
    periodic_point_coupling="CP_1"
)

# 10. Ground the floating potential of the passive ring on its master node.
selection.select(pytec3d.model.mesh, "MATERIAL", 24, location="elem")
selection.nsle(pytec3d.model.mesh)
# Take the first selected node as master node (lowest ID)
master_node_passive = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
selection.select(pytec3d.model.mesh, "ID", master_node_passive)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0.0)

# 11. Load the model; backward Euler. save_model() keeps a copy of the prepared model (overwritten at the end).
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 1.0
pytec3d.save_model()
# 12. Time loop (Newton iterations at each step): CURR = 15 A with VOLT free, then VOLT = 10 kV.
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 15)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", np.nan)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 13. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

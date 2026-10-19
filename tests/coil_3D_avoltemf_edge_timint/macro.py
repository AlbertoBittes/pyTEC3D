# =====================================================================================================================
# coil_3D_avoltemf_edge_timint - edge-element coil sector with time-integrated coil unknowns (AV3D_edge_timint)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The edge-element coil sector of coil_3D_avoltemf_edge (read that macro first: the script is the same) with the
# element type AV3D_edge_timint, where every electric unknown is a time integral, so that all the electric terms are in
# the damping matrix C and the A-V system is symmetric:
#   * coil terminals: VOLT -> int_0^t VOLT dt [V s] and EMF -> int_0^t EMF dt = Psi, the linked flux [Wb];
#   * passive ring: V is the time-integrated potential (in AV3D_edge it is the potential itself).
# CURR is still the coil current [A]; _VOLT = 10000 in step 2 fixes the time integral of the terminal voltage
# (10000 V s at t = 2e-5 s), not the voltage. Compare with the nodal time-integrated case coil_3D_avoltemf_timint.
#
# Geometry and mesh (coil_3D_avoltemf_edge_timint.vtu): 30 degree sector about y, 0 <= r <= 1 m, -1 <= y <= 1 m, 3800 HEX20
#   + 200 WEDGE15 cells; coil TYPE 23 at r = 0.5 m, y = -0.5 m; passive ring TYPE 24 at r = 0.5 m, y = +0.5 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 5 = AV3D_edge_timint / NON_CONDUCTING            / 2 / 0   air: edge A
#   * 23: 5 = AV3D_edge_timint / STRANDED_COIL             / 2 / 0   coil: edge A + int VOLT + int EMF
#   * 24: 5 = AV3D_edge_timint / SOLID_CONDUCTOR           / 2 / 0   passive ring: edge A + V (time-integrated)
#   Materials: 21 air; 23 coil; 24 rho = 2e-8 Ohm m (mur = 1).
#   Real constants 1 (set on every cell by the macro): cross_section = 0.04 m^2, turns = 2000, resistance = 20 Ohm,
#     volume = 0.12566 m^3 (whole coil), symmetry = 12, T = (0, 1, 0) in the cylindrical element system 1.
#   DOF couplings: 1 -> ["VOLT", "EMF"]; periodic pairs added by the macro in CP_1.
#
# Boundary conditions:  tangential A = 0 on r = 1 and y = +-1; tree gauge; cyclic symmetry
#   theta = 0 <-> 30 deg for A and V; V = 0 on the master node of the passive ring.
# Loads and time history (theta = 1, dt = 1e-5 s): step 1 CURR = 100 A (VOLT free); step 2 CURR = 0, _VOLT = 10000 (V s).
# Results: output/coil_3D_avoltemf_edge_timint_out_<n>.vtu, coil_3D_avoltemf_edge_timint_steps.csv, coil_3D_avoltemf_edge_timint.log
#
# How to run:  cd tests/coil_3D_avoltemf_edge_timint && python macro.py
# Documentation: "Element library" (AV3D_edge_timint), "User manual - gauging edge elements", "Examples" (3D coils).
# =====================================================================================================================


from pytec3d import mesh, selection, Solver, tree_gauge
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and mesh.
pytec3d = Solver()

filename = "coil_3D_avoltemf_edge_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)
# Every cell uses the real constants row 1 (the CONSTANT cell array selects the row of data["constants"]).
selection.select_all(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CONSTANT", 1, location="elem")
pytec3d.model.read_json("model_data.json")

# Far boundary, coil size, sector angle (degrees), selection tolerance.
size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

# Empty point array for the periodic couplings.
mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

# 2. Element coordinate system of the coil: cylindrical about y, so that T = (0, 1, 0) is toroidal.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")


# 3. Coil terminals: VOLT and EMF of the winding coupled (CP_0 = 1); master = lowest node id.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

# 4. Cyclic symmetry: the nodes of the face theta = 0 (corner and midside nodes) are paired with their images on
#    theta = 30 degrees; the pairs couple A (the edge values, through the midside nodes) and V.
selection.select_pos(
    pytec3d.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A", "V"),
        )

# ---------------------------------------------------------------------------------------------------------------------
# 5. Far field and gauge
#    a) Tangential A = 0 on the faces r = 1, y = +1 and y = -1: component 2 of _A on the selected nodes (only the
#       midside nodes carry edge DOFs; the value on the corner nodes is not used by the edge field).
#    b) tree_gauge(): builds a spanning tree of the edges of the selected cells (all of them: select_all() selected
#       every cell), grown from the selected nodes (the far-field faces), and sets _A = 0 on the midside nodes of the
#       tree edges. The periodic slave edges of CP_1 are left to the periodic constraint.
# ---------------------------------------------------------------------------------------------------------------------
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

# 6. Ground the floating potential of the passive ring: V = 0 on its master (lowest id) node.
selection.select(pytec3d.model.mesh, "MATERIAL", 24, location="elem")
selection.nsle(pytec3d.model.mesh)
# Take the first selected node as master node (lowest ID)
master_node_passive = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
selection.select(pytec3d.model.mesh, "ID", master_node_passive)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0.0)

# 7. Load the model; backward Euler. save_model() here keeps a copy of the prepared model (overwritten at the end).
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 1.0
pytec3d.save_model()

# 8. Time loop: CURR = 100 A with VOLT explicitly free (NaN), then VOLT = 10 kV.
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", np.nan)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

# =====================================================================================================================
# coil_3D_avoltemf_timint - coil_3D_avoltemf with time-integrated coil unknowns (AV3D_nodal_timint)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The analysis of coil_3D_avoltemf (read that macro first: the script is the same) with the element type
# AV3D_nodal_timint. In this formulation every electric unknown is a time integral, so that all the electric terms
# are in the damping matrix C and the coupled system is symmetric:
#   * coil terminals: VOLT -> int_0^t VOLT dt [V s] and EMF -> int_0^t EMF dt = Psi, the linked flux [Wb];
#   * passive ring: V is the time-integrated scalar potential (in AV3D_nodal it is the potential itself).
# Consequences for the input: CURR is still the coil current [A]; _VOLT = 10000 in step 2 fixes the time integral of
# the terminal voltage (10000 V s at t = 2e-5 s), not the voltage. On a static step the time-integrated DOFs, which
# have no stiffness row, are fixed to zero.
#
# Geometry and mesh (coil_3D_avoltemf_timint.vtu): as coil_3D_avoltemf: 30 degree sector about the y axis,
#   0 <= r <= 1 m, -1 <= y <= 1 m, 3800 HEX8 + 200 WEDGE6 cells; coil TYPE 23 at r = 0.5 m, y = -0.5 m; passive ring
#   TYPE 24 at r = 0.5 m, y = +0.5 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 4 = AV3D_nodal_timint / NON_CONDUCTING            / 1 / 0   air
#   * 23: 4 = AV3D_nodal_timint / STRANDED_COIL             / 1 / 0   coil: A + int VOLT + int EMF
#   * 24: 4 = AV3D_nodal_timint / SOLID_CONDUCTOR           / 1 / 0   passive ring: A + V (time-integrated)
#   Materials: 21 air; 23 coil; 24 rho = 2e-8 Ohm m (mur = 1).
#   Real constants 1: cross_section = 0.04 m^2, turns = 2000, resistance = 20 Ohm, volume = 0.12566 m^3 (whole coil),
#     symmetry = 12, T = (0, 1, 0) in the cylindrical element coordinate system 1 (toroidal).
#   DOF couplings: 1 -> ["VOLT", "EMF"]; periodic pairs (A, V) added by the macro in CP_1.
#
# Boundary conditions:  far field A = 0 (y = +-1, r = 1); axis A_x = A_z = 0; cyclic symmetry theta = 0 <-> 30 deg
#   (A, V; nodal frames of the 30 degree face rotated); V = 0 on node 2220 of the passive ring.
# Loads and time history (theta = 1, dt = 1e-5 s): step 1 CURR = 100 A; step 2 CURR = 0 and _VOLT = 10000 (V s).
# Results: output/coil_3D_avoltemf_timint_out_<n>.vtu, coil_3D_avoltemf_timint_steps.csv, .log
#
# How to run:  cd tests/coil_3D_avoltemf_timint && python macro.py
# Documentation: "Element library" (AV3D_nodal_timint), "Theory - time-integrated potentials, stranded coils".
# =====================================================================================================================


from pytec3d import mesh, selection, Solver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import vtk

# 1. Solver and model files.
pytec3d = Solver()

filename = "coil_3D_avoltemf_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# Geometric parameters of the model: far boundary (size), coil size, sector angle (degrees), selection tolerance.
size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

# A fresh, empty point array for the periodic couplings: CP_1 (CP_0 is used for the coil terminals). Several CP_k
# arrays let a node belong to more than one coupling.
mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

# ---------------------------------------------------------------------------------------------------------------------
# 2. Far field: A = 0 on the faces y = +1, y = -1 (axial coordinate, index 2 of the Y-axis cylindrical system) and on
#    the outer face r = 1 (index 0). mode="a" adds to the selection.
# ---------------------------------------------------------------------------------------------------------------------
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


# 3. Axis r = 0: the components of A normal to the axis (global x and z) vanish there.
selection.select_pos(
    pytec3d.model.mesh, 0, (-tol, tol), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)


# ---------------------------------------------------------------------------------------------------------------------
# 4. Coupling for stranded coils
#    a) Element coordinate system: the coil cells use the cylindrical system 1, so that the vector real constant
#       T = (0, 1, 0) is the toroidal direction e_theta in every cell (vector and tensor parameters are rotated from
#       the element system to the global one at the element centre).
#    b) Terminals: CP_0 = 1 on the coil nodes couples VOLT and EMF (dof_couplings["1"]); master = lowest node id.
# ---------------------------------------------------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------------------------------------------------
# 5. Cyclic symmetry
#    a) Select the nodes of the face theta = 0 (index 1 of the cylindrical system). The element selection that follows
#       does not change this NODE selection, which periodic_boundaries() uses.
#    b) periodic_boundaries(): each selected node is turned by dy = 30 degrees about the y axis and matched with the
#       mesh node there; both get a new CP id in CP_1, coupling A and V (added to dof_couplings).
#    c) node_angles: the nodal frames of the face theta = 30 degrees are rotated by 30 degrees about y, so that the
#       coupled components of A are the same cylindrical components (r, axial, theta) on the two faces.
# ---------------------------------------------------------------------------------------------------------------------
selection.select_pos(
    pytec3d.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
selection.select(pytec3d.model.mesh, "MATERIAL", 24, location="elem")
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

# 6. Ground the potential of the passive ring: V = 0 on one of its nodes (node 2220).
selection.select(pytec3d.model.mesh, "ID", 2220)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 7. Load the model into the solver session; backward Euler.
pytec3d.load_model_into_solver()

selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 1.0

# ---------------------------------------------------------------------------------------------------------------------
# 8. Time loop. CURR on the master node = current of the whole coil (symmetry and volume refer to the whole coil).
# ---------------------------------------------------------------------------------------------------------------------
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        # Current driven: I = 100 A.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        # Voltage driven: VOLT = 10 kV.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

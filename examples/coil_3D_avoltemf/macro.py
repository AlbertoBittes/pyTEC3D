# =====================================================================================================================
# coil_3D_avoltemf - 3D stranded coil (VOLT + EMF) and passive ring on a 30 degree cyclic sector, nodal A-V elements
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The 3D counterpart of coil_2D_avoltemf: the same two rings, revolved about the y axis, modelled on a 30 degree sector
# with cyclic-symmetry (periodic) conditions. The lower ring is a multi-turn stranded coil fed through its terminals
# (circuit unknowns VOLT and EMF); the upper ring is a passive massive conductor whose eddy currents are computed with
# the electric scalar potential V. It shows:
#   * selections in a cylindrical coordinate system (select_pos with coordinate_system=...);
#   * an ELEMENT coordinate system (COORDINATE_SYSTEM = 1, cylindrical about y) in which the winding direction T of
#     the coil is given: T = (0, 1, 0) = the toroidal direction everywhere in the coil;
#   * far-field, axis and cyclic-symmetry boundary conditions on a nodal vector potential, with the nodal frames of
#     one sector face rotated (node_angles) so that the coupled components match;
#   * Model.periodic_boundaries(): pairs of nodes of the two sector faces coupled for A and V (array CP_1);
#   * the grounding of the floating potential V of the passive conductor (one node with _V = 0);
#   * the symmetry real constant of a coil: the 30 degree sector is 1/12 of the coil (symmetry = 12), while the
#     loads (CURR, VOLT) refer to the whole coil.
#
# Geometry and mesh (coil_3D_avoltemf.vtu)
# ----------------------------------------
#   Global axes: y is the axis of revolution; in the Y-axis cylindrical system r = sqrt(x^2 + z^2), theta measured from
#   +x towards -z, axial coordinate = y. Sector 0 <= r <= 1 m, -1 <= y <= 1 m, 0 <= theta <= 30 degrees.
#   3800 HEX8 and 200 WEDGE6 (on the axis) first-order cells.
#   * Coil (TYPE / MATERIAL 23): 0.2 m x 0.2 m cross-section at r = 0.5 m, y = -0.5 m.
#   * Passive ring (TYPE / MATERIAL 24): same cross-section at r = 0.5 m, y = +0.5 m.
#   * Air (TYPE / MATERIAL 21).
#
# Model data (model_data.json)
# ----------------------------
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 1 = AV3D_nodal / NON_CONDUCTING            / 1 / 0   air: nodal A (3 components) + Coulomb gauge term
#   * 23: 1 = AV3D_nodal / STRANDED_COIL             / 1 / 0   coil: A + VOLT + EMF
#   * 24: 1 = AV3D_nodal / SOLID_CONDUCTOR           / 1 / 0   passive ring: A + V
#   Materials: 21 air (mur = 1); 23 coil (mur = 1, rho unused); 24 mur = 1, rho = 2e-8 Ohm m.
#   Real constants 1 (all cells; read by the coil):
#     cross_section S = 0.04 m^2, turns N = 2000, resistance R = 20 Ohm, volume = 0.12566 m^3 (volume of the WHOLE
#     coil, 2 pi 0.5 S), symmetry = 12 (whole coil / modelled part), T = (0, 1, 0) in the element coordinate system.
#   DOF couplings (set by the macro, add_dof_couplings): 1 -> ["VOLT", "EMF"] (coil terminals). The periodic pairs (A, V) are added by this macro.
#   Coordinate systems: 1 = Y-axis cylindrical, added by this macro and assigned to the coil cells.
#
# Boundary conditions
# -------------------
#   * Far field: A = 0 (all components) on the top and bottom faces (y = +-1 m) and on the outer face (r = 1 m).
#   * Axis (r = 0): A_x = A_z = 0 (the components normal to the axis).
#   * Cyclic symmetry: nodes of the face theta = 0 coupled to their images on theta = 30 degrees for A and V.
#   * V = 0 on one node of the passive ring (its potential is otherwise defined up to a constant).
#   * Step 2: coil terminal voltage VOLT = 10 kV.
#
# Loads and time history (theta = 1, backward Euler; dt = 1e-5 s)
# ---------------------------------------------------------------
#   * Step 1: coil current CURR = 100 A (whole coil).
#   * Step 2: CURR = 0, VOLT = 10 kV: the current follows from VOLT = R I + EMF.
#
# Results: output/coil_3D_avoltemf_out_<n>.vtu (A, V, VOLT, EMF, B, H, Jt), coil_3D_avoltemf_steps.csv, .log
#
# How to run:  cd examples/coil_3D_avoltemf && python macro.py
# Documentation: "Element library" (AV3D_nodal), "User manual - boundary conditions, DOF couplings, coordinate
# systems", "Theory - stranded coils".
# =====================================================================================================================


from pytec3d import mesh, selection, Solver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import vtk

# 1. Solver and model files.
pytec3d = Solver()

filename = "coil_3D_avoltemf"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")
# DOF couplings (set here, not in the model data): the nodes with CP_0 = 1 share one VOLT and one EMF value.
pytec3d.model.add_dof_couplings({1: ["VOLT", "EMF"]})

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
#    b) Terminals: CP_0 = 1 on the coil nodes couples VOLT and EMF (coupling 1, add_dof_couplings); master = lowest node id.
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

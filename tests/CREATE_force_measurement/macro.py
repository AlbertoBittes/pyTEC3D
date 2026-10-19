# =====================================================================================================================
# CREATE_force_measurement - static force measurement setup, edge elements
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A laboratory electromagnet (static): a racetrack STRANDED COIL around a conducting, saturating iron yoke, on a half
# model (symmetry = 2). It shows:
#   * a racetrack coil described piecewise by real constants rows: straight parts with a global winding direction T,
#     rounded ends with T = (0, 1, 0) in their own cylindrical ELEMENT coordinate systems (91, 92: Z-axis
#     cylindrical systems centred on the arc centres), so that T is the azimuthal direction of each arc;
#   * a nonlinear B-H curve on a conducting region (the yoke, SOLID_CONDUCTOR): Newton iterations together
#     with eddy currents;
#   * edge elements: tangential A = 0 on the external boundary and a tree-cotree gauge grown from it; the nodes of the plane y = 0
#     are left out of the tangential condition (symmetry plane where B is tangent);
#
# Geometry and mesh (CREATE_force_measurement.vtu): 0 <= x <= 0.2, 0 <= y <= 0.154, -0.1 <= z <= 0.177 m, 170040 HEX20 cells. Named cell arrays: AIR, COIL, COIL1 ... COIL5 (coil pieces), YOKE,
#   YOKE_LEFT, YOKE_RIGHT, YOKE_BRIDGE, SUPPORT_PLATE.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 1: 2 = AV3D_edge / NON_CONDUCTING            / 1 / 0   air and support plate
#   * 2: 2 = AV3D_edge / SOLID_CONDUCTOR           / 1 / 0   yoke: A + V, B-H curve
#   * 3: 2 = AV3D_edge / STRANDED_COIL             / 1 / 0   coil: A + VOLT + EMF
#   Materials: 2 yoke (rho = 1.09e-7 Ohm m, B-H curve added by the macro), 3 support plate (rho = 10, not used in a
#     non-conducting type), 10 coil, 11 air (mur = 1).
#   Real constants 1-5 (coil pieces): cross_section = 8.16e-5 m^2, turns = 127, volume = 1.27e-5 m^3 (whole coil), resistance = 1.144 Ohm, symmetry = 2, T per piece.
#   DOF couplings: 1 -> ["VOLT", "EMF"] (coil terminals).   Nonlinear laws: bh_curve on MATERIAL 2.
#
# Loads: one static step (t = 1e-5 s): coil current CURR = 5 A on the master node of the coil terminal (the current of
#   the whole coil, symmetry = 2).
# Results: output/CREATE_force_measurement_out_<n>.vtu (A, V, VOLT, EMF, B, H, Jt), CREATE_force_measurement_steps.csv, CREATE_force_measurement.log
#
# How to run:  cd tests/CREATE_force_measurement && python macro.py
# Documentation: "User manual - coordinate systems, nonlinear laws", "Element library" (stranded coil).
# =====================================================================================================================


from pytec3d import selection, Solver, tree_gauge
import os
from pytec3d.coordinate_system import CylindricalCoordinateSystem
from pytec3d.coordinate_system import ZAxisCylindricalCoordinateSystem
import vtk.util.numpy_support as vtknp


# from convert_mesh import convert_mesh_from_cdb

# 1. Solver and model files.
pytec3d = Solver()

# convert_mesh_from_cdb(
#     r"D:\pyTEC3D\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "CREATE_force_measurement"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# pytec3d.model.bh_curves[1] = [
#     [0.0, 0.0],
#     [98.3939940905, 0.0136155688067],
#     [146.889262617, 0.0213264689263],
#     [195.260662615, 0.0301054825108],
#     [243.695873428, 0.037431387609],
#     [291.607327062, 0.0486194287338],
#     [340.039423683, 0.0593478758073],
#     [388.464813606, 0.0708664329316],
#     [436.934832399, 0.0857570165318],
#     [485.208938096, 0.104183595294],
#     [581.617477627, 0.141325229405],
#     [673.458043004, 0.213527400286],
#     [774.711891476, 0.35653800222],
#     [865.993105257, 0.52966105117],
#     [967.219641459, 0.730254454256],
#     [1208.9308908, 1.12327379865],
#     [1450.97912972, 1.36148126568],
#     [1932.65944554, 1.61564633153],
#     [2411.44331821, 1.74225759844],
#     [2888.06786852, 1.82079053713],
#     [3359.66102603, 1.87224231279],
#     [3829.45853321, 1.91281144911],
#     [4293.81186332, 1.94450313269],
#     [4784.32897564, 1.97306058698],
#     [5671.19913662, 2.01774432657],
#     [6583.62023089, 2.05526569854],
#     [7434.39789268, 2.08888586741],
#     [8303.31592013, 2.11657067391],
# ]

# 2. B-H curve of the yoke (rows [H, B]): selects the bh_curve law on MATERIAL 2. The last point (B = 5 T) gives the
#    slope of the curve in deep saturation, where it is extended linearly.
pytec3d.model.set_bh_curve(2, [
    [0.0, 0.0],
    [100.0, 0.31],
    [200.0, 0.62],
    [255.0, 0.78],
    [299.086072, 0.918583873],
    [373.2104478, 1.072199867],
    [447.5127567, 1.17934968],
    [521.2510599, 1.239929871],
    [595.2477801, 1.29400392],
    [743.2890135, 1.367439528],
    [1113.801547, 1.46938531],
    [1484.740333, 1.531760307],
    [1855.209656, 1.565908439],
    [2226.874464, 1.601744566],
    [3661.029294, 1.666132249],
    [4396.08343, 1.69668764],
    [5135.739786, 1.720378902],
    [5870.655683, 1.744339419],
    [7345.089878, 1.785094394],
    [8819.744418, 1.822284428],
    [10289.25037, 1.852114984],
    [11749.12601, 1.885801436],
    [13201.51859, 1.917177494],
    [14644.04232, 1.941534839],
    [16091.24444, 1.968217538],
    [17528.87526, 1.992286782],
    [18959.20013, 2.011332763],
    [20396.79471, 2.033416969],
    [22468.45355, 2.055501175],
    [29881.25538, 2.099669587],
    [2337890.865, 5],
])

# Material 3 is linear: a curve would select the bh_curve law on it (nonlinear_laws).
# pytec3d.model.set_bh_curve(3, [
#     [0.0, 0.0],
#     [156.167232069, 0.239888044992],
#     [210.113464193, 0.316039547109],
#     [268.172738721, 0.376125353665],
#     [328.323296843, 0.411896921926],
#     [385.48042051, 0.430604288452],
#     [439.93479209, 0.441320807285],
#     [493.674322039, 0.447584189692],
#     [546.856093148, 0.45147893673],
#     [815.950268409, 0.460889474772],
#     [1081.56016787, 0.46628751316],
#     [1342.08642302, 0.470285285821],
#     [1601.49817015, 0.473868234982],
#     [1854.3064012, 0.477043124008],
#     [2199.37439592, 0.483047707152],
#     [2451.14947669, 0.485160621417],
#     [2684.76597187, 0.486756019625],
#     [3166.7575261, 0.489818986148],
#     [3637.99523776, 0.492120559561],
#     [4264.282677, 0.494422134],
#     [6131.085178, 0.499025282],
# ])

# Coupling for stranded coils
# 3. Coil terminals: VOLT and EMF of all the coil nodes coupled (CP_0 = 1); master = lowest node id.
selection.select(pytec3d.model.mesh, "TYPE", 3, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

# Specify coil rounds' ESYS
# 4. Element coordinate systems of the two rounded ends of the racetrack: cylindrical systems centred on the arc
#    centres; the cells of the real constants rows 2 and 4 use them, so that T = (0, 1, 0) is azimuthal there.
xc0, yc0, zc0 = 0.01786, 0.012965, 0.04555
xc1, yc1, zc1 = 0.01786, 0.050605, 0.04555

pytec3d.model.set_coordinate_system(91, ZAxisCylindricalCoordinateSystem(
    origin=[xc0, yc0, zc0],
))
pytec3d.model.set_coordinate_system(92, ZAxisCylindricalCoordinateSystem(
    origin=[xc1, yc1, zc1],
))

selection.select(pytec3d.model.mesh, "CONSTANT", 2, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 91, location="elem")

selection.select(pytec3d.model.mesh, "CONSTANT", 4, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 92, location="elem")

# Az: Symmetry BCs + Far field + Tree Gauging
# 5. Tangential A = 0 on the external boundary (component 2 of _A on the midside nodes), then the tree-cotree gauge
#    over the edges of all the cells, grown from the boundary nodes.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 1, (-1e-6, 1e-6), mode="u")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 2)
tree_gauge.tree_gauge(pytec3d.model.mesh)

# 6. Load the model; one static step.
pytec3d.load_model_into_solver()

pytec3d.time_step["time"] = 0.00001
pytec3d.time_step["transient_effects"] = False

# 7. Coil current CURR = 5 A on the master node of the coil terminal, solve, write the output.
selection.select(pytec3d.model.mesh, "ID", master_node)
selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 5)

pytec3d.solve()
pytec3d.write_output()

# 8. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

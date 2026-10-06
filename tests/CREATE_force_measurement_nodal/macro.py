from pytec3d_solver.solver import pyTEC3DSolver
import os
from pytec3d.coordinate_system import CylindricalCoordinateSystem
from pytec3d.tree_gauge import TreeGauge
from pytec3d.selector import Selector
from pytec3d.coordinate_system import ZAxisCylindricalCoordinateSystem
import vtk.util.numpy_support as vtknp

pytec3d = pyTEC3DSolver()

filename = "CREATE_force_measurement_nodal"

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

pytec3d.model.set_bh_curve(3, [
    [0.0, 0.0],
    [156.167232069, 0.239888044992],
    [210.113464193, 0.316039547109],
    [268.172738721, 0.376125353665],
    [328.323296843, 0.411896921926],
    [385.48042051, 0.430604288452],
    [439.93479209, 0.441320807285],
    [493.674322039, 0.447584189692],
    [546.856093148, 0.45147893673],
    [815.950268409, 0.460889474772],
    [1081.56016787, 0.46628751316],
    [1342.08642302, 0.470285285821],
    [1601.49817015, 0.473868234982],
    [1854.3064012, 0.477043124008],
    [2199.37439592, 0.483047707152],
    [2451.14947669, 0.485160621417],
    [2684.76597187, 0.486756019625],
    [3166.7575261, 0.489818986148],
    [3637.99523776, 0.492120559561],
    [4264.282677, 0.494422134],
    [6131.085178, 0.499025282],
])

# Coupling for stranded coils
Selector.select(pytec3d.model.mesh, "TYPE", 3, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

# Specify coil rounds' ESYS
xc0, yc0, zc0 = 0.01786, 0.012965, 0.04555
xc1, yc1, zc1 = 0.01786, 0.050605, 0.04555

pytec3d.model.set_coordinate_system(91, ZAxisCylindricalCoordinateSystem(
)
    origin=[xc0, yc0, zc0],
)
pytec3d.model.set_coordinate_system(92, ZAxisCylindricalCoordinateSystem(
)
    origin=[xc1, yc1, zc1],
)

Selector.select(pytec3d.model.mesh, "CONSTANT", 2, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 91, type="elem")

Selector.select(pytec3d.model.mesh, "CONSTANT", 4, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 92, type="elem")

# Az: Symmetry BCs + Far field + Tree Gauging
Selector.select_pos(pytec3d.model.mesh, 0, (-1e-6, 1e-6))
# Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select_pos(pytec3d.model.mesh, 1, (-1e-6, 1e-6))
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

Selector.select_all(pytec3d.model.mesh)
Selector.nsel_ext(pytec3d.model.mesh)
Selector.select_pos(pytec3d.model.mesh, 1, (-1e-6, 1e-6), mode="u")
Selector.select_pos(pytec3d.model.mesh, 0, (-1e-6, 1e-6), mode="u")
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.time_step["time"] = 0.00001
pytec3d.time_step["transient_effects"] = False

Selector.select(pytec3d.model.mesh, "ID", master_node)
Selector.set_value_on_selection(pytec3d.model.mesh, "AMPS", 5)

pytec3d.solve()
pytec3d.write_output()

pytec3d.save_model()

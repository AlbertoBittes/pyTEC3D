from apec_m import APEC_M
import os
import numpy as np
from apec_m.coordinate_system import CylindricalCoordinateSystem
from apec_m.tree_gauge import TreeGauge
from apec_m.selector import Selector
from apec_m.coordinate_system import XAxisCylindricalCoordinateSystem
import vtk.util.numpy_support as vtknp


# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

filename = "CREATE_B_measurement_nodal"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

# apec.model.bh_curves[1] = [
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

apec.model.bh_curves[2] = [
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
]

apec.model.bh_curves[3] = [
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
]

# Coupling for stranded coils
Selector.select(apec.model.mesh, "TYPE", 3, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CP_0", 1)
master_node = vtknp.vtk_to_numpy(
    apec.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

# Specify coil rounds' ESYS
xc0, yc0, zc0 = 0.0645, 0.035, 0.05
xc1, yc1, zc1 = 0.0645, 0.035, 0.078

apec.model.coordinate_systems[91] = XAxisCylindricalCoordinateSystem(
    origin=[xc0, yc0, zc0],
)
apec.model.coordinate_systems[92] = XAxisCylindricalCoordinateSystem(
    origin=[xc1, yc1, zc1],
)

Selector.select(apec.model.mesh, "CONSTANT", 2, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 91, type="elem")

Selector.select(apec.model.mesh, "CONSTANT", 4, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 92, type="elem")

Selector.select(apec.model.mesh, "CONSTANT", 2, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 91, type="elem")

Selector.select(apec.model.mesh, "CONSTANT", 4, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 92, type="elem")

# Az: Symmetry BCs + Far field
Selector.select_pos(apec.model.mesh, 1, (0.05-1e-6, 0.05+1e-6))
# Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select_all(apec.model.mesh)
Selector.nsel_ext(apec.model.mesh)
Selector.select_pos(apec.model.mesh, 1, (0.05-1e-6, 0.05+1e-6), mode="u")
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

# Symmetry BCs
Selector.select_all(apec.model.mesh)
Selector.nsel_ext(apec.model.mesh)
Selector.select(apec.model.mesh, "TYPE", 2, type="elem")
Selector.nsle(apec.model.mesh, mode="r")
Selector.set_value_on_selection(apec.model.mesh, "_V", 0.0)

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "times_CREATE.txt"))
voltages = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "volt_CREATE.txt"))

# # Consider only the time steps every 5 in the list
# time_steps = np.concatenate((time_steps[0:1], time_steps[1::2]))
# voltages = np.concatenate((voltages[0:1], voltages[1::2]))

apec.save_model()

apec.load_model_into_solver()

apec.model.time_step["theta"] = 1.0
for step_idx in range(len(time_steps)):
    apec.model.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.model.time_step["transient_effects"] = False
        apec.model.time_step["time"] = 1e-10
    else:
        apec.model.time_step["transient_effects"] = True

    Selector.select(apec.model.mesh, "ID", master_node)
    Selector.set_value_on_selection(apec.model.mesh, "_VOLT", voltages[step_idx])
    apec.solve()
    if step_idx % 10 == 0:
        apec.write_output()

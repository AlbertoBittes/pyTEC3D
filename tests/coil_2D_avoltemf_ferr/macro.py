from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
import os
import vtk
import vtk.util.numpy_support as vtknp

apec = APEC_M()

filename = "coil_2D_avoltemf"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

ModelEdit.add_array(apec.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

Selector.select_all(apec.model.mesh)
Selector.select_pos(apec.model.mesh, 0, (0.4, 0.6), type="elem")
Selector.select_pos(apec.model.mesh, 1, (-0.6, 0.6), type="elem", mode="r")
Selector.inv(apec.model.mesh)
Selector.select_pos(apec.model.mesh, 0, (0.25, 0.75), type="elem", mode="r")
Selector.select_pos(apec.model.mesh, 1, (-0.75, 0.75), type="elem", mode="r")
Selector.set_value_on_selection(apec.model.mesh, "MATERIAL", 26, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "TYPE", 26, type="elem")

apec.model.bh_curves[26] = [
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
]

Selector.select_all(apec.model.mesh)
Selector.nsel_ext(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(apec.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

Selector.select(apec.model.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CP_0", 2)
master_node_ground = vtknp.vtk_to_numpy(apec.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]
Selector.select(apec.model.mesh, "ID", master_node_ground)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

apec.load_model_into_solver()
apec.model.time_step["theta"] = 1.0

Selector.select_all(apec.model.mesh)
for i in range(1, 3):
    apec.model.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.model.mesh, "CURR", 100)
    else:
        Selector.select(apec.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.model.mesh, "CURR", 0)
        Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 10000)
    apec.solve()
    apec.write_output()

apec.save_model()

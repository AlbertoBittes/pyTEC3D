from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
from pytec3d.model_edit import ModelEdit
import os
import vtk
import vtk.util.numpy_support as vtknp

pytec3d = pyTEC3DSolver()

filename = "coil_2D_avoltemf_ferr"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_excel("input_template.xlsx")

ModelEdit.add_array(pytec3d.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

Selector.select_all(pytec3d.model.mesh)
Selector.select_pos(pytec3d.model.mesh, 0, (0.4, 0.6), type="elem")
Selector.select_pos(pytec3d.model.mesh, 1, (-0.6, 0.6), type="elem", mode="r")
Selector.inv(pytec3d.model.mesh)
Selector.select_pos(pytec3d.model.mesh, 0, (0.25, 0.75), type="elem", mode="r")
Selector.select_pos(pytec3d.model.mesh, 1, (-0.75, 0.75), type="elem", mode="r")
Selector.set_value_on_selection(pytec3d.model.mesh, "MATERIAL", 26, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "TYPE", 26, type="elem")

pytec3d.model.bh_curves[26] = [
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

Selector.select_all(pytec3d.model.mesh)
Selector.nsel_ext(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

Selector.select(pytec3d.model.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "CP_0", 2)
master_node_ground = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]
Selector.select(pytec3d.model.mesh, "ID", master_node_ground)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

pytec3d.load_model_into_solver()
pytec3d.model.time_step["theta"] = 1.0

Selector.select_all(pytec3d.model.mesh)
for i in range(1, 3):
    pytec3d.model.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(pytec3d.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        Selector.select(pytec3d.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

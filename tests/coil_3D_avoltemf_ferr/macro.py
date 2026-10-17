from pytec3d import mesh, selection, Solver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import vtk
import numpy as np

pytec3d = Solver()

filename = "coil_3D_avoltemf_ferr"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

selection.select_all(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.4, 0.6), location="elem", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 2, (-0.6, 0.6), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.inv(pytec3d.model.mesh)
selection.select_pos(pytec3d.model.mesh, 0, (0.25, 0.75), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 2, (-0.75, 0.75), location="elem", mode="r", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "MATERIAL", 26, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "TYPE", 26, location="elem")

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


selection.select_pos(
    pytec3d.model.mesh, 0, (-tol, tol), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)


# Coupling for stranded coils
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

selection.select(pytec3d.model.mesh, "ID", 2220)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

pytec3d.load_model_into_solver()

selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 1.0

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

pytec3d.save_model()

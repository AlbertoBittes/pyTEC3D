from pytec3d import mesh, selection, Solver, tree_gauge
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

pytec3d = Solver()

filename = "coil_3D_avoltemf_edge"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)
selection.select_all(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CONSTANT", 1, location="elem")
pytec3d.model.read_json("model_data.json")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")


selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
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

selection.select(pytec3d.model.mesh, "MATERIAL", 24, location="elem")
selection.nsle(pytec3d.model.mesh)
# Take the first selected node as master node (lowest ID)
master_node_passive = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
selection.select(pytec3d.model.mesh, "ID", master_node_passive)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0.0)

pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 1.0
pytec3d.save_model()
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

pytec3d.save_model()

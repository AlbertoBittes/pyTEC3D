from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
from pytec3d.model_edit import ModelEdit
import os
import vtk
import vtk.util.numpy_support as vtknp

pytec3d = pyTEC3DSolver()

filename = "coil_2D_avolt_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

ModelEdit.add_array(pytec3d.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

Selector.select_all(pytec3d.model.mesh)
Selector.nsel_ext(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

pytec3d.load_model_into_solver()

Selector.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 0.5

for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
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

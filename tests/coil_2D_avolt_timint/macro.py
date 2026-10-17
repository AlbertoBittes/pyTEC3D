from pytec3d import mesh, selection, Solver
import os
import vtk
import vtk.util.numpy_support as vtknp

pytec3d = Solver()

filename = "coil_2D_avolt_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

mesh.add_array(pytec3d.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

pytec3d.load_model_into_solver()

selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 0.5

for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

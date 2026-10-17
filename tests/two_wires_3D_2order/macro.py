from pytec3d import selection, Solver
import os

pytec3d = Solver()

filename = "two_wires_3D_2order"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

selection.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

selection.select(pytec3d.model.mesh, "TYPE", 24, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

selection.select(pytec3d.model.mesh, "TYPE", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1, 2, location="elem")

pytec3d.time_step["time"] = 1e-3
pytec3d.time_step["theta"] = 0.5
pytec3d.time_step["transient_effects"] = 1
pytec3d.save_model()
pytec3d.load_model_into_solver()

pytec3d.solve()
pytec3d.write_output()

pytec3d.save_model()

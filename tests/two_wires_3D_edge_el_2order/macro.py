from pytec3d import selection, Solver, tree_gauge
import os

pytec3d = Solver()

filename = "two_wires_3D_edge_el_2order"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

selection.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select(pytec3d.model.mesh, "LOWER_FACE_ROD", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)
selection.select(pytec3d.model.mesh, "UPPER_FACE_ROD", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

tree_gauge.tree_gauge(pytec3d.model.mesh)

selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1, 2, location="elem")

pytec3d.time_step["time"] = 1e-3
pytec3d.time_step["theta"] = 0.75
pytec3d.time_step["transient_effects"] = 0

pytec3d.load_model_into_solver()

pytec3d.solve()
pytec3d.write_output()

selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 100, 2, location="elem")

pytec3d.time_step["time"] = 2e-3
pytec3d.time_step["theta"] = 0.75
pytec3d.time_step["transient_effects"] = 1

pytec3d.solve()
pytec3d.write_output()

selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 2, location="elem")

pytec3d.time_step["time"] = 3e-3
pytec3d.time_step["theta"] = 0.65
pytec3d.time_step["transient_effects"] = 1

pytec3d.solve()
pytec3d.write_output()

pytec3d.save_model()

from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
from pytec3d.tree_gauge import TreeGauge
import os

pytec3d = pyTEC3DSolver()

filename = "two_wires_3D_edge_el_2order"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

Selector.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
Selector.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
Selector.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0, 0)

Selector.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
Selector.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
Selector.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0, 0)

Selector.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
Selector.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
Selector.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0, 0)

Selector.select(pytec3d.model.mesh, "LOWER_FACE_ROD", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)
Selector.select(pytec3d.model.mesh, "UPPER_FACE_ROD", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

TreeGauge.tree_gauge(pytec3d.model.mesh)

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", 1, 2, type="elem")

pytec3d.time_step["time"] = 1e-3
pytec3d.time_step["theta"] = 0.75
pytec3d.time_step["transient_effects"] = 0

pytec3d.load_model_into_solver()

pytec3d.solve()
pytec3d.write_output()

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", 100, 2, type="elem")

pytec3d.time_step["time"] = 2e-3
pytec3d.time_step["theta"] = 0.75
pytec3d.time_step["transient_effects"] = 1

pytec3d.solve()
pytec3d.write_output()

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 2, type="elem")

pytec3d.time_step["time"] = 3e-3
pytec3d.time_step["theta"] = 0.65
pytec3d.time_step["transient_effects"] = 1

pytec3d.solve()
pytec3d.write_output()

pytec3d.save_model()

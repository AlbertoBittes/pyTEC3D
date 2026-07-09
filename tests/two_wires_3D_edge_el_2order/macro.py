from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.tree_gauge import TreeGauge
import os

apec = APEC_M()

filename = "two_wires_3D_edge_el_2order"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.select_pos(apec.mesh, 0, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 0, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "_Az", 0, 0)

Selector.select_pos(apec.mesh, 1, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 1, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "_Az", 0, 0)

Selector.select_pos(apec.mesh, 2, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 2, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "_Az", 0, 0)

Selector.select(apec.mesh, "LOWER_FACE_ROD", 1)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)
Selector.select(apec.mesh, "UPPER_FACE_ROD", 1)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

TreeGauge.tree_gauge(apec.mesh)

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", 1, 2, type="elem")

apec.time_step["time"] = 1e-3
apec.time_step["theta"] = 0.75
apec.time_step["transient_effects"] = 0

apec.load_mesh_into_solver()

apec.solve()
apec.write_output()

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", 100, 2, type="elem")

apec.time_step["time"] = 2e-3
apec.time_step["theta"] = 0.75
apec.time_step["transient_effects"] = 1

apec.solve()
apec.write_output()

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", 0, 2, type="elem")

apec.time_step["time"] = 3e-3
apec.time_step["theta"] = 0.65
apec.time_step["transient_effects"] = 1

apec.solve()
apec.write_output()

apec.save_mesh()

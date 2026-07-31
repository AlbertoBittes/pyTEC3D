from apec_m import APEC_M
from apec_m.selector import Selector
import os

apec = APEC_M()

filename = "two_wires_3D_edge"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.select_pos(apec.mesh, 0, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 0, (-0.501, -0.49999), mode="a")
Selector.set_value_on_selection(apec.mesh, "_Az", 0, 0)

Selector.select_pos(apec.mesh, 1, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 1, (-0.501, -0.49999), mode="a")
Selector.set_value_on_selection(apec.mesh, "_Az", 0, 0)

Selector.select_pos(apec.mesh, 2, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 2, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "_Az", 0, 0)

Selector.select(apec.mesh, "TYPE", 24, type="elem")
Selector.nsle(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

Selector.select(apec.mesh, "TYPE", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", 1, 2, type="elem")

apec.time_step["time"] = 1e-3
apec.time_step["theta"] = 0.5
apec.time_step["transient_effects"] = 1
apec.save_mesh()
apec.load_mesh_into_solver()

apec.solve()
apec.write_output()

apec.save_mesh()

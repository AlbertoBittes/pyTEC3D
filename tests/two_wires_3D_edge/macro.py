from apec_m import APEC_M
from apec_m.selector import Selector
import os

apec = APEC_M()

filename = "two_wires_3D_edge"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

Selector.select_pos(apec.model.mesh, 0, (0.499, 0.50001))
Selector.select_pos(apec.model.mesh, 0, (-0.501, -0.49999), mode="a")
Selector.set_value_on_selection(apec.model.mesh, "_Az", 0, 0)

Selector.select_pos(apec.model.mesh, 1, (0.499, 0.50001))
Selector.select_pos(apec.model.mesh, 1, (-0.501, -0.49999), mode="a")
Selector.set_value_on_selection(apec.model.mesh, "_Az", 0, 0)

Selector.select_pos(apec.model.mesh, 2, (0.499, 0.50001))
Selector.select_pos(apec.model.mesh, 2, (-0.501, -0.49999), mode="a")
Selector.select(apec.model.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.model.mesh, "_Az", 0, 0)

Selector.select(apec.model.mesh, "TYPE", 24, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

Selector.select(apec.model.mesh, "TYPE", 23, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "Js", 1, 2, type="elem")

apec.model.time_step["time"] = 1e-3
apec.model.time_step["theta"] = 0.5
apec.model.time_step["transient_effects"] = 1
apec.save_model()
apec.load_model_into_solver()

apec.solve()
apec.write_output()

apec.save_model()

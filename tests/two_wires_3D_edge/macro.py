from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
import os

pytec3d = pyTEC3DSolver()

filename = "two_wires_3D_edge"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

Selector.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
Selector.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0, 0)

Selector.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
Selector.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0, 0)

Selector.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
Selector.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
Selector.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0, 0)

Selector.select(pytec3d.model.mesh, "TYPE", 24, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

Selector.select(pytec3d.model.mesh, "TYPE", 23, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", 1, 2, type="elem")

pytec3d.time_step["time"] = 1e-3
pytec3d.time_step["theta"] = 0.5
pytec3d.time_step["transient_effects"] = 1
pytec3d.save_model()
pytec3d.load_model_into_solver()

pytec3d.solve()
pytec3d.write_output()

pytec3d.save_model()

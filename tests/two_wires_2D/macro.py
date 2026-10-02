from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
import os

pytec3d = pyTEC3DSolver()

filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_excel("input_template.xlsx")

Selector.nsel_ext(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.model.time_step["theta"] = 0.91

for i in range(1, 3):
    if i == 1:
        pytec3d.model.time_step["transient_effects"] = 0
        pytec3d.model.time_step["time"] = i * 0.00001
    else:
        pytec3d.model.time_step["transient_effects"] = 1
        pytec3d.model.time_step["time"] = i * 0.00002

    Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i * 1000, 2, type="elem")
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

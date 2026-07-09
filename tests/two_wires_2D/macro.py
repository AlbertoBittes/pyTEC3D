from apec_m import APEC_M
from apec_m.selector import Selector
import os

apec = APEC_M()

filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

apec.load_mesh_into_solver()

apec.time_step["theta"] = 0.91

for i in range(1, 3):
    if i == 1:
        apec.time_step["transient_effects"] = 0
        apec.time_step["time"] = i * 0.00001
    else:
        apec.time_step["transient_effects"] = 1
        apec.time_step["time"] = i * 0.00002

    Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", i * 1000, 2, type="elem")
    apec.solve()
    apec.write_output()

apec.save_mesh()

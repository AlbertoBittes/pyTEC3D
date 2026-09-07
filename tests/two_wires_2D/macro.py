from apec_m_solver import APEC_M
from apec_m.selector import Selector
import os

apec = APEC_M()

filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

Selector.nsel_ext(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

apec.load_model_into_solver()

apec.model.time_step["theta"] = 0.91

for i in range(1, 3):
    if i == 1:
        apec.model.time_step["transient_effects"] = 0
        apec.model.time_step["time"] = i * 0.00001
    else:
        apec.model.time_step["transient_effects"] = 1
        apec.model.time_step["time"] = i * 0.00002

    Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
    Selector.set_value_on_selection(apec.model.mesh, "Js", i * 1000, 2, type="elem")
    apec.solve()
    apec.write_output()

apec.save_model()

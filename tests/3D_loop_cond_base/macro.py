from apec_m import APEC_M
from apec_m.selector import Selector
import os

apec = APEC_M()

filename = "3D_loop_cond_base"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

apec.load_mesh_into_solver()

Selector.select(apec.mesh, "Z_BOUNDS_1", 1)
Selector.select(apec.mesh, "Z_BOUNDS_2", 1, mode="a")
Selector.select(apec.mesh, "Y_BOUNDS_2", 1, mode="a")
Selector.select(apec.mesh, "X_BOUNDS_2", 1, mode="a")
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

Selector.select(apec.mesh, "X_BOUNDS_1", 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

Selector.select(apec.mesh, "Y_BOUNDS_1", 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

Selector.select(apec.mesh, "VOLT_0_1", 1)
Selector.select(apec.mesh, "VOLT_0_2", 1, mode="a")
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

b = 2**0.5 / 2

for i in range(1, 3):
    apec.time_step["time"] = i * 0.0001
    if i == 1:
        apec.time_step["transient_effects"] = 0
    else:
        apec.time_step["transient_effects"] = 1

    Selector.select(apec.mesh, "COIL_NOCOND_1", 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", 1 / i, 1, type="elem")
    Selector.select(apec.mesh, "COIL_NOCOND_2", 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", -1 / i, 0, type="elem")
    Selector.select(apec.mesh, "COIL_NOCOND_3", 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", -b / i, 0, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", b / i, 1, type="elem")
    apec.solve()
    apec.write_output()

apec.save_mesh()

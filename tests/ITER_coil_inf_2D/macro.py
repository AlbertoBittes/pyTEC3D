from apec_m import APEC_M
from apec_m.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_inf_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh, "s")
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

apec.load_mesh_into_solver()

apec.time_step["theta"] = 0.5

for i in range(1, 7):
    apec.time_step["time_step_number"] = i
    if i == 1:
        apec.time_step["transient_effects"] = 0
        apec.time_step["time"] = i * 0.0001
    elif i == 5:
        apec.time_step["transient_effects"] = 1
        apec.time_step["time"] = i * 0.0001
        apec.time_step["theta"] = 0.85
    elif i == 6:
        apec.time_step["transient_effects"] = 0
        apec.time_step["time"] = i * 0.001
    else:
        apec.time_step["transient_effects"] = 1
        apec.time_step["time"] = i * 0.0001

    Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", 1 / i, 2, type="elem")
    apec.solve()
    apec.write_output()

apec.save_mesh()

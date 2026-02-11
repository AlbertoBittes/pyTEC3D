from apec_m import APEC_M
from selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "3D_loop_cond_base"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

apec.write_model_data()

Selector.select(apec.mesh, "Z_BOUNDS_1", 1)
Selector.select(apec.mesh, "Z_BOUNDS_2", 1, mode="a")
Selector.select(apec.mesh, "Y_BOUNDS_2", 1, mode="a")
Selector.select(apec.mesh, "X_BOUNDS_2", 1, mode="a")
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "X_BOUNDS_1", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "Y_BOUNDS_1", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "VOLT_0_1", 1)
Selector.select(apec.mesh, "VOLT_0_2", 1, mode="a")
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

b = 2**0.5 / 2

apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 1.0
apec.time_parameters["write_substep_results"] = 1

for i in range(1, 3):
    apec.time_parameters["time_step_number"] = i
    apec.time_parameters["time"] = i * 0.0001
    if i == 1:
        apec.time_parameters["transient_effects"] = 0
        apec.time_parameters["number_of_substeps"] = 1
    else:
        apec.time_parameters["transient_effects"] = 1
        apec.time_parameters["number_of_substeps"] = 4

    Selector.select(apec.mesh, "COIL_NOCOND_1", 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", 1 / i, 1, type="elem")
    Selector.select(apec.mesh, "COIL_NOCOND_2", 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", -1 / i, 0, type="elem")
    Selector.select(apec.mesh, "COIL_NOCOND_3", 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", -b / i, 0, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", b / i, 1, type="elem")
    apec.write_load_step(i)

apec.save_mesh()

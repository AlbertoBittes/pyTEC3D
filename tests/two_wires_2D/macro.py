from apec_m import APEC_M
from selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

apec.write_model_data()

Selector.select(apec.mesh, "BOUNDARY_NODES", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

apec.time_parameters["ramped_loads"] = 0
apec.time_parameters["theta"] = 0.91
apec.time_parameters["write_substep_results"] = 1

for i in range(1, 3):
    apec.time_parameters["time_step_number"] = i
    if i == 1:
        apec.time_parameters["transient_effects"] = 0
        apec.time_parameters["number_of_substeps"] = 1
        apec.time_parameters["time"] = i * 0.00001
    else:
        apec.time_parameters["transient_effects"] = 1
        apec.time_parameters["number_of_substeps"] = 3
        apec.time_parameters["time"] = i * 0.00002

    Selector.select(apec.mesh, "ansys_material_type", 23, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", i * 1000, 2, type="elem")
    apec.write_load_step(i)

apec.save_mesh()

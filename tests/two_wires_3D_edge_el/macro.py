from apec_m import APEC_M
from apec_m.tree_gauge import TreeGauge
from apec_m.selector import Selector
import os
# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )

filename = "two_wires_3D_edge_el"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

apec.write_model_data()

Selector.select_pos(apec.mesh, 0, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 0, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "Az", 0, 0)

Selector.select_pos(apec.mesh, 1, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 1, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "Az", 0, 0)

Selector.select_pos(apec.mesh, 2, (0.499, 0.50001))
Selector.select_pos(apec.mesh, 2, (-0.501, -0.49999), mode="a")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "Az", 0, 0)

# Selector.nsel_ext(apec.mesh)
# Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
# Selector.set_value_on_selection(apec.mesh, "Az", 0, 0)

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", 1, 2, type="elem")
Selector.select(apec.mesh, "LOWER_FACE_ROD", 1)
Selector.select(apec.mesh, "MidsideNodes", 0, mode="r")
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)
Selector.select(apec.mesh, "UPPER_FACE_ROD", 1)
Selector.select(apec.mesh, "MidsideNodes", 0, mode="r")
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)
# Selector.select(apec.mesh, "UPPER_FACE_ROD", 1)
# Selector.set_value_on_selection(apec.mesh, "AMPS", 1)

TreeGauge.tree_gauge(apec.mesh)

apec.time_parameters["time_step_number"] = 1
apec.time_parameters["time"] = 1e-3
apec.time_parameters["theta"] = 0.5
apec.time_parameters["transient_effects"] = 1

apec.write_load_step(1)

apec.save_mesh()

from apec_m import APEC_M
from mesh_edit import MeshEdit
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

MeshEdit.reset_vtk_arrays_to_default((apec.mesh))

apec.read_excel("input_template.xlsx")

apec.write_model_data()

apec.select("BOUNDARY_NODES", 1)
apec.set_value_on_selection("A", 0, 2)

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

    apec.select("ansys_material_type", 23, type="elem")
    apec.set_value_on_selection("Js", i * 1000, 2, type="elem")
    apec.write_load_step(i)

apec.save_mesh()

from apec_m import APEC_M
from vtk_grid_helper import VTKGridHelper
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_inf_2D"

filepath = os.path.join(os.getcwd(),"tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

VTKGridHelper.reset_vtk_arrays_to_default((apec.mesh))

apec.read_excel("input_template.xlsx")

apec.write_model_data()

apec.select_pos(0, 0)
apec.set_value_on_selection("A", 0, 2)
apec.select("BASE_INF_NODES", 1)
apec.set_value_on_selection("A", 0, 2)

b = 2**0.5 / 2

apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 0.5
apec.time_parameters["write_substep_results"] = 1

for i in range(1, 7):
    apec.time_parameters["time_step_number"] = i
    if i == 1:
        apec.time_parameters["ramped_loads"] = 1
        apec.time_parameters["transient_effects"] = 0
        apec.time_parameters["number_of_substeps"] = 1
        apec.time_parameters["time"] = i * 0.0001
    elif i == 5:
        apec.time_parameters["ramped_loads"] = 0
        apec.time_parameters["transient_effects"] = 1
        apec.time_parameters["number_of_substeps"] = 3
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["theta"] = 0.85
    elif i == 6:
        apec.time_parameters["ramped_loads"] = 0
        apec.time_parameters["transient_effects"] = 0
        apec.time_parameters["number_of_substeps"] = 1
        apec.time_parameters["time"] = i * 0.001
    else:
        apec.time_parameters["ramped_loads"] = 1
        apec.time_parameters["transient_effects"] = 1
        apec.time_parameters["number_of_substeps"] = 4
        apec.time_parameters["time"] = i * 0.0001

    apec.select("ansys_material_type", 23, type="elem")
    apec.set_value_on_selection("Js", 1 / i, 2, type="elem")
    apec.write_load_step(i)

apec.save_mesh()

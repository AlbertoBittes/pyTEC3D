from apec_m import APEC_M
from vtk_grid_helper import VTKGridHelper
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "3D_loop_cond_base"

filepath = os.path.join(os.getcwd(),"tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

VTKGridHelper.reset_vtk_arrays_to_default((apec.mesh))

apec.read_excel("input_template.xlsx")

apec.write_model_data()

apec.select("Z_BOUNDS_1", 1)
apec.select("Z_BOUNDS_2", 1, mode="a")
apec.select("Y_BOUNDS_2", 1, mode="a")
apec.select("X_BOUNDS_2", 1, mode="a")
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 1)
apec.set_value_on_selection("A", 0, 2)

apec.select("X_BOUNDS_1", 1)
apec.set_value_on_selection("A", 0, 1)
apec.set_value_on_selection("A", 0, 2)

apec.select("Y_BOUNDS_1", 1)
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 2)

apec.select("VOLT_0_1", 1)
apec.select("VOLT_0_2", 1, mode="a")
apec.set_value_on_selection("VOLT", 0)

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

    apec.select("COIL_NOCOND_1", 1, cell_or_point="cell_data")
    apec.set_value_on_selection("Js", 1 / i, 1, cell_or_point="cell_data")
    apec.select("COIL_NOCOND_2", 1, cell_or_point="cell_data")
    apec.set_value_on_selection("Js", -1 / i, 0, cell_or_point="cell_data")
    apec.select("COIL_NOCOND_3", 1, cell_or_point="cell_data")
    apec.set_value_on_selection("Js", -b / i, 0, cell_or_point="cell_data")
    apec.set_value_on_selection("Js", b / i, 1, cell_or_point="cell_data")
    apec.write_load_step(i)

apec.save_mesh()

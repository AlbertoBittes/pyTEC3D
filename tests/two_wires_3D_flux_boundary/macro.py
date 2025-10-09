from apec_m import APEC_M
from vtk_grid_helper import VTKGridHelper

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )

apec.read_mesh(
    r"D:\PlasmaArc\simulations\two_wires_3D_flux_boundary\two_wires_3D_flux_boundary.vtu"
)

VTKGridHelper.reset_vtk_arrays_to_default((apec.mesh))

apec.read_excel("input_template.xlsx")

apec.select_pos(0, (0.499, 0.50001))
apec.select_pos(0, (-0.501, -0.49999), mode="a")
apec.set_value_on_selection("flux_normal", 1)
apec.set_value_on_selection("flux_parallel", 1)

apec.select_pos(1, (0.499, 0.50001))
apec.select_pos(1, (-0.501, -0.49999), mode="a")
apec.set_value_on_selection("flux_normal", 1)
apec.set_value_on_selection("flux_parallel", 1)

apec.select_pos(2, (0.499, 0.50001))
apec.select_pos(2, (-0.501, -0.49999), mode="a")
apec.set_value_on_selection("flux_parallel", 1)

# apec.nsel_ext()
# apec.select("MidsideNodes", 1, mode="r")
# apec.set_value_on_selection("Az", 0, 0)


apec.select("LOWER_FACE_ROD", 1)
apec.set_value_on_selection("VOLT", 0)
# apec.select("UPPER_FACE_ROD", 1)
# apec.set_value_on_selection("VOLT", 0)


apec.select("UPPER_FACE_ROD", 1)
apec.set_value_on_selection("AMPS", 1)
# apec.select("ansys_material_type", 23, cell_or_point="cell_data")
# apec.set_value_on_selection("Js", 1, 2, cell_or_point="cell_data")

apec.time_parameters["time_step_number"] = 1
apec.time_parameters["time"] = 1e-3
apec.time_parameters["number_of_substeps"] = 2
apec.time_parameters["theta"] = 0.55
apec.time_parameters["transient_effects"] = 1
apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["write_substep_results"] = 1

apec.write_load_step("load_steps\\load_step_1.bin")


apec.select("UPPER_FACE_ROD", 1)
apec.set_value_on_selection("AMPS", 100)
# apec.select("ansys_material_type", 23, cell_or_point="cell_data")
# apec.set_value_on_selection("Js", 100, 2, cell_or_point="cell_data")

apec.time_parameters["time_step_number"] = 1
apec.time_parameters["time"] = 2e-3
apec.time_parameters["number_of_substeps"] = 5
apec.time_parameters["theta"] = 0.85
apec.time_parameters["transient_effects"] = 1
apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["write_substep_results"] = 1

apec.write_load_step("load_steps\\load_step_2.bin")

apec.select("UPPER_FACE_ROD", 1)
apec.set_value_on_selection("AMPS", 0)
# apec.select("ansys_material_type", 23, cell_or_point="cell_data")
# apec.set_value_on_selection("Js", 0, 2, cell_or_point="cell_data")

apec.time_parameters["time_step_number"] = 1
apec.time_parameters["time"] = 3e-3
apec.time_parameters["number_of_substeps"] = 5
apec.time_parameters["theta"] = 0.91
apec.time_parameters["transient_effects"] = 1
apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["write_substep_results"] = 1

apec.write_load_step("load_steps\\load_step_3.bin")

apec.save_mesh()

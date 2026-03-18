from apec_m import APEC_M
from selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

aplasma = 0

filename = "ITER_2D_APEC_coils"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.apec_settings["target_central_axis_x"] = 6.52136102486433
apec.apec_settings["target_central_axis_y"] = 0.536537145137914
apec.apec_settings["target_central_axis_z"] = 0.0
apec.apec_settings["R0"] = 6.2
apec.apec_settings["max_picard_iterations"] = 70
apec.apec_settings["min_distance"] = 0.6
apec.apec_settings["cp_gain"] = 0.250
apec.apec_settings["tolerance_value"] = 1e-6
apec.apec_settings["bradx_value"] = -3.873415861e-08
apec.apec_settings["brady_value"] = 5.313857037e-09
apec.apec_settings["bverx_value"] = 1.822467692e-10
apec.apec_settings["bvery_value"] = 1.133287835e-07
apec.apec_settings["alpha_mp"] = 1.0001
apec.apec_settings["alpha_np"] = 0.83197
apec.apec_settings["plasma_current"] = 15e6
apec.apec_settings["beta_p"] = 6.8216e-001
apec.apec_settings["print_full_output"] = 0

i_CS3U = -4.8204e06
i_CS2U = -4.0480e06
i_CS1U = -2.0804e07
i_CS1L = -2.0804e07
i_CS2L = -1.0297e07
i_CS3L = 5.3012e06

i_PF1 = 5.7325e06
i_PF2 = -2.8834e06
i_PF3 = -5.8084e06
i_PF4 = -4.7762e06
i_PF5 = -7.8133e06
i_PF6 = 1.6941e07

a_PF1 = 0.94375190
a_PF2 = 0.41453946
a_PF3 = 0.66413094
a_PF4 = 0.60871516
a_PF5 = 0.77496250
a_PF6 = 1.72659250
a_CS = 1.55608
a_plasma = 28.6027095

apec.read_excel("input_template.xlsx")

apec.write_model_data()

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

apec.time_parameters["ramped_loads"] = 0
apec.time_parameters["theta"] = 0.5
apec.time_parameters["write_substep_results"] = 1
apec.time_parameters["transient_effects"] = 0
apec.time_parameters["time_step_number"] = 1

Selector.select(apec.mesh, "CS3L", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_CS3L / a_CS, 2, type="elem")
Selector.select(apec.mesh, "CS2L", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_CS2L / a_CS, 2, type="elem")
Selector.select(apec.mesh, "CS1L", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_CS1L / a_CS, 2, type="elem")
Selector.select(apec.mesh, "CS1U", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_CS1U / a_CS, 2, type="elem")
Selector.select(apec.mesh, "CS2U", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_CS2U / a_CS, 2, type="elem")
Selector.select(apec.mesh, "CS3U", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_CS3U / a_CS, 2, type="elem")

Selector.select(apec.mesh, "PFC_1", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_PF1 / a_PF1, 2, type="elem")
Selector.select(apec.mesh, "PFC_2", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_PF2 / a_PF2, 2, type="elem")
Selector.select(apec.mesh, "PFC_3", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_PF3 / a_PF3, 2, type="elem")
Selector.select(apec.mesh, "PFC_4", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_PF4 / a_PF4, 2, type="elem")
Selector.select(apec.mesh, "PFC_5", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_PF5 / a_PF5, 2, type="elem")
Selector.select(apec.mesh, "PFC_6", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", i_PF6 / a_PF6, 2, type="elem")

Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem")
Selector.set_value_on_selection(apec.mesh, "Js", 15e6 / a_plasma, 2, type="elem")

apec.write_load_step(1)

apec.save_mesh()

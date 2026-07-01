from apec_m import APEC_M
from apec_m.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

aplasma = 0

filename = "ITER_2D_APEC_carma"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.algorithm.input_parameters["target_central_axis_x"] = 6.52136102486433
apec.algorithm.input_parameters["target_central_axis_y"] = 0.536537145137914
apec.algorithm.input_parameters["target_central_axis_z"] = 0.0
apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["control_coils"] = [
	{"component": "PFC_2", "radial_gain": 0.8, "vertical_gain": 0.0},
	{"component": "PFC_3", "radial_gain": 0.0, "vertical_gain": 0.8},
	{"component": "PFC_4", "radial_gain": 0.0, "vertical_gain": 0.8},
	{"component": "PFC_5", "radial_gain": -0.8, "vertical_gain": 0.0},
]
apec.algorithm.input_parameters["tolerance_value"] = 1e-6
apec.algorithm.input_parameters["bradx_value"] = -3.873415861e-08
apec.algorithm.input_parameters["brady_value"] = 5.313857037e-09
apec.algorithm.input_parameters["bverx_value"] = 1.822467692e-10
apec.algorithm.input_parameters["bvery_value"] = 1.133287835e-07
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["plasma_current"] = 15e6
apec.algorithm.input_parameters["beta_p"] = 6.8216e-001
apec.algorithm.input_parameters["print_full_output"] = 0

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

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

apec.time_step["theta"] = 0.5
apec.time_step["transient_effects"] = 0
apec.time_step["time_step_number"] = 1

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

apec.load_mesh_into_solver()

apec.algorithm.equilibrium_solve()

apec.write_output()

apec.save_mesh()

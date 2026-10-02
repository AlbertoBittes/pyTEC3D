from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

pytec3d = pyTEC3DSolver()

aplasma = 0

filename = "ITER_2D_coils"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.algorithm.input_parameters["target_central_axis_x"] = 6.52136102486433
pytec3d.algorithm.input_parameters["target_central_axis_y"] = 0.536537145137914
pytec3d.algorithm.input_parameters["target_central_axis_z"] = 0.0
pytec3d.algorithm.input_parameters["R0"] = 6.2
pytec3d.algorithm.input_parameters["max_picard_iterations"] = 70
pytec3d.algorithm.input_parameters["control_coils"] = [
	{"component": "PFC_2", "radial_gain": 0.8, "vertical_gain": 0.0},
	{"component": "PFC_3", "radial_gain": 0.0, "vertical_gain": 0.8},
	{"component": "PFC_4", "radial_gain": 0.0, "vertical_gain": 0.8},
	{"component": "PFC_5", "radial_gain": -0.8, "vertical_gain": 0.0},
]
pytec3d.algorithm.input_parameters["tolerance_value"] = 1e-6
pytec3d.algorithm.input_parameters["bradx_value"] = -3.873415861e-08
pytec3d.algorithm.input_parameters["brady_value"] = 5.313857037e-09
pytec3d.algorithm.input_parameters["bverx_value"] = 1.822467692e-10
pytec3d.algorithm.input_parameters["bvery_value"] = 1.133287835e-07
pytec3d.algorithm.input_parameters["alpha_mp"] = 1.0001
pytec3d.algorithm.input_parameters["alpha_np"] = 0.83197
pytec3d.algorithm.input_parameters["plasma_current"] = 15e6
pytec3d.algorithm.input_parameters["beta_p"] = 6.8216e-001

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

pytec3d.model.read_excel("input_template.xlsx")

Selector.select_all(pytec3d.model.mesh)
Selector.nsel_ext(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.model.time_step["time"] = 1.0
pytec3d.model.time_step["theta"] = 0.5
pytec3d.model.time_step["transient_effects"] = 0

Selector.select(pytec3d.model.mesh, "CS3L", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS3L / a_CS, 2, type="elem")
Selector.select(pytec3d.model.mesh, "CS2L", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS2L / a_CS, 2, type="elem")
Selector.select(pytec3d.model.mesh, "CS1L", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS1L / a_CS, 2, type="elem")
Selector.select(pytec3d.model.mesh, "CS1U", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS1U / a_CS, 2, type="elem")
Selector.select(pytec3d.model.mesh, "CS2U", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS2U / a_CS, 2, type="elem")
Selector.select(pytec3d.model.mesh, "CS3U", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS3U / a_CS, 2, type="elem")

Selector.select(pytec3d.model.mesh, "PFC_1", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF1 / a_PF1, 2, type="elem")
Selector.select(pytec3d.model.mesh, "PFC_2", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF2 / a_PF2, 2, type="elem")
Selector.select(pytec3d.model.mesh, "PFC_3", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF3 / a_PF3, 2, type="elem")
Selector.select(pytec3d.model.mesh, "PFC_4", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF4 / a_PF4, 2, type="elem")
Selector.select(pytec3d.model.mesh, "PFC_5", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF5 / a_PF5, 2, type="elem")
Selector.select(pytec3d.model.mesh, "PFC_6", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF6 / a_PF6, 2, type="elem")

Selector.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "Js", 15e6 / a_plasma, 2, type="elem")

pytec3d.load_model_into_solver()

pytec3d.algorithm.equilibrium_solve()

pytec3d.write_output()

pytec3d.save_model()

from apec_m import APEC_M
from vtk_grid_helper import VTKGridHelper
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

aplasma=0

filename = "ITER_2D_APEC"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

VTKGridHelper.reset_vtk_arrays_to_default((apec.mesh))

apec.apec_settings = {
    "target_central_axis_x": 6.52136102486433,
    "target_central_axis_y": 0.536537145137914,
    "target_central_axis_z": 0.0,
    "target_central_axis_element_id": 651,
    "R0": 6.2,
    "max_picard_iterations": 70,
    "min_distance": 0.6,
    "cp_gain": 0.250,
    "tolerance_value": 1e-3,
    "bradx_value":  -3.873415861E-08,
    "brady_value":  5.313857037E-09,
    "bverx_value": 1.822467692E-10,
    "bvery_value": 1.133287835E-07,
    "alpha_value": 8.3197e-01,
    "plasma_current": 15e6,
    "beta_p": 6.8216e-001,
    "print_full_output": 1,
}

i_CS3U=-4.8204e+06
i_CS2U=-4.0480e+06
i_CS1U=-2.0804e+07
i_CS1L=-2.0804e+07
i_CS2L=-1.0297e+07
i_CS3L=5.3012e+06

i_PF1=5.7325e+06
i_PF2=-2.8834e+06
i_PF3=-5.8084e+06
i_PF4=-4.7762e+06
i_PF5=-7.8133e+06
i_PF6=1.6941e+07

a_PF1=0.94375190
a_PF2=0.41453946
a_PF3=0.66413094
a_PF4=0.60871516
a_PF5=0.77496250
a_PF6=1.72659250
a_CS=1.55608
a_plasma=28.6027095

apec.read_excel("input_template.xlsx")

apec.write_model_data()

apec.select_all()
apec.nsel_ext()
apec.set_value_on_selection("A", 0, 2)

apec.time_parameters["ramped_loads"] = 0
apec.time_parameters["theta"] = 0.5
apec.time_parameters["write_substep_results"] = 1
apec.time_parameters["transient_effects"] = 0
apec.time_parameters["time_step_number"] = 1

apec.select("CS3L", 1, type="elem")
apec.set_value_on_selection("Js", i_CS3L / a_CS, 2, type="elem")
apec.select("CS2L", 1, type="elem")
apec.set_value_on_selection("Js", i_CS2L / a_CS, 2, type="elem")
apec.select("CS1L", 1, type="elem")
apec.set_value_on_selection("Js", i_CS1L / a_CS, 2, type="elem")
apec.select("CS1U", 1, type="elem")
apec.set_value_on_selection("Js", i_CS1U / a_CS, 2, type="elem")
apec.select("CS2U", 1, type="elem")
apec.set_value_on_selection("Js", i_CS2U / a_CS, 2, type="elem")
apec.select("CS3U", 1, type="elem")
apec.set_value_on_selection("Js", i_CS3U / a_CS, 2, type="elem")

apec.select("PFC_1", 1, type="elem")
apec.set_value_on_selection("Js", i_PF1 / a_PF1, 2, type="elem")
apec.select("PFC_2", 1, type="elem")
apec.set_value_on_selection("Js", i_PF2 / a_PF2, 2, type="elem")
apec.select("PFC_3", 1, type="elem")
apec.set_value_on_selection("Js", i_PF3 / a_PF3, 2, type="elem")
apec.select("PFC_4", 1, type="elem")
apec.set_value_on_selection("Js", i_PF4 / a_PF4, 2, type="elem")
apec.select("PFC_5", 1, type="elem")
apec.set_value_on_selection("Js", i_PF5 / a_PF5, 2, type="elem")
apec.select("PFC_6", 1, type="elem")
apec.set_value_on_selection("Js", i_PF6 / a_PF6, 2, type="elem")

apec.select("VACUUM_CHAMBER", 1, type="elem")
apec.set_value_on_selection("Js", 15e6 / a_plasma, 2, type="elem")

apec.write_load_step(1)

apec.save_mesh()

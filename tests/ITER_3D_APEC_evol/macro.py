from apec_m import APEC_M
import vtk.util.numpy_support as vtknp
import numpy as np
from selector import Selector
from coordinate_system import CylindricalCoordinateSystem, SphericalCoordinateSystem
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_equilibrium_output = APEC_M()

eq_filename = "ITER_2D_APEC"
filename = "ITER_3D_APEC_evol"

eq_filepath = os.path.join(os.getcwd(), "tests", eq_filename, eq_filename + ".vtu")
filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")
static_mesh_filepath = os.path.join(os.getcwd(), "tests", filename, filename + "_1.vtu")

# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    eq_filename,
    eq_filename,
    "output",
    eq_filename + "_out_0.pvtu",
)
apec_equilibrium_output.read_mesh(output_eq_mesh_filepath, False)

eps = 1e-1
theta_eps = 1e-3
# Read the static mesh
apec.read_mesh(static_mesh_filepath, False)

# Set Boundary conditions
apec.select_pos(
    0, (0.0 - eps, 0.0 + eps), coordinate_system=CylindricalCoordinateSystem()
)
# Az bound to 0 on vertical axis
apec.set_value_on_selection("A", 0, 2)

apec.select_pos(
    0, (100.0 - eps, 100.0 + eps), coordinate_system=SphericalCoordinateSystem()
)
# Far field zero boundary conditions
apec.set_value_on_selection("A", 0, 2)

# Initial loads: Project the current density from the equilibrium output to the static mesh
apec.select_all()
apec_equilibrium_output.select_all()
apec_equilibrium_output.select("subdomain_id", 23, type="elem")
apec.toroidal_projection(
    apec_equilibrium_output.mesh, "Jt_z", 0, "Js", 2, location="elem"
)

apec.write_loads_bcs()
apec.save_mesh()

# Now start the evolution analysis
apec.read_mesh(filepath, False)
apec.remove_array("boundary_id_0")
apec.remove_array("boundary_id_1")
apec.remove_array("boundary_id_2")
apec.add_boundary_array()
apec.add_boundary_array()

# Define periodic boundaries
apec.select_all()
apec.nsel_ext()
apec.select("ansys_material_type", 24, type="elem")
apec.nsle(mode="r")
apec.select_pos(
    1,
    (0.0 - theta_eps, 0.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
apec.set_value_on_selection("boundary_id_1", 4)

apec.select_all()
apec.nsel_ext()
apec.select("ansys_material_type", 24, type="elem")
apec.nsle(mode="r")
apec.select_pos(
    1,
    (2 * np.pi / 180.0 - theta_eps, 2 * np.pi / 180.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
apec.set_value_on_selection("boundary_id_1", 5)

apec.select_all()
apec.nsel_ext()
apec.select_pos(
    1,
    (0.0 - theta_eps, 0.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
apec.select_pos(
    0,
    (0.0 - eps, 0.0 + eps),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
apec.select_pos(
    0,
    (100.0 - eps, 100.0 + eps),
    mode="u",
    coordinate_system=SphericalCoordinateSystem(),
)
apec.set_value_on_selection("boundary_id_0", 2)

apec.select_all()
apec.nsel_ext()
apec.select_pos(
    1,
    (2 * np.pi / 180.0 - theta_eps, 2 * np.pi / 180.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
apec.select_pos(
    0,
    (0.0 - eps, 0.0 + eps),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
apec.select_pos(
    0,
    (100.0 - eps, 100.0 + eps),
    mode="u",
    coordinate_system=SphericalCoordinateSystem(),
)
apec.set_value_on_selection("boundary_id_0", 3)

# Set Boundary conditions
apec.select_pos(
    0,
    (0.0 - eps, 0.0 + eps),
    coordinate_system=CylindricalCoordinateSystem(),
)
# Ax and Az bound to 0 on vertical axis
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 2)

apec.select_all()
apec.select_pos(
    0,
    (100.0 - eps, 100.0 + eps),
    coordinate_system=SphericalCoordinateSystem(),
)
# Far field zero boundary conditions
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 1)
apec.set_value_on_selection("A", 0, 2)

apec.write_loads_bcs()
apec.save_mesh()
# Minus 2 degrees rotation
matrix_2deg = [
    [np.cos(np.deg2rad(-2)), 0.0, np.sin(np.deg2rad(-2)), 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [-np.sin(np.deg2rad(-2)), 0.0, np.cos(np.deg2rad(-2)), 0.0],
    [0.0, 0.0, 0.0, 1.0],
]

apec.periodic_boundaries.append(
    {
        "main_boundary": 3,
        "paired_boundary": 2,
        "periodic_matrix": matrix_2deg,
        "variable": "A",
    }
)

apec.periodic_boundaries.append(
    {
        "main_boundary": 5,
        "paired_boundary": 4,
        "periodic_matrix": matrix_2deg,
        "variable": "VOLT",
    }
)

apec.read_excel("input_template.xlsx")
apec.write_model_data()

apec.apec_settings = {
    "R0": 6.2,
    "max_picard_iterations": 70,
    "tolerance_value": 1e-3,
    "alpha_value": 8.3197e-01,
    "plasma_current": 15e6,
    "beta_p": 6.8216e-001,
    "evolution_mesh": "ITER_3D_APEC_evol",
    "equilibrium_mesh": "ITER_3D_APEC_evol_1",
    "print_full_output": 1,
}


# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt("time_steps.txt")
currents = np.loadtxt("plasma_current.txt")
betas = np.loadtxt("beta_p.txt")

# Consider only the time steps every 5 in the list
time_steps = time_steps[::5][1:]
currents = currents[::5][1:]
betas = betas[::5][1:]

for step_idx in range(len(time_steps)):
    apec.time_parameters["time"] = time_steps[step_idx]
    apec.apec_settings["plasma_current"] = currents[step_idx]
    apec.apec_settings["beta_p"] = betas[step_idx]
    apec.write_time_step(step_idx + 1)

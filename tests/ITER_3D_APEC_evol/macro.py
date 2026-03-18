from apec_m import APEC_M
import vtk.util.numpy_support as vtknp
import numpy as np
from selector import Selector
from model_edit import ModelEdit
from cyclic import Cyclic
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
max_rad = 100.0
# Read the static mesh
apec_equilibrium = APEC_M()
apec_equilibrium.read_mesh(static_mesh_filepath)

# Set Boundary conditions
apec_equilibrium.select_pos(
    0, (0.0 - eps, 0.0 + eps), coordinate_system=CylindricalCoordinateSystem()
)
# Az bound to 0 on vertical axis
apec_equilibrium.set_value_on_selection("A", 0, 2)
apec_equilibrium.select_pos(
    0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem()
)
# Far field zero boundary conditions
apec_equilibrium.set_value_on_selection("A", 0, 2)

apec_equilibrium.write_loads_bcs()
apec_equilibrium.save_mesh()

# Now start the evolution analysis
apec.read_mesh(filepath)

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.select(apec.mesh, "ansys_material_type", 24, type="elem")
apec.nsle(mode="r")
Selector.select_pos(apec.mesh, 
    1,
    (0.0 - theta_eps, 0.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "boundary_id_1", 4)

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.select(apec.mesh, "ansys_material_type", 24, type="elem")
apec.nsle(mode="r")
Selector.select_pos(apec.mesh, 
    1,
    (2 * np.pi / 180.0 - theta_eps, 2 * np.pi / 180.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "boundary_id_1", 5)

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.select_pos(apec.mesh, 
    1,
    (0.0 - theta_eps, 0.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(apec.mesh, 
    0,
    (0.0 - eps, 0.0 + eps),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(apec.mesh, 
    0,
    (max_rad - eps, max_rad + eps),
    mode="u",
    coordinate_system=SphericalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "boundary_id_0", 2)

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.select_pos(apec.mesh, 
    1,
    (2 * np.pi / 180.0 - theta_eps, 2 * np.pi / 180.0 + theta_eps),
    mode="r",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(apec.mesh, 
    0,
    (0.0 - eps, 0.0 + eps),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(apec.mesh, 
    0,
    (max_rad - eps, max_rad + eps),
    mode="u",
    coordinate_system=SphericalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "boundary_id_0", 3)

# Set Boundary conditions
Selector.select_pos(apec.mesh, 
    0,
    (0.0 - eps, 0.0 + eps),
    coordinate_system=CylindricalCoordinateSystem(),
)
# Ax and Az bound to 0 on vertical axis
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select_all(apec.mesh)
Selector.select_pos(apec.mesh, 
    0,
    (max_rad - eps, max_rad + eps),
    coordinate_system=SphericalCoordinateSystem(),
)
# Far field zero boundary conditions
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "ID", 16170)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

Selector.select_all(apec.mesh)
apec_equilibrium.select_all()
apec.apec_settings = {
    "R0": 6.2,
    "max_picard_iterations": 70,
    "tolerance_value": 1e-3,
    "alpha_value": 8.3197e-01,
    "plasma_current": 15e6,
    "beta_p": 6.8216e-001,
    "evolution_mesh": "ITER_3D_APEC_evol",
    "static_mesh": "ITER_3D_APEC_evol_static",
    "equilibrium_mesh": "ITER_3D_APEC_evol_equilibrium",
    "print_full_output": 1,
    "toroidal_slices_nodes": Cyclic.build_toroidal_slices(
        apec_equilibrium.mesh, apec.mesh, "node", n_slices=3, tolerance=1e-3
    ),
    "toroidal_slices_elems": Cyclic.build_toroidal_slices(
        apec_equilibrium.mesh, apec.mesh, "elem", n_slices=2, tolerance=1e-3
    ),
}

# Save the toroidal projection maps to files
np.savetxt("toroidal_slices_nodes.txt", apec.apec_settings["toroidal_slices_nodes"].astype(int), fmt="%d")
np.savetxt("toroidal_slices_elems.txt", apec.apec_settings["toroidal_slices_elems"].astype(int), fmt="%d")

# Initial loads: Project the current density from the equilibrium output to the static mesh
Selector.select(apec.mesh, "ansys_material_type", 23, type="elem")
apec_equilibrium_output.select_all()
apec_equilibrium_output.select("subdomain_id", 23, type="elem")
Cyclic.toroidal_projection(
    apec_equilibrium_output.mesh,
    "Jt_z",
    0,
    apec.mesh,
    "Js",
    1,
    location="elem",
    n_slices=2,
    tolerance=1e-3,
)
ModelEdit.map_vector_field_to_global_components(apec.mesh, "Js", CylindricalCoordinateSystem(), location="elem")

apec.write_loads_bcs()
apec.save_mesh()
apec.read_excel("input_template.xlsx")
apec.write_model_data()


# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt("time_steps.txt")
currents = np.loadtxt("plasma_current.txt")
betas = np.loadtxt("beta_p.txt")

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::10]))
currents = np.concatenate((currents[0:1], currents[1::10]))
betas = np.concatenate((betas[0:1], betas[1::10]))

for step_idx in range(len(time_steps)):
    apec.time_parameters["time"] = time_steps[step_idx]
    apec.apec_settings["plasma_current"] = currents[step_idx]
    apec.apec_settings["beta_p"] = betas[step_idx]
    apec.write_time_step(step_idx + 1)
    if step_idx == 0:
        apec.apec_settings.pop("toroidal_slices_nodes", None)
        apec.apec_settings.pop("toroidal_slices_elems", None)

from apec_m import APEC_M
from mesh_edit import MeshEdit
import vtk.util.numpy_support as vtknp
import numpy as np
from selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_equilibrium_output = APEC_M()

eq_filename = "ITER_2D_APEC"
filename = "ITER_2D_APEC_evol"
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

# Read the static mesh
apec_equilibrium = APEC_M()
apec_equilibrium.read_mesh(static_mesh_filepath, False)
apec_equilibrium.read_mesh(static_mesh_filepath)
apec_equilibrium.read_excel("input_template.xlsx")
apec_equilibrium.write_model_data()

# # Set conductors to non conductors
# apec.select("ansys_material_type", 24, type="elem")
# apec.set_value_on_selection("ansys_material_type", 23, type="elem")
# apec.select("ansys_etype", 24, type="elem")
# apec.set_value_on_selection("ansys_etype", 23, type="elem")

# Set Boundary conditions
apec_equilibrium.select_all()
apec_equilibrium.nsel_ext()
apec_equilibrium.set_value_on_selection("A", 0, 2)

# Project the current density from the equilibrium output to the static mesh
apec_equilibrium.select_all()
apec_equilibrium_output.select_all()
apec_equilibrium_output.select("subdomain_id", 23, type="elem")
apec_equilibrium.toroidal_projection(
    apec_equilibrium_output.mesh, "Jt_z", 0, "Js", 2, location="elem", n_slices=1, tolerance=1e-3
)

# Now start the evolution analysis
apec.read_mesh(filepath)
apec.select_all()

apec.apec_settings = {
    "R0": 6.2,
    "max_picard_iterations": 70,
    "tolerance_value": 1e-3,
    "alpha_value": 8.3197e-01,
    "plasma_current": 15e6,
    "beta_p": 6.8216e-001,
    "evolution_mesh": "ITER_2D_APEC_evol",
    "equilibrium_mesh": "ITER_2D_APEC_evol_1",
    "print_full_output": 1,
    "toroidal_projection_map_nodes": apec.build_toroidal_projection_map(apec_equilibrium.mesh, "node", n_slices=1, tolerance=1e-3),
    "toroidal_projection_map_elems": apec.build_toroidal_projection_map(apec_equilibrium.mesh, "elem", n_slices=1, tolerance=1e-3),

}

# Write loads and BCs and save the static mesh
apec_equilibrium.write_loads_bcs()
apec_equilibrium.save_mesh()

# Save the toroidal projection maps to files
np.savetxt("toroidal_projection_map_nodes.txt", apec.apec_settings["toroidal_projection_map_nodes"].astype(int), fmt="%d")
np.savetxt("toroidal_projection_map_elems.txt", apec.apec_settings["toroidal_projection_map_elems"].astype(int), fmt="%d")

apec.select_all()
apec.nsel_ext()
apec.set_value_on_selection("A", 0, 2)

apec.select_all()
apec_equilibrium_output.select_all()
apec_equilibrium_output.select("subdomain_id", 23, type="elem")
apec.toroidal_projection(
    apec_equilibrium_output.mesh, "Jt_z", 0, "Js", 2, location="elem"
)

apec.write_loads_bcs()

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
    if step_idx == 0:
        apec.apec_settings.pop("toroidal_projection_map_nodes", None)
        apec.apec_settings.pop("toroidal_projection_map_elems", None)
apec.save_mesh()

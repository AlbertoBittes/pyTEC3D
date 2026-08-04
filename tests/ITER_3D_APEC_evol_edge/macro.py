from apec_m import APEC_M
import vtk.util.numpy_support as vtknp
import numpy as np
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.tree_gauge import TreeGauge
from apec_m.coordinate_system import CylindricalCoordinateSystem, SphericalCoordinateSystem
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_evolution_mesh = "ITER_3D_APEC_evol_edge"
apec_evolution_mesh_notrans = "ITER_3D_APEC_evol_edge_1"
apec_equilibrium_output = "ITER_2D_APEC_out_0"


evol_filepath = os.path.join(os.getcwd(),"tests", apec_evolution_mesh, apec_evolution_mesh + ".vtu")
static_mesh_filepath = os.path.join(os.getcwd(),"tests", apec_evolution_mesh, apec_evolution_mesh_notrans + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    "ITER_2D_APEC",
    "ITER_2D_APEC",
    "output",
    apec_equilibrium_output + ".pvtu",
)


eps = 1e-1
theta_eps = 1e-3
max_rad = 100
# Read the static mesh
apec.read_mesh(static_mesh_filepath)

# Set Boundary conditions
Selector().select_pos(
    apec.mesh, 0, (0.0 - eps, 0.0 + eps), coordinate_system=CylindricalCoordinateSystem()
)
# Az bound to 0 on vertical axis
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)
Selector.select_pos(
    apec.mesh, 0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem()
)
# Far field zero boundary conditions
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

apec.write_loads_bcs()
apec.save_mesh()

# Now start the evolution analysis
apec.read_mesh(evol_filepath)

# Set Boundary conditions
Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "Az", 0, 0)

TreeGauge.tree_gauge(apec.mesh)

Selector.select(apec.mesh, "ID", 16170)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

Selector.select_all(apec.mesh)
apec.apec_settings = {
    "R0": 6.2,
    "max_picard_iterations": 70,
    "tolerance_value": 1e-3,
    "alpha_value": 8.3197e-01,
    "plasma_current": 15e6,
    "beta_p": 6.8216e-001,
    "evolution_mesh": "ITER_3D_APEC_evol_edge",
    "equilibrium_mesh": "ITER_3D_APEC_evol_edge_1",
    "print_full_output": 1,
    "toroidal_slices_nodes": apec.algorithm.build_toroidal_slices(
        apec.get_mesh(apec_evolution_mesh_notrans), apec.mesh, "node", n_slices=3, tolerance=1e-3
    ),
    "toroidal_slices_elems": apec.algorithm.build_toroidal_slices(
        apec.get_mesh(apec_evolution_mesh_notrans), apec.mesh, "elem", n_slices=2, tolerance=1e-3
    ),
}

# Save the toroidal projection maps to files
np.savetxt("toroidal_slices_nodes.txt", apec.apec_settings["toroidal_slices_nodes"].astype(int), fmt="%d")
np.savetxt("toroidal_slices_elems.txt", apec.apec_settings["toroidal_slices_elems"].astype(int), fmt="%d")

# Initial loads: Project the current density from the equilibrium output to the static mesh
apec.read_mesh(output_eq_mesh_filepath)
apec.set_active_mesh(apec_evolution_mesh)
Selector.select_all(apec.mesh)
Selector.select_all(apec.get_mesh(apec_equilibrium_output))
Selector.select(apec.get_mesh(apec_equilibrium_output), "subdomain_id", 23, type="elem")
Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
apec.algorithm.equilibrium_mesh = apec.get_mesh(apec_equilibrium_output)
apec.algorithm._toroidal_slices_nodes = apec.apec_settings["toroidal_slices_nodes"].astype(np.int64)
apec.algorithm._toroidal_slices_elems = apec.apec_settings["toroidal_slices_elems"].astype(np.int64)
apec.algorithm._toroidal_projection_current()
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
time_steps = time_steps[::5][1:]
currents = currents[::5][1:]
betas = betas[::5][1:]

for step_idx in range(len(time_steps)):
    apec.time_parameters["time"] = time_steps[step_idx]
    apec.apec_settings["plasma_current"] = currents[step_idx]
    apec.apec_settings["beta_p"] = betas[step_idx]
    apec.write_time_step(step_idx + 1)
    if step_idx == 0:
        apec.apec_settings.pop("toroidal_slices_nodes", None)
        apec.apec_settings.pop("toroidal_slices_elems", None)

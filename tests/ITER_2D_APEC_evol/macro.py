from apec_m import APEC_M
from model_edit import ModelEdit
from cyclic import Cyclic
import vtk.util.numpy_support as vtknp
import numpy as np
from selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_evolution_mesh = "ITER_2D_APEC_evol"
apec_evolution_mesh_notrans = "ITER_2D_APEC_evol_equilibrium"
apec_equilibrium_output = "ITER_2D_APEC_out_0"


evol_filepath = os.path.join(os.getcwd(), "tests", apec_evolution_mesh, apec_evolution_mesh + ".vtu")
static_mesh_filepath = os.path.join(os.getcwd(), "tests", apec_evolution_mesh, apec_evolution_mesh_notrans + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    "ITER_2D_APEC",
    "ITER_2D_APEC",
    "output",
    apec_equilibrium_output + ".pvtu",
)


apec.read_mesh(evol_filepath)
apec.read_mesh(static_mesh_filepath)
apec.read_mesh(output_eq_mesh_filepath, False)

apec.read_excel("input_template.xlsx")
apec.write_model_data()


# Set Boundary conditions
apec.set_active_mesh(apec_evolution_mesh_notrans)

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)
Selector.select_all(apec.mesh)

# Write loads and BCs and save the static mesh
apec.write_loads_bcs()

apec.save_mesh()

# Now start the evolution analysis
apec.set_active_mesh(apec_evolution_mesh)
Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

# 1. Build source->destination element map in VACUUM_CHAMBER using toroidal slices
Selector.select(apec.get_mesh(apec_equilibrium_output), "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.select(apec.get_mesh(apec_equilibrium_output), "COILS", 1, type="elem", mode="a")
# Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")

elem_id_map = Cyclic.build_toroidal_slices(
    apec.get_mesh(apec_equilibrium_output),
    apec.mesh,
    location="elem",
    n_slices=1,
    tolerance=1e-3,
)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
jt_z_array = vtknp.vtk_to_numpy(apec.get_mesh(apec_equilibrium_output).GetCellData().GetArray("Jt_z"))
js_array = vtknp.vtk_to_numpy(apec.mesh.GetCellData().GetArray("Js"))

src_elem_ids = elem_id_map[:, 0]
dst_elem_ids = elem_id_map[:, 1]
js_array[dst_elem_ids, 2] = jt_z_array[src_elem_ids]
apec.mesh.GetCellData().GetArray("Js").Modified()

apec.write_loads_bcs()

Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.nsle(apec.mesh)
apec.apec_settings["target_central_axis_x"] = 6.52136102486433
apec.apec_settings["target_central_axis_y"] = 0.536537145137914
apec.apec_settings["target_central_axis_z"] = 0.0
apec.apec_settings["R0"] = 6.2
apec.apec_settings["max_picard_iterations"] = 70
apec.apec_settings["min_distance"] = 0.1
apec.apec_settings["cp_gain"] = 0.250
apec.apec_settings["tolerance_value"] = 1e-3
apec.apec_settings["bradx_value"] = -3.873415861e-08
apec.apec_settings["brady_value"] = 5.313857037e-09
apec.apec_settings["bverx_value"] = 1.822467692e-10
apec.apec_settings["bvery_value"] = 1.133287835e-07
apec.apec_settings["alpha_mp"] = 1.0001
apec.apec_settings["alpha_np"] = 0.83197
apec.apec_settings["equilibrium_mesh"] = apec_evolution_mesh_notrans
apec.apec_settings["static_mesh"] = apec_evolution_mesh_notrans
apec.apec_settings["evolution_mesh"] = apec_evolution_mesh
apec.apec_settings["AF0"] = 5.3 * 6.2
# apec.apec_settings["plasma_current"] = 15e6
# apec.apec_settings["beta_p"] = 6.8216e-001
apec.apec_settings["print_full_output"] = 0
apec.apec_settings["toroidal_slices_nodes"] = Cyclic.build_toroidal_slices(
    apec.mesh, apec.mesh, "node", n_slices=1, tolerance=1e-3
)
apec.apec_settings["toroidal_slices_elems"] = Cyclic.build_toroidal_slices(
    apec.mesh, apec.mesh, "elem", n_slices=1, tolerance=1e-3
)

# Save the toroidal projection maps to files
np.savetxt("toroidal_slices_nodes.txt", apec.apec_settings["toroidal_slices_nodes"].astype(int), fmt="%d")
np.savetxt("toroidal_slices_elems.txt", apec.apec_settings["toroidal_slices_elems"].astype(int), fmt="%d")

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt("time_steps.txt")
currents = np.loadtxt("plasma_current.txt")
betas = np.loadtxt("beta_p.txt")

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::10]))
currents = np.concatenate((currents[0:1], currents[1::10]))
betas = np.concatenate((betas[0:1], betas[1::10]))
apec.time_parameters["theta"] = 1.0
for step_idx in range(len(time_steps)):
    apec.time_parameters["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.time_parameters["transient_effects"] = 0
        apec.time_parameters["time"] = 1e-10
    else:
        apec.time_parameters["transient_effects"] = 1
    apec.apec_settings["plasma_current"] = currents[step_idx]
    apec.apec_settings["beta_p"] = betas[step_idx]
    apec.write_time_step(step_idx + 1)
    if step_idx == 0:
        apec.apec_settings.pop("toroidal_slices_nodes", None)
        apec.apec_settings.pop("toroidal_slices_elems", None)

apec.save_mesh()

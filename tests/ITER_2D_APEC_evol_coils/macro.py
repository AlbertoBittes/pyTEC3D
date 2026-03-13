from apec_m import APEC_M
from model_edit import ModelEdit
from cyclic import Cyclic
import vtk.util.numpy_support as vtknp
import numpy as np
from selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_evolution_mesh = "ITER_2D_APEC_evol_coils"
apec_evolution_mesh_notrans = "ITER_2D_APEC_evol_coils_equilibrium"
apec_equilibrium_output = "ITER_2D_APEC_coils_out_0"


evol_filepath = os.path.join(os.getcwd(), "tests", apec_evolution_mesh, apec_evolution_mesh + ".vtu")
static_mesh_filepath = os.path.join(os.getcwd(), "tests", apec_evolution_mesh, apec_evolution_mesh_notrans + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    "ITER_2D_APEC_coils",
    "ITER_2D_APEC_coils",
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

# 1. Select COILS and VACUUM_CHAMBER in output mesh
Selector.select(apec.get_mesh(apec_equilibrium_output), "COILS", 1, type="elem", mode="s")
Selector.select(apec.get_mesh(apec_equilibrium_output), "VACUUM_CHAMBER", 1, type="elem", mode="a")

# 2. Extract selection mask, Jt_z, and libmesh_elem_id arrays from output mesh
output_mesh = apec.get_mesh(apec_equilibrium_output)
selection_mask = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("selection_cell"))
jt_z_array = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("Jt_z"))
libmesh_elem_id_out = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("libmesh_elem_id"))

# 3. Get selected libmesh ids and corresponding Jt_z values
selected = selection_mask == 1
selected_libmesh_ids = libmesh_elem_id_out[selected]
selected_jt_z = jt_z_array[selected]

# The libmesh_elem_id from output mesh corresponds directly to VTK cell ids in apec.mesh
js_array = vtknp.vtk_to_numpy(apec.mesh.GetCellData().GetArray("Js"))
# Only update valid VTK cell ids (must be within js_array shape)
valid_mask = (selected_libmesh_ids >= 0) & (selected_libmesh_ids < js_array.shape[0])
js_array[selected_libmesh_ids[valid_mask], 2] = selected_jt_z[valid_mask]
apec.mesh.GetCellData().GetArray("Js").Modified()

# Write loads and BCs and save the static mesh
apec.write_loads_bcs()

apec.save_mesh()

# Now start the evolution analysis
apec.set_active_mesh(apec_evolution_mesh)
Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

for k, comp in enumerate(["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L"]):
    Selector.select(apec.mesh, comp, 1, type="elem", mode="s")
    Selector.nsle(apec.mesh)
    Selector.set_value_on_selection(apec.mesh, "CP_0", 11 + k)
    master_node_cs = vtknp.vtk_to_numpy(apec.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]
    Selector.select(apec.mesh, "ID", master_node_cs)
    mask = apec.real_constants["ansys_real_constant"] == 11 + k
    n_turns = apec.real_constants[mask]["NC"].values[0]
    resistance = apec.real_constants[mask]["R"].values[0]
    Selector.select(apec.get_mesh(apec_equilibrium_output), comp, 1, type="elem", mode="s")
    output_mesh = apec.get_mesh(apec_equilibrium_output)
    selection_mask = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("selection_cell"))
    jt_z_array = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("Jt_z"))
    jt_z = jt_z_array[selection_mask == 1][0]
    sc = apec.real_constants[mask]["SC"].values[0]
    current = sc * jt_z
    voltage = - resistance * current / n_turns
    Selector.set_value_on_selection(apec.mesh, "V", voltage)

for k, comp in enumerate(["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6"]):
    Selector.select(apec.mesh, comp, 1, type="elem", mode="s")
    Selector.nsle(apec.mesh)
    Selector.set_value_on_selection(apec.mesh, "CP_0", 1 + k)
    master_node_pfcs = vtknp.vtk_to_numpy(apec.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]
    Selector.select(apec.mesh, "ID", master_node_pfcs)
    mask = apec.real_constants["ansys_real_constant"] == 1 + k
    n_turns = apec.real_constants[mask]["NC"].values[0]
    resistance = apec.real_constants[mask]["R"].values[0]
    Selector.select(apec.get_mesh(apec_equilibrium_output), comp, 1, type="elem", mode="s")
    output_mesh = apec.get_mesh(apec_equilibrium_output)
    selection_mask = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("selection_cell"))
    jt_z_array = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("Jt_z"))
    jt_z = jt_z_array[selection_mask == 1][0]
    sc = apec.real_constants[mask]["SC"].values[0]
    current = sc * jt_z
    voltage = - resistance * current / n_turns
    Selector.set_value_on_selection(apec.mesh, "V", voltage)

# 1. Select COILS and VACUUM_CHAMBER in output mesh
Selector.select(apec.get_mesh(apec_equilibrium_output), "VACUUM_CHAMBER", 1, type="elem", mode="s")

# 2. Extract selection mask, Jt_z, and libmesh_elem_id arrays from output mesh
output_mesh = apec.get_mesh(apec_equilibrium_output)
selection_mask = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("selection_cell"))
jt_z_array = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("Jt_z"))
libmesh_elem_id_out = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("libmesh_elem_id"))

# 3. Get selected libmesh ids and corresponding Jt_z values
selected = selection_mask == 1
selected_libmesh_ids = libmesh_elem_id_out[selected]
selected_jt_z = jt_z_array[selected]

# The libmesh_elem_id from output mesh corresponds directly to VTK cell ids in apec.mesh
js_array = vtknp.vtk_to_numpy(apec.mesh.GetCellData().GetArray("Js"))
# Only update valid VTK cell ids (must be within js_array shape)
valid_mask = (selected_libmesh_ids >= 0) & (selected_libmesh_ids < js_array.shape[0])
js_array[selected_libmesh_ids[valid_mask], 2] = selected_jt_z[valid_mask]
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

# Consider only the time steps every 10 in the list
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

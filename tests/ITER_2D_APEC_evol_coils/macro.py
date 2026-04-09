from apec_m import APEC_M
from apec_m.model_edit import ModelEdit
from apec_m.cyclic import Cyclic
import vtk.util.numpy_support as vtknp
import numpy as np
from apec_m.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_evolution_mesh = "ITER_2D_APEC_evol_coils"
apec_evolution_mesh_notrans = "ITER_2D_APEC_coils_out_1"

evol_filepath = os.path.join(os.getcwd(), "tests", apec_evolution_mesh, apec_evolution_mesh + ".vtu")
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    "ITER_2D_APEC_coils",
    "ITER_2D_APEC_coils",
    "output",
    "ITER_2D_APEC_coils_out_1.vtu",
)


apec.read_mesh(evol_filepath)
apec.read_mesh(static_mesh_filepath, reset_arrays=False)

apec.read_excel("input_template.xlsx")

# Set Boundary conditions
apec.set_active_mesh(apec_evolution_mesh_notrans)

# Selector.select_all(apec.mesh)
# Selector.nsel_ext(apec.mesh)
# Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)
# Selector.select_all(apec.mesh)

# Write loads and BCs and save the static mesh
apec.load_mesh_into_solver()
apec.save_mesh()

# Now start the evolution analysis
apec.set_active_mesh(apec_evolution_mesh)
Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

# 1. Build source->destination element map in VACUUM_CHAMBER using toroidal slices
# Selector.select(apec.get_mesh(apec_evolution_mesh_notrans), "VACUUM_CHAMBER", 1, type="elem", mode="s")
# Selector.select(apec.get_mesh(apec_evolution_mesh_notrans), "COILS", 1, type="elem", mode="a")
# Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
jt_z_array = vtknp.vtk_to_numpy(apec.get_mesh(apec_evolution_mesh_notrans).GetCellData().GetArray("Js"))
js_array = vtknp.vtk_to_numpy(apec.mesh.GetCellData().GetArray("Js"))

js_array[:] = jt_z_array
apec.mesh.GetCellData().GetArray("Js").Modified()

for k, comp in enumerate(["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L"]):
    Selector.select(apec.mesh, comp, 1, type="elem", mode="s")
    Selector.nsle(apec.mesh)
    Selector.set_value_on_selection(apec.mesh, "CP_0", 11 + k)
    master_node_cs = vtknp.vtk_to_numpy(
        apec.mesh.GetPointData().GetArray("selection_point")
    ).nonzero()[0][0]
    Selector.select(apec.mesh, "ID", master_node_cs)
    mask = apec.constants["CONSTANTS_ID"] == 11 + k
    n_turns = apec.constants[mask]["TURNS"].values[0]
    resistance = apec.constants[mask]["R"].values[0]
    Selector.select(
        apec.get_mesh(apec_evolution_mesh_notrans), comp, 1, type="elem", mode="s"
    )
    output_mesh = apec.get_mesh(apec_evolution_mesh_notrans)
    selection_mask = vtknp.vtk_to_numpy(
        output_mesh.GetCellData().GetArray("selection_cell")
    )
    jt_array = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("Jt"))
    jt_z = jt_array[selection_mask == 1, 2][0]
    cs = apec.constants[mask]["CS"].values[0]
    current = cs * jt_z
    voltage = -resistance * current / n_turns
    Selector.set_value_on_selection(apec.mesh, "_VOLT", voltage)

for k, comp in enumerate(["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6"]):
    Selector.select(apec.mesh, comp, 1, type="elem", mode="s")
    Selector.nsle(apec.mesh)
    Selector.set_value_on_selection(apec.mesh, "CP_0", 1 + k)
    master_node_pfcs = vtknp.vtk_to_numpy(
        apec.mesh.GetPointData().GetArray("selection_point")
    ).nonzero()[0][0]
    Selector.select(apec.mesh, "ID", master_node_pfcs)
    mask = apec.constants["CONSTANTS_ID"] == 1 + k
    n_turns = apec.constants[mask]["TURNS"].values[0]
    resistance = apec.constants[mask]["R"].values[0]
    Selector.select(
        apec.get_mesh(apec_evolution_mesh_notrans), comp, 1, type="elem", mode="s"
    )
    output_mesh = apec.get_mesh(apec_evolution_mesh_notrans)
    selection_mask = vtknp.vtk_to_numpy(
        output_mesh.GetCellData().GetArray("selection_cell")
    )
    jt_array = vtknp.vtk_to_numpy(output_mesh.GetCellData().GetArray("Jt"))
    jt_z = jt_array[selection_mask == 1, 2][0]
    cs = apec.constants[mask]["CS"].values[0]
    current = cs * jt_z
    voltage = -resistance * current / n_turns
    Selector.set_value_on_selection(apec.mesh, "_VOLT", voltage)


Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.nsle(apec.mesh)
apec.algorithm.input_parameters["target_central_axis_x"] = 6.52136102486433
apec.algorithm.input_parameters["target_central_axis_y"] = 0.536537145137914
apec.algorithm.input_parameters["target_central_axis_z"] = 0.0
apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["min_distance"] = 0.1
apec.algorithm.input_parameters["cp_gain"] = 0.250
apec.algorithm.input_parameters["tolerance_value"] = 1e-3
apec.algorithm.input_parameters["bradx_value"] = -3.873415861e-08
apec.algorithm.input_parameters["brady_value"] = 5.313857037e-09
apec.algorithm.input_parameters["bverx_value"] = 1.822467692e-10
apec.algorithm.input_parameters["bvery_value"] = 1.133287835e-07
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["equilibrium_mesh"] = apec_evolution_mesh_notrans
apec.algorithm.input_parameters["static_mesh"] = apec_evolution_mesh_notrans
apec.algorithm.input_parameters["evolution_mesh"] = apec_evolution_mesh
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2
# apec.algorithm.input_parameters["plasma_current"] = 15e6
# apec.algorithm.input_parameters["beta_p"] = 6.8216e-001
apec.algorithm.input_parameters["toroidal_slices_nodes"] = Cyclic.build_toroidal_slices(
    apec.mesh, apec.mesh, "node", n_slices=1, tolerance=1e-3
)
apec.algorithm.input_parameters["toroidal_slices_elems"] = Cyclic.build_toroidal_slices(
    apec.mesh, apec.mesh, "elem", n_slices=1, tolerance=1e-3
)
apec.save_mesh()
apec.load_mesh_into_solver()
# Save the toroidal projection maps to files
# np.savetxt("toroidal_slices_nodes.txt", apec.algorithm.input_parameters["toroidal_slices_nodes"].astype(int), fmt="%d")
# np.savetxt("toroidal_slices_elems.txt", apec.algorithm.input_parameters["toroidal_slices_elems"].astype(int), fmt="%d")

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::10]))
currents = np.concatenate((currents[0:1], currents[1::10]))
betas = np.concatenate((betas[0:1], betas[1::10]))
apec.time_step["theta"] = 1.0
for step_idx in range(len(time_steps)):
    apec.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.time_step["transient_effects"] = False
        apec.time_step["time"] = 1e-10
    else:
        apec.time_step["transient_effects"] = True
    apec.algorithm.input_parameters["plasma_current"] = currents[step_idx]
    apec.algorithm.input_parameters["beta_p"] = betas[step_idx]
    apec.algorithm.evolution_solve()
    apec.write_output()

apec.save_mesh()

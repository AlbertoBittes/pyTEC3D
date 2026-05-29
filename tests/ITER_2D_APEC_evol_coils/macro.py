from apec_m import APEC_M
from apec_m.model_edit import ModelEdit
import vtk.util.numpy_support as vtknp
import numpy as np
from apec_m.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_evolution_mesh = "ITER_2D_APEC_evol_coils"

evol_filepath = os.path.join(os.getcwd(),"tests", apec_evolution_mesh, apec_evolution_mesh + ".vtu")
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_APEC_coils",
    "output",
    "ITER_2D_APEC_coils_out_1.vtu",
)


apec.read_mesh(evol_filepath)
apec.algorithm.equilibrium_mesh = ModelEdit.read_mesh(static_mesh_filepath, reset_arrays=False)

apec.read_excel("input_template.xlsx")

# Now start the evolution analysis
Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

# 1. Build source->destination element map in VACUUM_CHAMBER using toroidal slices
# Selector.select(apec.get_mesh(apec_evolution_mesh_notrans), "VACUUM_CHAMBER", 1, type="elem", mode="s")
# Selector.select(apec.get_mesh(apec_evolution_mesh_notrans), "COILS", 1, type="elem", mode="a")
# Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
jt_z_array = vtknp.vtk_to_numpy(apec.algorithm.equilibrium_mesh.GetCellData().GetArray("Js"))
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
        apec.algorithm.equilibrium_mesh, comp, 1, type="elem", mode="s"
    )
    selection_mask = vtknp.vtk_to_numpy(
        apec.algorithm.equilibrium_mesh.GetCellData().GetArray("selection_cell")
    )
    jt_array = vtknp.vtk_to_numpy(apec.algorithm.equilibrium_mesh.GetCellData().GetArray("Jt"))
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
        apec.algorithm.equilibrium_mesh, comp, 1, type="elem", mode="s"
    )
    selection_mask = vtknp.vtk_to_numpy(
        apec.algorithm.equilibrium_mesh.GetCellData().GetArray("selection_cell")
    )
    jt_array = vtknp.vtk_to_numpy(apec.algorithm.equilibrium_mesh.GetCellData().GetArray("Jt"))
    jt_z = jt_array[selection_mask == 1, 2][0]
    cs = apec.constants[mask]["CS"].values[0]
    current = cs * jt_z
    voltage = -resistance * current / n_turns
    Selector.set_value_on_selection(apec.mesh, "_VOLT", voltage)


Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.nsle(apec.mesh)
apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["tolerance_value"] = 1e-4
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))
apec.save_mesh()

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
apec.time_step["theta"] = 1.0

apec.load_mesh_into_solver()
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

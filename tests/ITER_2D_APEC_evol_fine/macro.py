from apec_m import APEC_M
from apec_m.apec_m_model import APEC_M_Model
import vtk.util.numpy_support as vtknp
import numpy as np
from apec_m.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

apec_evolution_mesh = "ITER_2D_APEC_evol_fine"
apec_evolution_mesh_notrans = "ITER_2D_APEC_fine_out_1"

evol_filepath = apec_evolution_mesh + ".vtu"
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_APEC_fine",
    "output",
    "ITER_2D_APEC_fine_out_1.vtu",
)


apec.read_model(evol_filepath)
output_apec_equilibrium = APEC_M_Model.read_model(static_mesh_filepath)
apec.algorithm.vacuum_chamber_mesh = output_apec_equilibrium.mesh

apec.model.read_excel("input_template.xlsx")

# Write loads and BCs and save the static mesh
apec.load_model_into_solver()

# Now start the evolution analysis
Selector.select_all(apec.model.mesh)
Selector.nsel_ext(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
jt_z_array = vtknp.vtk_to_numpy(output_apec_equilibrium.mesh.GetCellData().GetArray("Jt"))
js_array = vtknp.vtk_to_numpy(apec.model.mesh.GetCellData().GetArray("Js"))

js_array[:] = jt_z_array
apec.model.mesh.GetCellData().GetArray("Js").Modified()

Selector.select(apec.model.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.nsle(apec.model.mesh)
apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["tolerance_value"] = 1e-3
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
apec.model.time_step["theta"] = 1.0
for step_idx in range(len(time_steps)):
    apec.model.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.model.time_step["transient_effects"] = False
        apec.model.time_step["time"] = 1e-10
    else:
        apec.model.time_step["transient_effects"] = True
    apec.algorithm.input_parameters["plasma_current"] = currents[step_idx]
    apec.algorithm.input_parameters["beta_p"] = betas[step_idx]
    apec.algorithm.evolution_solve()
    apec.write_output()

apec.save_model()

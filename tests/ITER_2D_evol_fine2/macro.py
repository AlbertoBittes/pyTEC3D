from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.model import pyTEC3D_Model
import vtk.util.numpy_support as vtknp
import numpy as np
from pytec3d.selector import Selector
import os

# from convert_mesh import convert_mesh_from_cdb

pytec3d = pyTEC3DSolver()

pytec3d_evolution_mesh = "ITER_2D_evol_fine2"
pytec3d_evolution_mesh_notrans = "ITER_2D_fine2_out_1"

evol_filepath = os.path.join(os.getcwd(), pytec3d_evolution_mesh + ".vtu")
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_fine2",
    "output",
    "ITER_2D_fine2_out_1.vtu",
)


pytec3d.read_model(evol_filepath)
output_pytec3d_equilibrium = pyTEC3D_Model.read_model(static_mesh_filepath)
pytec3d.algorithm.vacuum_chamber_mesh = output_pytec3d_equilibrium.mesh

pytec3d.model.read_excel("input_template.xlsx")

# Write loads and BCs and save the static mesh
pytec3d.load_model_into_solver()

# Now start the evolution analysis
Selector.select_all(pytec3d.model.mesh)
Selector.nsel_ext(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
jt_z_array = vtknp.vtk_to_numpy(output_pytec3d_equilibrium.mesh.GetCellData().GetArray("Jt"))
js_array = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetCellData().GetArray("Js"))

js_array[:] = jt_z_array
pytec3d.model.mesh.GetCellData().GetArray("Js").Modified()

Selector.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.nsle(pytec3d.model.mesh)
pytec3d.algorithm.input_parameters["R0"] = 6.2
pytec3d.algorithm.input_parameters["max_picard_iterations"] = 70
pytec3d.algorithm.input_parameters["tolerance_value"] = 1e-4
pytec3d.algorithm.input_parameters["alpha_mp"] = 1.0001
pytec3d.algorithm.input_parameters["alpha_np"] = 0.83197
pytec3d.algorithm.input_parameters["AF0"] = 5.3 * 6.2

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
pytec3d.model.time_step["theta"] = 1.0
for step_idx in range(len(time_steps)):
    pytec3d.model.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        pytec3d.model.time_step["transient_effects"] = False
        pytec3d.model.time_step["time"] = 1e-10
    else:
        pytec3d.model.time_step["transient_effects"] = True
    pytec3d.algorithm.input_parameters["plasma_current"] = currents[step_idx]
    pytec3d.algorithm.input_parameters["beta_p"] = betas[step_idx]
    pytec3d.algorithm.evolution_solve()
    pytec3d.write_output()

pytec3d.save_model()

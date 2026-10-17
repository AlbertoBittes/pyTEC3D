from pytec3d import Algorithm, Model, selection, Solver
import vtk.util.numpy_support as vtknp
import numpy as np
import os

# from convert_mesh import convert_mesh_from_cdb

pytec3d = Solver()
algorithm = Algorithm(pytec3d)

pytec3d_evolution_mesh = "ITER_2D_evol_coils"

evol_filepath = pytec3d_evolution_mesh + ".vtu"
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_coils",
    "output",
    "ITER_2D_coils_out_1.vtu",
)


pytec3d.read_model(evol_filepath)
output_pytec3d_equilibrium = Model.read_model(static_mesh_filepath)
algorithm.vacuum_chamber_mesh = output_pytec3d_equilibrium.mesh

pytec3d.model.read_json("model_data.json")

# Now start the evolution analysis
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
jt_z_array = vtknp.vtk_to_numpy(output_pytec3d_equilibrium.mesh.GetCellData().GetArray("Js"))
js_array = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetCellData().GetArray("Js"))

js_array[:] = jt_z_array
pytec3d.model.mesh.GetCellData().GetArray("Js").Modified()

for k, comp in enumerate(["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L"]):
    selection.select(pytec3d.model.mesh, comp, 1, location="elem", mode="s")
    selection.nsle(pytec3d.model.mesh)
    selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 11 + k)
    master_node_cs = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetPointData().GetArray("selection_point")
    ).nonzero()[0][0]
    selection.select(pytec3d.model.mesh, "ID", master_node_cs)
    constants_id = 11 + k
    n_turns = pytec3d.model.real_constant(constants_id, "turns")
    resistance = pytec3d.model.real_constant(constants_id, "resistance")
    selection.select(
        output_pytec3d_equilibrium.mesh, comp, 1, location="elem", mode="s"
    )
    selection_mask = vtknp.vtk_to_numpy(
        output_pytec3d_equilibrium.mesh.GetCellData().GetArray("selection_cell")
    )
    jt_array = vtknp.vtk_to_numpy(output_pytec3d_equilibrium.mesh.GetCellData().GetArray("Jt"))
    jt_z = jt_array[selection_mask == 1, 2][0]
    cs = pytec3d.model.real_constant(constants_id, "cross_section")
    current = cs * jt_z
    voltage = -resistance * current / n_turns
    selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", voltage)

for k, comp in enumerate(["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6"]):
    selection.select(pytec3d.model.mesh, comp, 1, location="elem", mode="s")
    selection.nsle(pytec3d.model.mesh)
    selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1 + k)
    master_node_pfcs = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetPointData().GetArray("selection_point")
    ).nonzero()[0][0]
    selection.select(pytec3d.model.mesh, "ID", master_node_pfcs)
    constants_id = 1 + k
    n_turns = pytec3d.model.real_constant(constants_id, "turns")
    resistance = pytec3d.model.real_constant(constants_id, "resistance")
    selection.select(
        output_pytec3d_equilibrium.mesh, comp, 1, location="elem", mode="s"
    )
    selection_mask = vtknp.vtk_to_numpy(
        output_pytec3d_equilibrium.mesh.GetCellData().GetArray("selection_cell")
    )
    jt_array = vtknp.vtk_to_numpy(output_pytec3d_equilibrium.mesh.GetCellData().GetArray("Jt"))
    jt_z = jt_array[selection_mask == 1, 2][0]
    cs = pytec3d.model.real_constant(constants_id, "cross_section")
    current = cs * jt_z
    voltage = -resistance * current / n_turns
    selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", voltage)


selection.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, location="elem", mode="s")
selection.nsle(pytec3d.model.mesh)
algorithm.input_parameters["R0"] = 6.2
algorithm.input_parameters["max_picard_iterations"] = 70
algorithm.input_parameters["tolerance_value"] = 1e-4
algorithm.input_parameters["alpha_mp"] = 1.0001
algorithm.input_parameters["alpha_np"] = 0.83197
algorithm.input_parameters["AF0"] = 5.3 * 6.2

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))
pytec3d.save_model()

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
pytec3d.time_step["theta"] = 1.0

pytec3d.load_model_into_solver()
for step_idx in range(len(time_steps)):
    pytec3d.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        pytec3d.time_step["transient_effects"] = False
        pytec3d.time_step["time"] = 1e-10
    else:
        pytec3d.time_step["transient_effects"] = True
    algorithm.input_parameters["plasma_current"] = currents[step_idx]
    algorithm.input_parameters["beta_p"] = betas[step_idx]
    algorithm.evolution_solve()
    pytec3d.write_output()

pytec3d.save_model()

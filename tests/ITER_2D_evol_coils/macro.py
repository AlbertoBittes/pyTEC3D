# =====================================================================================================================
# ITER_2D_evol_coils - ITER vertical displacement (beta_p drop) in 2D with the coils as voltage-driven STRANDED COILS (circuit coupling)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The vertical displacement of ITER_2D_evol (read that macro first), with the central-solenoid modules and the poloidal
# field coils modelled as stranded coils connected to voltage sources instead of fixed current densities: each coil
# has its own circuit equation VOLT = R i + EMF, so the coil currents change during the transient as the plasma moves
# and induces voltages in the windings. It shows:
#   * one DOF coupling per coil (CP ids 1-6 for PF_1 ... PF_6, 11-16 for the CS modules), coupling VOLT and EMF;
#   * reading real constants from the model data in Python (Model.real_constant) to compute the loads;
#   * the terminal voltage of each coil set to the steady-state value that reproduces the equilibrium current:
#     VOLT = R i, with the coil current i = -(S Jt_z) / N from the equilibrium current density (the minus sign because
#     the winding direction is T = (0, 0, -1)).
#
# Geometry and mesh (ITER_2D_evol_coils.vtu): the 2D axisymmetric ITER mesh of ITER_2D_coils (x = R, y = Z),
#   17818 QUAD8 + 16 TRI6 cells; TYPE 25 on the coil cells (440), 24 on the passive structures, 23 on the plasma
#   region, 21 air. Named cell arrays CS1U ... CS3L, PFC_1 ... PFC_6, VACUUM_CHAMBER.
#
# Model data (model_data.json)
#   Element types (AV2D_nodal, option2 = 2, option3 = 1): 21 NON_CONDUCTING, 23 STRANDED_CONDUCTOR (plasma region, Js
#   imposed), 24 SOLID_CONDUCTOR (structures, eddy currents), 25 STRANDED_COIL (coils: A + VOLT + EMF).
#   Materials: mur = 1; rho = 8e-7 Ohm m for the structures.
#   Real constants: rows 1-6 (PF coils) and 11-16 (CS modules): cross_section S, turns N, radius, resistance R,
#     symmetry = 1, T = (0, 0, -1).
#   DOF couplings: 1-6 and 11-16 -> ["VOLT", "EMF"] (one terminal per coil).
#
# Boundary conditions:  A_z = 0 on the external boundary; VOLT of each coil = its steady-state voltage R i.
# Loads:  Js of the plasma region from the equilibrium (copied from the Js array of the equilibrium output).
# Time history: as ITER_2D_evol (tests/files waveforms, one point every 5, static first step at t = 1e-10 s, then
#   backward-Euler transient steps with one evolution solve each).
# Algorithm parameters: R0 = 6.2 m, AF0 = 5.3 * 6.2, alpha_mp = 1.0001, alpha_np = 0.83197, tolerance 1e-4.
#
# Results: output/ITER_2D_evol_coils_out_<n>.vtu (with VOLT, EMF of every coil), ITER_2D_evol_coils_steps.csv, .log
#
# How to run:  run ITER_2D_coils first (it writes ../ITER_2D_coils/output), then
#              cd tests/ITER_2D_evol_coils && python macro.py      (a long run)
# Documentation: "User manual - plasma equilibrium and evolution", "Theory - stranded coils (circuit coupling)".
# =====================================================================================================================


from pytec3d import Algorithm, Model, selection, Solver
import vtk.util.numpy_support as vtknp
import numpy as np
import glob
import os

# from convert_mesh import convert_mesh_from_cdb

# 1. The solver and the plasma algorithm.
pytec3d = Solver()
algorithm = Algorithm(pytec3d)

# 2. Files: the evolution mesh and the equilibrium output of ITER_2D_coils.
pytec3d_evolution_mesh = "ITER_2D_evol_coils"

evol_filepath = pytec3d_evolution_mesh + ".vtu"
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_coils",
    "output",
    "ITER_2D_coils_out_*.vtu",
)
# The equilibrium case writes its result once, after its Picard loop, as output/<name>_out_<n>.vtu (n = number of
# solves): read the most recent file of that series.
static_mesh_filepath = max(glob.glob(static_mesh_filepath), key=os.path.getmtime)


# 3. Read the evolution model; read the equilibrium as a Model and give its vacuum-chamber cells to the algorithm.
pytec3d.read_model(evol_filepath)
output_pytec3d_equilibrium = Model.read_model(static_mesh_filepath)
algorithm.vacuum_chamber_mesh = output_pytec3d_equilibrium.mesh

# 4. Model data of the evolution (coils as stranded coils, conducting structures).
pytec3d.model.read_json("model_data.json")

# Now start the evolution analysis
# 5. Boundary condition: A_z = 0 on the external boundary.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
# 6. Plasma source current density: the Js of the equilibrium copied cell by cell (the coil cells ignore Js:
#    a stranded coil takes its current from its circuit).
jt_z_array = vtknp.vtk_to_numpy(output_pytec3d_equilibrium.mesh.GetCellData().GetArray("Js"))
js_array = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetCellData().GetArray("Js"))

js_array[:] = jt_z_array
pytec3d.model.mesh.GetCellData().GetArray("Js").Modified()

# 7. Central-solenoid modules: for each module, couple VOLT and EMF of its nodes (CP id 11 + k), find the master
#    node, read N, R and S of its real constants row (11 + k), and impose the steady-state terminal voltage that
#    gives the equilibrium current of the module.
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

# 8. Poloidal field coils: the same, CP ids and real constants rows 1 + k.
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


# 9. Vacuum-chamber cells and nodes of the solver mesh (mapped by the algorithm), then the evolution parameters.
selection.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, location="elem", mode="s")
selection.nsle(pytec3d.model.mesh)
algorithm.input_parameters["R0"] = 6.2
algorithm.input_parameters["max_picard_iterations"] = 70
algorithm.input_parameters["tolerance_value"] = 1e-4
algorithm.input_parameters["alpha_mp"] = 1.0001
algorithm.input_parameters["alpha_np"] = 0.83197
algorithm.input_parameters["AF0"] = 5.3 * 6.2

# Read the numpy arrays of time, current and beta from file
# 10. Waveforms of the scenario: times, plasma current [A], poloidal beta; save the prepared model.
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))
pytec3d.save_model()

# Consider only the time steps every 5 in the list
# Keep the first point and one point every 5; backward Euler.
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
pytec3d.time_step["theta"] = 1.0

# 11. Load the model, then the time loop: static initial step, transient steps, one evolution solve each.
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

# 12. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

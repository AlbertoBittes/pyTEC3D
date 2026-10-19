# =====================================================================================================================
# ITER_2D_evol_fine2 - vertical displacement of ITER (beta_p drop) in 2D axisymmetry, starting from the ITER_2D_fine2 equilibrium
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The plasma EVOLUTION (vertical displacement) after the equilibrium computed by ITER_2D_fine2: a transient analysis in which the
# poloidal beta drops from 0.68 to 0.24 at about 15 MA (prescribed waveforms): the plasma moves up and becomes limited
# around t = 0.36 s, a vertical displacement event, while the induced (eddy) currents in the
# passive structures react to it. At each time step Algorithm.evolution_solve() iterates (Picard) on the plasma
# current density: plasma boundary from the flux, Grad-Shafranov profile scaled to the plasma current of the step,
# field solve, until PSI converges. It shows:
#   * how to start from the output of another analysis: the equilibrium mesh is read with Model.read_model() and given
#     to the algorithm as its vacuum-chamber mesh; its total current density Jt becomes the source Js of the
#     evolution (the two meshes must have the same cells, in the same order);
#   * boundary conditions and loads written AFTER load_model_into_solver(): array values can change at any time;
#   * the selection of the vacuum-chamber cells and nodes on the solver mesh, which the algorithm maps to the
#     vacuum-chamber mesh;
#   * a long time loop driven by waveforms read from text files.
#
# Geometry and mesh (ITER_2D_evol_fine2.vtu): the 2D axisymmetric ITER mesh of ITER_2D_fine2 (x = R, y = Z, R <= 100 m, |Z| <= 100 m),
#   17860 QUAD4 + 16 TRI3 cells, with the passive structures as conductors (TYPE 24).
#
# Model data (model_data.json)
#   Element types (AV2D_nodal, option2 = 1, option3 = 1): 21 NON_CONDUCTING (air), 23 STRANDED_CONDUCTOR (coils and plasma
#   region: Js imposed), 24 SOLID_CONDUCTOR (vacuum vessel and passive structures: eddy currents, no
#   applied voltage), 25 STRANDED_COIL (row not used by the cells of this mesh).
#   Materials: mur = 1; rho = 8e-7 Ohm m for the structures (MATERIAL 24).
#
# Boundary conditions:  A_z = 0 on the external boundary.
# Loads:  Js = Jt of the equilibrium (coil currents kept constant, plasma current density updated by the algorithm).
# Time history (tests/files): time_steps.txt (0 to 0.636 s, dt = 2e-4 s), plasma_current.txt (about 15 MA)
#   and beta_p.txt; the macro keeps the first point and one point every 5 (dt = 1e-3 s). Step 0 is a STATIC solve at
#   t = 1e-10 s (the initial state, consistent with the equilibrium); the following steps are transient, backward
#   Euler (theta = 1). At each step evolution_solve() iterates on the plasma shape for the prescribed plasma current
#   and beta_p, without control currents, while the eddy currents of the conducting structures are computed.
# Algorithm parameters:  R0 = 6.2 m, AF0 = 5.3 * 6.2 (R0 B_phi0, vacuum toroidal field function), alpha_mp = 1.0001,
#   alpha_np = 0.83197, tolerance 1e-4, up to 70 Picard iterations per step.
#
# Results: output/ITER_2D_evol_fine2_out_<n>.vtu after every step (A, PSI, Js, Jt, PLASMA, B, H), ITER_2D_evol_fine2_steps.csv,
#   ITER_2D_evol_fine2.log, __pytec3d__/pytec3d_iteration_diagnostics.csv (one row per Picard iteration).
#
# How to run:  run ITER_2D_fine2 first (it writes ../ITER_2D_fine2/output), then  cd tests/ITER_2D_evol_fine2 && python macro.py  (a long run).
# Documentation: "User manual - plasma equilibrium and evolution (Algorithm)", "Examples" (ITER).
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

# 2. Files: the evolution mesh (this case) and the equilibrium output of the equilibrium case.
pytec3d_evolution_mesh = "ITER_2D_evol_fine2"
pytec3d_evolution_mesh_notrans = "ITER_2D_fine2_out_*"

evol_filepath = os.path.join(os.getcwd(), pytec3d_evolution_mesh + ".vtu")
# Read the output of the equilibrium analysis
static_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_fine2",
    "output",
    "ITER_2D_fine2_out_*.vtu",
)
# The equilibrium case writes its result once, after its Picard loop, as output/<name>_out_<n>.vtu (n = number of
# solves): read the most recent file of that series.
static_mesh_filepath = max(glob.glob(static_mesh_filepath), key=os.path.getmtime)


# 3. Read the evolution model, read the equilibrium as a separate Model (not solved: only its arrays are used),
#    and give its vacuum-chamber cells (with PSI, Js and PLASMA) to the algorithm.
pytec3d.read_model(evol_filepath)
output_pytec3d_equilibrium = Model.read_model(static_mesh_filepath)
algorithm.vacuum_chamber_mesh = output_pytec3d_equilibrium.mesh

# 4. Model data of the evolution (conducting structures).
pytec3d.model.read_json("model_data.json")

# Write loads and BCs and save the static mesh
# 5. Load the model into the solver session (the arrays can still change after this point).
pytec3d.load_model_into_solver()

# Now start the evolution analysis
# 6. Boundary condition: A_z = 0 on the external boundary.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] through map ids
# 7. Source current density: the total current density Jt of the equilibrium, copied cell by cell into Js
#    (NumPy views of the VTK arrays; Modified() tells VTK the array changed).
jt_z_array = vtknp.vtk_to_numpy(output_pytec3d_equilibrium.mesh.GetCellData().GetArray("Jt"))
js_array = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetCellData().GetArray("Js"))

js_array[:] = jt_z_array
pytec3d.model.mesh.GetCellData().GetArray("Js").Modified()

# 8. Select the vacuum-chamber cells and their nodes: the region of the solver mesh mapped to the vacuum-chamber
#    mesh of the algorithm. Then the parameters of the evolution.
selection.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, location="elem", mode="s")
selection.nsle(pytec3d.model.mesh)
algorithm.input_parameters["R0"] = 6.2
algorithm.input_parameters["max_picard_iterations"] = 70
algorithm.input_parameters["tolerance_value"] = 1e-4
algorithm.input_parameters["alpha_mp"] = 1.0001
algorithm.input_parameters["alpha_np"] = 0.83197
algorithm.input_parameters["AF0"] = 5.3 * 6.2

# Read the numpy arrays of time, current and beta from file
# 9. Waveforms of the scenario: times, plasma current [A], poloidal beta.
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
# Keep the first point and one point every 5.
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
# 10. Time loop: static initial step at t = 1e-10 s, then transient backward-Euler steps; one evolution solve
#     (Picard iterations) and one output file per step.
pytec3d.time_step["theta"] = 1.0
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

# 11. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

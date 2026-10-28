# =====================================================================================================================
# ITER_3D_evol_coils_edge - ITER vertical displacement on a 2 degree 3D sector: edge elements, tree gauge, stranded coils
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The 3D evolution of the ITER vertical displacement (workflow: documentation, "Examples - ITER") with edge elements
# and the CS modules and PF coils modelled as STRANDED COILS fed by constant voltage sources: the coil currents react
# to the plasma displacement.
# It shows the TYPE of a region changed by the macro (COILS cells -> TYPE 35, stranded coils), two CP arrays (CP_0 for
# the coil terminals, CP_1 for the periodic pairs) and the coil voltages computed from the 2D equilibrium: VOLT = R i,
# with i = S Jt_z / N the equilibrium current of the coil.
#
# Geometry and mesh: ITER_3D_evol_coils_edge.vtu, the quadratic ITER_2D_carma mesh revolved over 2 degrees in 4 cells
#   (68404 HEX20 + 400 WEDGE15 on the axis, 463974 nodes), written by make_mesh.py (python make_mesh.py --sector 2
#   --divisions 4 rebuilds it, and updates the coil volumes of model_data.json).
# Model data (model_data.json): AV3D_edge: 31 NON_CONDUCTING, 33 STRANDED_CONDUCTOR (plasma region), 34
#   SOLID_CONDUCTOR (vessel and structures), 35 STRANDED_COIL (coils: A + VOLT + EMF). Coil real constants:
#   volume = volume of the WHOLE coil, symmetry = 180 (360 / 2 degrees), T = (0, 1, 0) in the cylindrical element
#   system 1 (toroidal). DOF couplings 1-6 (PF) and 11-16 (CS) -> ["VOLT", "EMF"],
#   set by the macro (add_dof_couplings).
# Boundary conditions:  cyclic symmetry theta = 0 <-> 2 degrees (A, V; edge values need no frame
#   rotation); tangential A = 0 on the outer sphere (midside nodes); tree-cotree gauge over all the cells, aware of
#   the periodic pairs; V = 0 on one node of the vessel (the lowest node id of the passive structures).
# Time history (examples/files): time_steps.txt (0 to 0.636 s), plasma_current.txt (about 15 MA), beta_p.txt;
#   the first point and one point every 5 (dt = 1e-3 s). Step 0 static at t = 1e-10 s, then transient backward-Euler
#   steps, one evolution solve (Picard iterations on the plasma) and one output file per step.
# Algorithm parameters:  R0 = 6.2 m, AF0 = 5.3 * 6.2, alpha_mp = 1.0001, alpha_np = 0.83197, tolerance 1e-4;
#   extrusion_angle and n_toroidal_divisions describe the 3D sector (toroidal maps between the 2D plasma-region mesh
#   and the 3D mesh); halo currents (AMPS on the vessel wall) from t0_halo on, with halo_delta1 = 0.04.
#
# Results: output/ITER_3D_evol_coils_edge_out_<n>.vtu after every step, ITER_3D_evol_coils_edge_steps.csv, ITER_3D_evol_coils_edge.log,
#   output/ITER_3D_evol_coils_edge_iteration_diagnostics.csv.
#
# How to run:  run ITER_2D_carma first (it writes ../ITER_2D_carma/output), then
#   cd examples/ITER_3D_evol_coils_edge && python macro.py  (a long run).
# Documentation: "User manual - plasma equilibrium and evolution, periodic boundaries", "Examples" (ITER).
# =====================================================================================================================


from pytec3d import Algorithm, mesh, selection, Solver, tree_gauge
import numpy as np
import vtk
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import glob
import os
from vtkmodules.util import numpy_support as vtknp

# 1. The solver and the plasma algorithm.
pytec3d = Solver()
algorithm = Algorithm(pytec3d)
filename = "ITER_3D_evol_coils_edge"
pytec3d_equilibrium_output_filename = "ITER_2D_carma_out_*"

# 2. Files: the 3D evolution mesh and the 2D equilibrium output; geometric parameters (tolerances, radius of the
#    far boundary, sector angle in degrees).
filepath = os.path.join(os.getcwd(), filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_carma",
    "output",
    pytec3d_equilibrium_output_filename + ".vtu",
)
# The equilibrium case writes its result once, after its Picard loop, as output/<name>_out_1.vtu (the outputs are
# numbered by step): read the most recent file of that series.
output_eq_mesh_filepath = max(glob.glob(output_eq_mesh_filepath), key=os.path.getmtime)

eps = 1e-1
rad_eps = 10
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 2.0

# # Read the evolution mesh
# 3. Read the 3D model; the 2D equilibrium gives the plasma-region mesh of the algorithm.
pytec3d.read_model(filepath)
# the full equilibrium output: the algorithm keeps a copy of its PLASMA_REGION cells only, the coil current
# densities are read from the full mesh below
equilibrium = mesh.read_mesh(output_eq_mesh_filepath, reset_arrays=False)
algorithm.plasma_region_mesh = equilibrium

# 4. Model data of the 3D evolution.
pytec3d.model.read_json("model_data.json")
# DOF couplings of the coil terminals (set here, not in the model data): the nodes of each coil (CP_0 = 1-6 PF,
# 11-16 CS) share one VOLT and one EMF.
pytec3d.model.add_dof_couplings({cp: ["VOLT", "EMF"] for cp in (1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16)})

# 5. The coil cells become stranded coils (TYPE 35) and use the cylindrical element system 1 (T toroidal); the
#    plasma region keeps the global system (the algorithm writes its Js in global components).
selection.select(pytec3d.model.mesh, "COILS", 1, location="elem", mode="s")
selection.set_value_on_selection(pytec3d.model.mesh, "TYPE", 35, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())

selection.select(pytec3d.model.mesh, "PLASMA_REGION", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, location="elem")

# 1. Set coils equilibrium voltages from the equilibrium output mesh
# 6. Coil terminals in a fresh CP_0: for each coil, couple VOLT and EMF of its nodes, read N, R and S of its real
#    constants row, and impose the steady-state voltage of the equilibrium current on the master node.
mesh.remove_array(pytec3d.model.mesh, "CP_0")
mesh.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
selection.select(pytec3d.model.mesh, "COILS", 1, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")

for k, comp in enumerate(["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L"]):
    selection.select(pytec3d.model.mesh, comp, 1, location="elem", mode="s")
    selection.nsle(pytec3d.model.mesh)
    selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 11 + k)
    master_node_cs = selection.get_selection_values(pytec3d.model.mesh, "ID")[0]
    selection.select(pytec3d.model.mesh, "ID", master_node_cs)
    constants_id = 11 + k
    n_turns = pytec3d.model.real_constant(constants_id, "turns")
    resistance = pytec3d.model.real_constant(constants_id, "resistance")
    selection.select(
        equilibrium, comp, 1, location="elem", mode="s"
    )
    jt_z = selection.get_selection_values(equilibrium, "Jt", 2, location="elem")[0]
    cs = pytec3d.model.real_constant(constants_id, "cross_section")
    current = cs * jt_z
    voltage = resistance * current / n_turns
    selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", voltage)

for k, comp in enumerate(["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6"]):
    selection.select(pytec3d.model.mesh, comp, 1, location="elem", mode="s")
    selection.nsle(pytec3d.model.mesh)
    selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1 + k)
    master_node_pfcs = selection.get_selection_values(pytec3d.model.mesh, "ID")[0]
    selection.select(pytec3d.model.mesh, "ID", master_node_pfcs)
    constants_id = 1 + k
    n_turns = pytec3d.model.real_constant(constants_id, "turns")
    resistance = pytec3d.model.real_constant(constants_id, "resistance")
    selection.select(equilibrium, comp, 1, location="elem", mode="s")
    jt_z = selection.get_selection_values(equilibrium, "Jt", 2, location="elem")[0]
    cs = pytec3d.model.real_constant(constants_id, "cross_section")
    current = cs * jt_z
    voltage = resistance * current / n_turns
    selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", voltage)


# Cyclic periodic boundaries
# 7. Cyclic periodic boundaries in CP_1: the face theta = 0 paired with theta = extr_angle, coupling A and V.
mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)
selection.select_pos(pytec3d.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_1", variables=("A", "V"))

# Far-field boundary conditions and tree gauge constraints
# 8. Far field: tangential A = 0 on the outer sphere, then the tree-cotree gauge (periodic pairs in CP_1).
selection.select_all(pytec3d.model.mesh)
selection.select_pos(
    pytec3d.model.mesh, 0, (max_rad - rad_eps, max_rad + rad_eps), coordinate_system=SphericalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 2)
tree_gauge.tree_gauge(
    pytec3d.model.mesh,
    periodic_point_coupling="CP_1",
)

# Ground Vacuum Vessel
# 9. Ground the potential V of the vacuum vessel on one node.
# the lowest node id of the passive structures
selection.select(pytec3d.model.mesh, "PASSIVE_STRUCTURES", 1, location="elem")
selection.nsle(pytec3d.model.mesh)
ground_node = int(selection.selected_ids(pytec3d.model.mesh, "node")[0])
selection.select(pytec3d.model.mesh, "ID", ground_node)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 10. Save the prepared model (overwrites the input .vtu).
pytec3d.save_model()

# Evolution parameters of the algorithm (see the header).
algorithm.input_parameters["R0"] = 6.2
algorithm.input_parameters["max_picard_iterations"] = 70
algorithm.input_parameters["tolerance_value"] = 1e-4
algorithm.input_parameters["alpha_mp"] = 1.0001
algorithm.input_parameters["alpha_np"] = 0.83197
algorithm.input_parameters["AF0"] = 5.3 * 6.2
algorithm.input_parameters["halo_delta1"] = 0.04
algorithm.input_parameters["halo_delta2"] = 0
algorithm.input_parameters["t0_halo"] = 0.41
algorithm.input_parameters["extrusion_angle"] = extr_angle
algorithm.input_parameters["n_toroidal_divisions"] = 4

# Read the numpy arrays of time, current and beta from file
# Waveforms of the scenario: times, plasma current [A], poloidal beta (first point + one every 5).
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
pytec3d.time_step["theta"] = 1.0

# Load the model, then the time loop: static initial step, transient steps, one evolution solve and one output each.
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

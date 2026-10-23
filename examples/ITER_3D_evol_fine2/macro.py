# =====================================================================================================================
# ITER_3D_evol_fine2 - ITER vertical displacement (beta_p drop) on a 1 degree 3D sector (nodal, time-integrated V)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The 3D evolution of the ITER vertical displacement (workflow: documentation, "Examples - ITER") on a 1 degree,
# 1 cell sector of the fine2 mesh, with halo currents from the limiter phase. The plasma-region mesh is read
# with mesh.read_mesh(reset_arrays=False), which keeps every array of the equilibrium output.
#
# Geometry and mesh: ITER_3D_evol_fine2.vtu, the fine2 ITER mesh revolved over 1 degree with 1 cell (17772 HEX8 +
#   104 WEDGE6), z <= 0.
#
# Model data (model_data.json)
#   Element types (option2 = 1, option3 = 0): 21 NON_CONDUCTING, 23 STRANDED_CONDUCTOR, 24 SOLID_CONDUCTOR,
#   all AV3D_nodal_timint (element_type 4).
#   Materials: mur = 1; rho = 8e-7 Ohm m for the vacuum vessel and structures (MATERIAL 24).
#   Real constants: one row per coil and rows for the other regions (as in the 2D models).
#
# Boundary conditions
#   * Axis (r < 0.1 m): A_x = A_z = 0.   * Outer sphere r = 100 m: A = 0.
#   * Cyclic symmetry: the face theta = 0 coupled to theta = 1 degrees (A and V, CP_0), the nodal frames of the
#     second face rotated by the sector angle (node_angles).
#   * V = 0 on one node of the vacuum vessel (node 9664): grounds the floating
#     potential of the conductors.
# Loads:  Js of the coils from the 2D equilibrium (toroidal component, in the cylindrical element coordinate system 1
#   of the COILS cells); the plasma current density is projected by the algorithm onto the vacuum chamber.
# Time history (examples/files): time_steps.txt (0 to 0.636 s), plasma_current.txt (about 15 MA), beta_p.txt;
#   the first point and one point every 5 (dt = 1e-3 s). Step 0 static at t = 1e-10 s, then transient backward-Euler
#   steps (theta = 1), one evolution solve (Picard iterations on the plasma) and one output file per step.
# Algorithm parameters:  R0 = 6.2 m, AF0 = 5.3 * 6.2, alpha_mp = 1.0001, alpha_np = 0.83197, tolerance 1e-4, up to
#   70 Picard iterations per step; extrusion_angle and n_toroidal_divisions describe the 3D sector, so that the
#   algorithm maps the 2D plasma-region mesh onto the toroidal stations of the 3D mesh (projection of Js, toroidal
#   average of PSI); halo currents (AMPS on the vacuum-vessel wall nodes) from t0_halo on, with halo_delta1 = 0.04.
#
# Results: output/ITER_3D_evol_fine2_out_<n>.vtu after every step (A, V, Js, Jt, PLASMA, B, H), ITER_3D_evol_fine2_steps.csv, ITER_3D_evol_fine2.log,
#   output/ITER_3D_evol_fine2_iteration_diagnostics.csv.
#
# How to run:  run ITER_2D_fine2 first (it writes ../ITER_2D_fine2/output), then  cd examples/ITER_3D_evol_fine2 && python macro.py  (a long run).
# Documentation: "User manual - plasma equilibrium and evolution, periodic boundaries", "Examples" (ITER).
# =====================================================================================================================


from pytec3d import Algorithm, mesh, selection, Solver
import vtk.util.numpy_support as vtknp
import numpy as np
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import glob
import os
import vtk

# from convert_mesh import convert_mesh_from_cdb

# 1. The solver and the plasma algorithm.
pytec3d = Solver()
algorithm = Algorithm(pytec3d)

filename = "ITER_3D_evol_fine2"
pytec3d_equilibrium_output_filename = "ITER_2D_fine2_out_*"

# 2. Files: the 3D evolution mesh and the 2D equilibrium output; geometric parameters (tolerances, radius of the
#    far boundary, sector angle in degrees).
filepath = os.path.join(os.getcwd(), filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_fine2",
    "output",
    pytec3d_equilibrium_output_filename + ".vtu",
)
# The equilibrium case writes its result once, after its Picard loop, as output/<name>_out_1.vtu (the outputs are
# numbered by step): read the most recent file of that series.
output_eq_mesh_filepath = max(glob.glob(output_eq_mesh_filepath), key=os.path.getmtime)

eps = 1e-1
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 1.0

# # Read the evolution mesh
# 3. Read the 3D model; the 2D equilibrium gives the plasma-region mesh of the algorithm.
pytec3d.read_model(filepath)
# the full equilibrium output: the algorithm keeps a copy of its PLASMA_REGION cells only, the coil current
# densities are read from the full mesh below
equilibrium = mesh.read_mesh(output_eq_mesh_filepath, reset_arrays=False)
algorithm.plasma_region_mesh = equilibrium

# 4. Model data of the 3D evolution.
pytec3d.model.read_json("model_data.json")

# # Set Boundary conditions
# Axis boundary conditions
# 5. Axis (r < eps): A_x = A_z = 0.
selection.select_pos(
    pytec3d.model.mesh, 0, (-eps, +eps), coordinate_system=YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 2)

# Outer boundary far field zero boundary conditions
# 6. Far field: A = 0 on the outer sphere r = 100 m (spherical radius).
selection.select_pos(
    pytec3d.model.mesh, 0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 2)

# Cyclic periodic boundaries
# 7. Cyclic periodic boundaries: fresh CP_0, the face theta = 0 paired with its image at theta = extr_angle.
mesh.remove_array(pytec3d.model.mesh, "CP_0")
mesh.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
selection.select_pos(pytec3d.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_0", variables=("A", "V"))

# rotate nodal coordinate system to align with periodic plane
# 8. Rotate the nodal frames of the second face by the sector angle about y.
selection.select_pos(pytec3d.model.mesh, 1, (extr_angle - eps, extr_angle + eps), coordinate_system=YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", extr_angle, 2)

# Grounding vacuum vessel
# 9. Ground the vacuum vessel (one node of the structures).
selection.select(pytec3d.model.mesh, "ID", 9664)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

selection.select(pytec3d.model.mesh, "PLASMA_REGION", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, location="elem")

# Apply equilibrium coils current densities
# 10. Coil current densities of the equilibrium: toroidal component (1) of Js in the cylindrical element system 1.
selection.select(pytec3d.model.mesh, "COILS", 1, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")
coils = ["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
         "CS3L", "CS3U", "CS2L", "CS2U", "CS1L", "CS1U"]
for coil in coils:
    selection.select(equilibrium, coil, 1, location="elem", mode="s")
    js_z = selection.get_selection_values(equilibrium, "Js", 2, location="elem")[0]
    selection.select(pytec3d.model.mesh, coil, 1, location="elem", mode="s")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", js_z, 1, location="elem")

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
algorithm.input_parameters["t0_halo"] = 0.37
algorithm.input_parameters["extrusion_angle"] = extr_angle
algorithm.input_parameters["n_toroidal_divisions"] = 1

# Read the numpy arrays of time, current and beta from file
# Waveforms of the scenario: times, plasma current [A], poloidal beta.
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
# Keep the first point and one point every 5 (above); backward Euler.
pytec3d.time_step["theta"] = 1.0

pytec3d.save_model()

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


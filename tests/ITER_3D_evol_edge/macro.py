# =====================================================================================================================
# ITER_3D_evol_edge - ITER vertical displacement (beta_p drop) on a 2 degree sector with EDGE elements and a tree gauge
# =====================================================================================================================
#
# What this example shows
# -----------------------
# ITER_3D_evol (read that macro first) with Nedelec edge elements for A on the quadratic version of the 2 degree mesh.
# Compared with the nodal case: no axis condition and no frame rotation (edge DOFs are circulations along the edges),
# the far field is a tangential condition on the midside nodes, and the gauge is a tree-cotree decomposition built
# by tree_gauge(), aware of the periodic pairs.
#
# Geometry and mesh: ITER_3D_evol_edge.vtu, the quadratic ITER_3D_evol mesh (22272 HEX20 + 208 WEDGE15), 2 degrees,
#   2 cells along the torus.
#
# Model data (model_data.json)
#   Element types (option2 = 1: order of the Lagrange V, option3 = 0): 31 NON_CONDUCTING (air), 33 STRANDED_CONDUCTOR (coils and plasma region:
#   Js imposed), 34 SOLID_CONDUCTOR (vacuum vessel and structures: A + time-integrated V, nodal load AMPS =
#   halo currents), all AV3D_edge_timint (element_type 5).
#   Materials: mur = 1; rho = 8e-7 Ohm m for the vacuum vessel and structures (MATERIAL 24).
#   Real constants: one row per coil and rows for the other regions (as in the 2D models).
#
# Boundary conditions
#   * Cyclic symmetry: the face theta = 0 coupled to theta = 2 degrees (A, V, CP_0); edge values need no frame rotation.
#   * Outer sphere (r = 100 +- 10 m): tangential A = 0 (component 2 of _A on the midside nodes).
#   * Tree-cotree gauge over the edges of all the cells, grown from the nodes on the sphere r = 20 m, aware of the periodic pairs.
#   * V = 0 on node 16170 of the vacuum vessel.
# Loads:  Js of the coils from the 2D equilibrium (toroidal component, in the cylindrical element coordinate system 1
#   of the COILS cells); the plasma current density is projected by the algorithm onto the vacuum chamber.
# Time history (tests/files): time_steps.txt (0 to 0.636 s), plasma_current.txt (about 15 MA), beta_p.txt;
#   the first point and one point every 5 (dt = 1e-3 s). Step 0 static at t = 1e-10 s, then transient backward-Euler
#   steps (theta = 1), one evolution solve (Picard iterations on the plasma) and one output file per step.
# Algorithm parameters:  R0 = 6.2 m, AF0 = 5.3 * 6.2, alpha_mp = 1.0001, alpha_np = 0.83197, tolerance 1e-4, up to
#   70 Picard iterations per step; extrusion_angle and n_toroidal_divisions describe the 3D sector, so that the
#   algorithm maps the 2D vacuum-chamber mesh onto the toroidal stations of the 3D mesh (projection of Js, toroidal
#   average of PSI); halo currents (AMPS on the vacuum-vessel wall nodes) from t0_halo on, with halo_delta1 = 0.04.
#
# Results: output/ITER_3D_evol_edge_out_<n>.vtu after every step (A, V, Js, Jt, PLASMA, B, H), ITER_3D_evol_edge_steps.csv, ITER_3D_evol_edge.log,
#   __pytec3d__/pytec3d_iteration_diagnostics.csv.
#
# How to run:  run ITER_2D first (it writes ../ITER_2D/output), then  cd tests/ITER_3D_evol_edge && python macro.py  (a long run).
# Documentation: "User manual - plasma equilibrium and evolution, periodic boundaries", "Examples" (ITER).
# =====================================================================================================================


from pytec3d import Algorithm, mesh, Model, selection, Solver, tree_gauge
import numpy as np
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import glob
import os
import vtk
from vtkmodules.util import numpy_support as vtknp
# from convert_mesh import convert_mesh_from_cdb

# 1. The solver and the plasma algorithm.
pytec3d = Solver()
algorithm = Algorithm(pytec3d)
filename = "ITER_3D_evol_edge"
pytec3d_equilibrium_output_filename = "ITER_2D_out_*"

# 2. Files: the 3D evolution mesh and the 2D equilibrium output; geometric parameters (tolerances, radius of the
#    far boundary, sector angle in degrees).
filepath = os.path.join(os.getcwd(), filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D",
    "output",
    pytec3d_equilibrium_output_filename + ".vtu",
)
# The equilibrium case writes its result once, after its Picard loop, as output/<name>_out_<n>.vtu (n = number of
# solves): read the most recent file of that series.
output_eq_mesh_filepath = max(glob.glob(output_eq_mesh_filepath), key=os.path.getmtime)

eps = 1e-1
rad_eps = 10
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 2.0

# # Read the evolution mesh
# 3. Read the 3D model; the 2D equilibrium gives the vacuum-chamber mesh of the algorithm.
pytec3d.read_model(filepath)
output_pytec3d_equilibrium = Model.read_model(output_eq_mesh_filepath)
algorithm.vacuum_chamber_mesh = output_pytec3d_equilibrium.mesh

# 4. Model data of the 3D evolution.
pytec3d.model.read_json("model_data.json")

# Cyclic periodic boundaries
# 5. Cyclic periodic boundaries: fresh CP_0, the face theta = 0 paired with its image at theta = extr_angle.
mesh.remove_array(pytec3d.model.mesh, "CP_0")
mesh.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
selection.select_pos(pytec3d.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_0", variables=("A", "V"))

# BCs: Outer boundary far field zero boundary conditions + Tree gauge
# 6. Far field: tangential A = 0 on the outer sphere (midside nodes; all the cells in the global system), then the
#    tree-cotree gauge over the edges of all the cells, aware of the periodic pairs.
selection.select_all(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, location="elem")
selection.select_pos(
    pytec3d.model.mesh, 0, (max_rad - rad_eps, max_rad + rad_eps), coordinate_system=SphericalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 2)
selection.select_pos(
    pytec3d.model.mesh, 0, (20 - 0.1, 20 + 0.1), coordinate_system=SphericalCoordinateSystem()
)
tree_gauge.tree_gauge(
    pytec3d.model.mesh,
    periodic_point_coupling="CP_0",
)

# Ground Vacuum Vessel
# 7. Ground the vacuum vessel (see the header for the field this writes).
selection.select(pytec3d.model.mesh, "ID", 16170)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

selection.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, location="elem")

# 1. Build source->destination element map using toroidal slices
# 8. Coil current densities of the equilibrium: toroidal component (1) of Js in the cylindrical element system 1.
selection.select(pytec3d.model.mesh, "COILS", 1, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")
coils = ["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
         "CS3L", "CS3U", "CS2L", "CS2U", "CS1L", "CS1U"]
for coil in coils:
    selection.select(output_pytec3d_equilibrium.mesh, coil, 1, location="elem", mode="s")
    js_z = selection.get_selection_values(output_pytec3d_equilibrium.mesh, "Js", 2, location="elem")[0]
    selection.select(pytec3d.model.mesh, coil, 1, location="elem", mode="s")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", js_z, 1, location="elem")

# Evolution parameters of the algorithm (see the header).
algorithm.input_parameters["R0"] = 6.2
algorithm.input_parameters["max_picard_iterations"] = 70
algorithm.input_parameters["tolerance_value"] = 1e-4
algorithm.input_parameters["alpha_mp"] = 1.0001
algorithm.input_parameters["alpha_np"] = 0.83197
algorithm.input_parameters["AF0"] = 5.3 * 6.2
algorithm.input_parameters["halo_delta1"] = 0.04
algorithm.input_parameters["halo_delta2"] = 0
algorithm.input_parameters["t0_halo"] = 0.38
algorithm.input_parameters["extrusion_angle"] = extr_angle
algorithm.input_parameters["n_toroidal_divisions"] = 2

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

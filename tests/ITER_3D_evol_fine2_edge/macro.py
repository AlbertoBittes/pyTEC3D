from pytec3d import Algorithm, mesh, selection, Solver, tree_gauge
import numpy as np
import vtk
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
from vtkmodules.util import numpy_support as vtknp
# from convert_mesh import convert_mesh_from_cdb

pytec3d = Solver()
algorithm = Algorithm(pytec3d)
filename = "ITER_3D_evol_fine2_edge"
pytec3d_equilibrium_output_filename = "ITER_2D_fine2_out_1"

filepath = os.path.join(os.getcwd(), filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_fine2",
    "output",
    pytec3d_equilibrium_output_filename + ".vtu",
)

eps = 1e-1
rad_eps = 10
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 1.0

# # Read the evolution mesh
pytec3d.read_model(filepath)
algorithm.vacuum_chamber_mesh = mesh.read_mesh(output_eq_mesh_filepath, reset_arrays=False)

pytec3d.model.read_json("model_data.json")

# Cyclic periodic boundaries
mesh.remove_array(pytec3d.model.mesh, "CP_0")
mesh.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
selection.select_pos(pytec3d.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_0", variables=("A", "VOLT"))

# Grounding vacuum vessel
selection.select(pytec3d.model.mesh, "ID", 9664)
selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

# BCs: Outer boundary far field zero boundary conditions + Tree gauge
selection.select_all(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, location="elem")
selection.select_pos(
    pytec3d.model.mesh, 0, (max_rad - rad_eps, max_rad + rad_eps), coordinate_system=SphericalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0.0, 2)
tree_gauge.tree_gauge(
    pytec3d.model.mesh,
    periodic_point_coupling="CP_0",
)

selection.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, location="elem")

# Apply equilibrium coils current densities
selection.select(pytec3d.model.mesh, "COILS", 1, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")
coils = ["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
         "CS3L", "CS3U", "CS2L", "CS2U", "CS1L", "CS1U"]
for coil in coils:
    selection.select(algorithm.vacuum_chamber_mesh, coil, 1, location="elem", mode="s")
    js_z = selection.get_selection_values(algorithm.vacuum_chamber_mesh, "Js", 2, location="elem")[0]
    selection.select(pytec3d.model.mesh, coil, 1, location="elem", mode="s")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", js_z, 1, location="elem")

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
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
pytec3d.time_step["theta"] = 1.0

pytec3d.save_model()

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

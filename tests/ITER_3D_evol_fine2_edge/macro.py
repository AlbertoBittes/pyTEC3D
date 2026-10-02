from pytec3d_solver.solver import pyTEC3DSolver
import numpy as np
import vtk
from pytec3d.selector import Selector
from pytec3d.model_edit import ModelEdit
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import pandas as pd
from vtkmodules.util import numpy_support as vtknp
from pytec3d.tree_gauge import TreeGauge
# from convert_mesh import convert_mesh_from_cdb

pytec3d = pyTEC3DSolver()
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
pytec3d.algorithm.vacuum_chamber_mesh = ModelEdit.read_mesh(output_eq_mesh_filepath, reset_arrays=False)

pytec3d.model.read_excel("input_template.xlsx")

# Cyclic periodic boundaries
ModelEdit.remove_array(pytec3d.model.mesh, "CP_0")
ModelEdit.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
Selector.select_pos(pytec3d.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
pytec3d.model.dof_couplings = pd.concat([pytec3d.model.dof_couplings, ModelEdit.periodic_boundaries(pytec3d.model.mesh, YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_0", variables=("Az", "VOLT"))])

# Grounding vacuum vessel
Selector.select(pytec3d.model.mesh, "ID", 9664)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

# BCs: Outer boundary far field zero boundary conditions + Tree gauge
Selector.select_all(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, type="elem")
Selector.select_pos(
    pytec3d.model.mesh, 0, (max_rad - rad_eps, max_rad + rad_eps), coordinate_system=SphericalCoordinateSystem()
)
Selector.set_value_on_selection(pytec3d.model.mesh, "_Az", 0.0)
TreeGauge.tree_gauge(
    pytec3d.model.mesh,
    periodic_point_coupling="CP_0",
)

Selector.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, type="elem")
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 0, type="elem")

# Apply equilibrium coils current densities
Selector.select(pytec3d.model.mesh, "COILS", 1, type="elem")
pytec3d.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")
coils = ["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
         "CS3L", "CS3U", "CS2L", "CS2U", "CS1L", "CS1U"]
for coil in coils:
    Selector.select(pytec3d.algorithm.vacuum_chamber_mesh, coil, 1, type="elem", mode="s")
    js_z = Selector.get_selection_values(pytec3d.algorithm.vacuum_chamber_mesh, "Js", 2, type="elem")[0]
    Selector.select(pytec3d.model.mesh, coil, 1, type="elem", mode="s")
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", js_z, 1, type="elem")

pytec3d.algorithm.input_parameters["R0"] = 6.2
pytec3d.algorithm.input_parameters["max_picard_iterations"] = 70
pytec3d.algorithm.input_parameters["tolerance_value"] = 1e-4
pytec3d.algorithm.input_parameters["alpha_mp"] = 1.0001
pytec3d.algorithm.input_parameters["alpha_np"] = 0.83197
pytec3d.algorithm.input_parameters["AF0"] = 5.3 * 6.2
pytec3d.algorithm.input_parameters["halo_delta1"] = 0.04
pytec3d.algorithm.input_parameters["halo_delta2"] = 0
pytec3d.algorithm.input_parameters["t0_halo"] = 0.37
pytec3d.algorithm.input_parameters["extrusion_angle"] = extr_angle
pytec3d.algorithm.input_parameters["n_toroidal_divisions"] = 1

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
pytec3d.model.time_step["theta"] = 1.0

pytec3d.save_model()

pytec3d.load_model_into_solver()

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

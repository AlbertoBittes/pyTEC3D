from apec_m_solver import APEC_M
import numpy as np
import vtk
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m_coordinate_system.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import pandas as pd
from vtkmodules.util import numpy_support as vtknp
from apec_m.tree_gauge import TreeGauge
# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()
filename = "ITER_3D_APEC_evol_fine2_edge"
apec_equilibrium_output_filename = "ITER_2D_APEC_fine2_out_1"

filepath = os.path.join(os.getcwd(), filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_APEC_fine2",
    "output",
    apec_equilibrium_output_filename + ".vtu",
)

eps = 1e-1
rad_eps = 10
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 1.0

# # Read the evolution mesh
apec.read_model(filepath)
apec.algorithm.vacuum_chamber_mesh = ModelEdit.read_mesh(output_eq_mesh_filepath, reset_arrays=False)

apec.model.read_excel("input_template.xlsx")

# Cyclic periodic boundaries
ModelEdit.remove_array(apec.model.mesh, "CP_0")
ModelEdit.add_array(apec.model.mesh, "CP_0", dtype=vtk.VTK_INT)
Selector.select_pos(apec.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
apec.model.dof_couplings = pd.concat([apec.model.dof_couplings, ModelEdit.periodic_boundaries(apec.model.mesh, YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_0", variables=("Az", "VOLT"))])

# Grounding vacuum vessel
Selector.select(apec.model.mesh, "ID", 9664)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

# BCs: Outer boundary far field zero boundary conditions + Tree gauge
Selector.select_all(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 0, type="elem")
Selector.select_pos(
    apec.model.mesh, 0, (max_rad - rad_eps, max_rad + rad_eps), coordinate_system=SphericalCoordinateSystem()
)
Selector.set_value_on_selection(apec.model.mesh, "_Az", 0.0)
TreeGauge.tree_gauge(
    apec.model.mesh,
    periodic_point_coupling="CP_0",
)

Selector.select(apec.model.mesh, "VACUUM_CHAMBER", 1, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 0, type="elem")

# Apply equilibrium coils current densities
Selector.select(apec.model.mesh, "COILS", 1, type="elem")
apec.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")
coils = ["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
         "CS3L", "CS3U", "CS2L", "CS2U", "CS1L", "CS1U"]
for coil in coils:
    Selector.select(apec.algorithm.vacuum_chamber_mesh, coil, 1, type="elem", mode="s")
    js_z = Selector.get_selection_values(apec.algorithm.vacuum_chamber_mesh, "Js", 2, type="elem")[0]
    Selector.select(apec.model.mesh, coil, 1, type="elem", mode="s")
    Selector.set_value_on_selection(apec.model.mesh, "Js", js_z, 1, type="elem")

apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["tolerance_value"] = 1e-4
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2
apec.algorithm.input_parameters["halo_delta1"] = 0.04
apec.algorithm.input_parameters["halo_delta2"] = 0
apec.algorithm.input_parameters["t0_halo"] = 0.37
apec.algorithm.input_parameters["extrusion_angle"] = extr_angle
apec.algorithm.input_parameters["n_toroidal_divisions"] = 1

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
apec.model.time_step["theta"] = 1.0

apec.save_model()

apec.load_model_into_solver()

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

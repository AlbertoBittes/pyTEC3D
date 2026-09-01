from apec_m import APEC_M
import numpy as np
import vtk
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.tree_gauge import TreeGauge
from apec_m.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import pandas as pd
from vtkmodules.util import numpy_support as vtknp
# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()
filename = "ITER_3D_APEC_evol_coils_edge"
apec_equilibrium_output_filename = "ITER_2D_APEC_out_1"

filepath = os.path.join(os.getcwd(), filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "..",
    "ITER_2D_APEC",
    "output",
    apec_equilibrium_output_filename + ".vtu",
)

eps = 1e-1
rad_eps = 10
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 2.0

# # Read the evolution mesh
apec.read_model(filepath)
apec.algorithm.equilibrium_mesh = ModelEdit.read_mesh(output_eq_mesh_filepath, reset_arrays=False)

apec.model.read_excel("input_template.xlsx")

Selector.select(apec.model.mesh, "COILS", 1, type="elem", mode="s")
Selector.set_value_on_selection(apec.model.mesh, "TYPE", 35, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")
apec.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()

Selector.select(apec.model.mesh, "VACUUM_CHAMBER", 1, type="elem")
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 0, type="elem")

# 1. Set coils equilibrium voltages from the equilibrium output mesh
ModelEdit.remove_array(apec.model.mesh, "CP_0")
ModelEdit.add_array(apec.model.mesh, "CP_0", dtype=vtk.VTK_INT)
Selector.select(apec.model.mesh, "COILS", 1, type="elem")
apec.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")

for k, comp in enumerate(["CS3U", "CS2U", "CS1U", "CS1L", "CS2L", "CS3L"]):
    Selector.select(apec.model.mesh, comp, 1, type="elem", mode="s")
    Selector.nsle(apec.model.mesh)
    Selector.set_value_on_selection(apec.model.mesh, "CP_0", 11 + k)
    master_node_cs = Selector.get_selection_values(apec.model.mesh, "ID")[0]
    Selector.select(apec.model.mesh, "ID", master_node_cs)
    mask = apec.model.constants["CONSTANTS_ID"] == 11 + k
    n_turns = apec.model.constants[mask]["TURNS"].values[0]
    resistance = apec.model.constants[mask]["R"].values[0]
    Selector.select(
        apec.algorithm.equilibrium_mesh, comp, 1, type="elem", mode="s"
    )
    jt_z = Selector.get_selection_values(apec.algorithm.equilibrium_mesh, "Jt", 2, type="elem")[0]
    cs = apec.model.constants[mask]["CS"].values[0]
    current = cs * jt_z
    voltage = resistance * current / n_turns
    Selector.set_value_on_selection(apec.model.mesh, "_VOLT", voltage)

for k, comp in enumerate(["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6"]):
    Selector.select(apec.model.mesh, comp, 1, type="elem", mode="s")
    Selector.nsle(apec.model.mesh)
    Selector.set_value_on_selection(apec.model.mesh, "CP_0", 1 + k)
    master_node_pfcs = Selector.get_selection_values(apec.model.mesh, "ID")[0]
    Selector.select(apec.model.mesh, "ID", master_node_pfcs)
    mask = apec.model.constants["CONSTANTS_ID"] == 1 + k
    n_turns = apec.model.constants[mask]["TURNS"].values[0]
    resistance = apec.model.constants[mask]["R"].values[0]
    Selector.select(apec.algorithm.equilibrium_mesh, comp, 1, type="elem", mode="s")
    jt_z = Selector.get_selection_values(apec.algorithm.equilibrium_mesh, "Jt", 2, type="elem")[0]
    cs = apec.model.constants[mask]["CS"].values[0]
    current = cs * jt_z
    voltage = resistance * current / n_turns
    Selector.set_value_on_selection(apec.model.mesh, "_VOLT", voltage)


# Cyclic periodic boundaries
ModelEdit.remove_array(apec.model.mesh, "CP_1")
ModelEdit.add_array(apec.model.mesh, "CP_1", dtype=vtk.VTK_INT)
Selector.select_pos(apec.model.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
apec.model.dof_couplings = pd.concat([apec.model.dof_couplings, ModelEdit.periodic_boundaries(apec.model.mesh, YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_1", variables=("Az", "V"))])

# Far-field boundary conditions and tree gauge constraints
Selector.select_all(apec.model.mesh)
Selector.select_pos(
    apec.model.mesh, 0, (max_rad - rad_eps, max_rad + rad_eps), coordinate_system=SphericalCoordinateSystem()
)
Selector.set_value_on_selection(apec.model.mesh, "_Az", 0.0)
TreeGauge.tree_gauge(
    apec.model.mesh,
    periodic_point_coupling="CP_1",
)

# Ground Vacuum Vessel
Selector.select(apec.model.mesh, "ID", 16170)
Selector.set_value_on_selection(apec.model.mesh, "_V", 0)

apec.save_model()

apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["tolerance_value"] = 1e-4
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2
apec.algorithm.input_parameters["halo_delta1"] = 0.04
apec.algorithm.input_parameters["halo_delta2"] = 0
apec.algorithm.input_parameters["t0_halo"] = 0.41
apec.algorithm.input_parameters["extrusion_angle"] = extr_angle
apec.algorithm.input_parameters["n_toroidal_divisions"] = 2

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "time_steps.txt"))
currents = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "plasma_current.txt"))
betas = np.loadtxt(os.path.join(os.getcwd(), "..", "files", "beta_p.txt"))

# Consider only the time steps every 5 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::5]))
currents = np.concatenate((currents[0:1], currents[1::5]))
betas = np.concatenate((betas[0:1], betas[1::5]))
apec.model.time_step["theta"] = 1.0

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

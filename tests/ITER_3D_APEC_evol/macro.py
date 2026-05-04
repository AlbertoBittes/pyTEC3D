from apec_m import APEC_M
import vtk.util.numpy_support as vtknp
import numpy as np
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.cyclic import Cyclic
from apec_m.coordinate_system import CylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import pandas as pd

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

filename = "ITER_3D_APEC_evol"
apec_equilibrium_output_filename = "ITER_2D_APEC_out_1"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")
# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    "ITER_2D_APEC",
    "output",
    apec_equilibrium_output_filename + ".vtu",
)

eps = 1e-1
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 2.0

# # Read the evolution mesh
apec.read_mesh(filepath)
apec.algorithm._equilibrium_mesh = ModelEdit.read_mesh(output_eq_mesh_filepath, reset_arrays=False)

apec.read_excel("input_template.xlsx")

# # Set Boundary conditions
# Axis boundary conditions
Selector.select_pos(
    apec.mesh, 0, (-eps, +eps), coordinate_system=CylindricalCoordinateSystem()
)
Selector.set_value_on_selection(apec.mesh, "A", 0.0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0.0, 2)

# Outer boundary far field zero boundary conditions
Selector.select_pos(
    apec.mesh, 0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem()
)
Selector.set_value_on_selection(apec.mesh, "A", 0.0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0.0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0.0, 2)

# Cyclic periodic boundaries
Selector.select_pos(apec.mesh, 1, (-eps, eps), coordinate_system=CylindricalCoordinateSystem())
Selector.select_pos(
    apec.mesh, 0, (-eps, +eps), coordinate_system=CylindricalCoordinateSystem(), mode="u"
)
Selector.select_pos(
    apec.mesh, 0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem(), mode="u"
)
apec.dof_couplings = pd.concat([apec.dof_couplings, Cyclic.cpcyc(apec.mesh, CylindricalCoordinateSystem(), extr_angle, array_name="CP_0", variables=("A", "VOLT"))])

# rotate angles
Selector.select_pos(apec.mesh, 1, (extr_angle - eps, extr_angle + eps), coordinate_system=CylindricalCoordinateSystem())
Selector.select_pos(
    apec.mesh, 0, (-eps, +eps), coordinate_system=CylindricalCoordinateSystem(), mode="u"
)
Selector.select_pos(
    apec.mesh, 0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem(), mode="u"
)
Selector.set_value_on_selection(apec.mesh, "angles", extr_angle, 2)

Selector.select(apec.mesh, "ID", 16170)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

# 1. Build source->destination element map using toroidal slices
Selector.select(apec.algorithm._equilibrium_mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.select_all(apec.mesh, type="elem")

apec.algorithm.toroidal_projection_current()

apec.save_mesh()


apec.save_mesh(static_mesh_filepath)
apec.read_mesh(static_mesh_filepath, False)
Selector.select(apec.mesh, "PASSIVE_STRUCTURES", 1, type="elem", mode="s")
Selector.set_value_on_selection(apec.mesh, "MATERIAL", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "TYPE", 33, type="elem")
apec.write_loads_bcs()
apec.save_mesh()

apec.read_mesh(eq_filepath)
Selector.select(apec.mesh, "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.nsle(apec.mesh)
apec.set_active_mesh(filename)
Selector.select_all(apec.mesh)
apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["min_distance"] = 0.5
apec.algorithm.input_parameters["tolerance_value"] = 1e-3
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2
# apec.apec_settings["plasma_current"] = 15e6
# apec.apec_settings["beta_p"] = 6.8216e-001
apec.algorithm.input_parameters["toroidal_slices_nodes"] = Cyclic.build_toroidal_slices(
    apec.get_mesh(equilibrium_filename), apec.mesh, "node", n_slices=3, tolerance=1e-3
)
apec.algorithm.input_parameters["toroidal_slices_elems"] = Cyclic.build_toroidal_slices(
    apec.get_mesh(equilibrium_filename), apec.mesh, "elem", n_slices=2, tolerance=1e-3
)

# Save the toroidal projection maps to files
np.savetxt("toroidal_slices_nodes.txt", apec.apec_settings["toroidal_slices_nodes"].astype(int), fmt="%d")
np.savetxt("toroidal_slices_elems.txt", apec.apec_settings["toroidal_slices_elems"].astype(int), fmt="%d")

# Read the numpy arrays of time, current and beta from file
time_steps = np.loadtxt("time_steps.txt")
currents = np.loadtxt("plasma_current.txt")
betas = np.loadtxt("beta_p.txt")

# Consider only the time steps every 10 in the list
time_steps = np.concatenate((time_steps[0:1], time_steps[1::10]))
currents = np.concatenate((currents[0:1], currents[1::10]))
betas = np.concatenate((betas[0:1], betas[1::10]))

for step_idx in range(len(time_steps)):
    apec.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.time_step["transient_effects"] = 0
        apec.time_step["time"] = 1e-10
    else:
        apec.time_step["transient_effects"] = 1
    apec.algorithm.input_parameters["plasma_current"] = currents[step_idx]
    apec.algorithm.input_parameters["beta_p"] = betas[step_idx]


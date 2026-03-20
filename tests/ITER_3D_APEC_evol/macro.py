from apec_m import APEC_M
import vtk.util.numpy_support as vtknp
import numpy as np
from selector import Selector
from model_edit import ModelEdit
from cyclic import Cyclic
from coordinate_system import CylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import pandas as pd

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

filename = "ITER_3D_APEC_evol"
equilibrium_filename = "ITER_3D_APEC_evol_equilibrium"
static_filename = "ITER_3D_APEC_evol_static"
apec_equilibrium_output_filename = "ITER_2D_APEC_out_0"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")
eq_filepath = os.path.join(os.getcwd(), "tests", filename, equilibrium_filename + ".vtu")
static_mesh_filepath = os.path.join(os.getcwd(), "tests", filename, static_filename + ".vtu")

# Read the output of the equilibrium analysis
output_eq_mesh_filepath = os.path.join(
    os.getcwd(),
    "tests",
    "ITER_2D_APEC",
    "ITER_2D_APEC",
    "output",
    apec_equilibrium_output_filename + ".pvtu",
)

eps = 1e-1
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 2.0

# # Read the evolution mesh
apec.read_mesh(filepath)
apec.read_mesh(output_eq_mesh_filepath, False)
apec.set_active_mesh(filename)

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
Selector.select(apec.get_mesh(apec_equilibrium_output_filename), "VACUUM_CHAMBER", 1, type="elem", mode="s")
Selector.select(apec.get_mesh(apec_equilibrium_output_filename), "COILS", 1, type="elem", mode="a")
Selector.select_all(apec.mesh, type="elem")

elem_id_map = Cyclic.build_toroidal_slices(
    apec.get_mesh(apec_equilibrium_output_filename),
    apec.mesh,
    location="elem",
    n_slices=2,
    tolerance=1e-3,
)

# 2. Transfer Jt_z from source (equilibrium output) to destination Js[:,2] for all slices
jt_z_array = vtknp.vtk_to_numpy(apec.get_mesh(apec_equilibrium_output_filename).GetCellData().GetArray("Jt_z"))
js_array = vtknp.vtk_to_numpy(apec.mesh.GetCellData().GetArray("Js"))

src_elem_ids = elem_id_map[:, 0]
for slice_idx in range(1, elem_id_map.shape[1]):
    dst_elem_ids = elem_id_map[:, slice_idx]
    js_array[dst_elem_ids, 1] = jt_z_array[src_elem_ids]
apec.mesh.GetCellData().GetArray("Js").Modified()

ModelEdit.map_vector_field_to_global_components(apec.mesh, "Js", CylindricalCoordinateSystem(), location="elem")

apec.write_loads_bcs()
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
apec.apec_settings["R0"] = 6.2
apec.apec_settings["max_picard_iterations"] = 70
apec.apec_settings["min_distance"] = 0.5
apec.apec_settings["tolerance_value"] = 1e-3
apec.apec_settings["alpha_mp"] = 1.0001
apec.apec_settings["alpha_np"] = 0.83197
apec.apec_settings["equilibrium_mesh"] = equilibrium_filename
apec.apec_settings["static_mesh"] = static_filename
apec.apec_settings["evolution_mesh"] = filename
apec.apec_settings["AF0"] = 5.3 * 6.2
# apec.apec_settings["plasma_current"] = 15e6
# apec.apec_settings["beta_p"] = 6.8216e-001
apec.apec_settings["print_full_output"] = 0
apec.apec_settings["toroidal_slices_nodes"] = Cyclic.build_toroidal_slices(
    apec.get_mesh(equilibrium_filename), apec.mesh, "node", n_slices=3, tolerance=1e-3
)
apec.apec_settings["toroidal_slices_elems"] = Cyclic.build_toroidal_slices(
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
apec.time_parameters["theta"] = 1.0
for step_idx in range(len(time_steps)):
    apec.time_parameters["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.time_parameters["transient_effects"] = 0
        apec.time_parameters["time"] = 1e-10
    else:
        apec.time_parameters["transient_effects"] = 1
    apec.apec_settings["plasma_current"] = currents[step_idx]
    apec.apec_settings["beta_p"] = betas[step_idx]
    apec.write_time_step(step_idx + 1)
    if step_idx == 0:
        apec.apec_settings.pop("toroidal_slices_nodes", None)
        apec.apec_settings.pop("toroidal_slices_elems", None)

apec.write_model_data()
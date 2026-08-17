from apec_m import APEC_M
import numpy as np
import vtk
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import pandas as pd
from vtkmodules.util import numpy_support as vtknp
# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()
filename = "ITER_3D_APEC_evol"
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
theta_eps = 1e-3
max_rad = 100.0
extr_angle = 2.0

# # Read the evolution mesh
apec.read_mesh(filepath)
apec.algorithm.equilibrium_mesh = ModelEdit.read_mesh(output_eq_mesh_filepath, reset_arrays=False)

apec.read_excel("input_template.xlsx")

# # Set Boundary conditions
# Axis boundary conditions
Selector.select_pos(
    apec.mesh, 0, (-eps, +eps), coordinate_system=YAxisCylindricalCoordinateSystem()
)
Selector.set_value_on_selection(apec.mesh, "_A", 0.0, 0)
# Selector.set_value_on_selection(apec.mesh, "_A", 0.0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0.0, 2)

# Outer boundary far field zero boundary conditions
Selector.select_pos(
    apec.mesh, 0, (max_rad - eps, max_rad + eps), coordinate_system=SphericalCoordinateSystem()
)
Selector.set_value_on_selection(apec.mesh, "_A", 0.0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0.0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0.0, 2)


# Cyclic periodic boundaries
ModelEdit.remove_array(apec.mesh, "CP_0")
ModelEdit.add_array(apec.mesh, "CP_0", dtype=vtk.VTK_INT)
Selector.select_pos(apec.mesh, 1, (-eps, eps), coordinate_system=YAxisCylindricalCoordinateSystem())
apec.dof_couplings = pd.concat([apec.dof_couplings, ModelEdit.periodic_boundaries(apec.mesh, YAxisCylindricalCoordinateSystem(), dy=extr_angle, array_name="CP_0", variables=("A", "VOLT"))])

# rotate angles
Selector.select_pos(apec.mesh, 1, (extr_angle - eps, extr_angle + eps), coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "angles", extr_angle, 2)

Selector.select(apec.mesh, "ID", 16170)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

# 1. Build source->destination element map using toroidal slices
Selector.select(apec.mesh, "COILS", 1, type="elem", mode="a")
apec.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(apec.mesh, "COORDINATE_SYSTEM", 1, type="elem")
coils = ["PFC_1", "PFC_2", "PFC_3", "PFC_4", "PFC_5", "PFC_6",
         "CS3L", "CS3U", "CS2L", "CS2U", "CS1L", "CS1U"]
for coil in coils:
    Selector.select(apec.algorithm.equilibrium_mesh, coil, 1, type="elem", mode="s")
    eq_cell_data = apec.algorithm.equilibrium_mesh.GetCellData()
    first_elem = vtknp.vtk_to_numpy(
        eq_cell_data.GetArray("selection_cell")
    ).nonzero()[0][0]
    js_z = float(eq_cell_data.GetArray("Js").GetTuple(first_elem)[2])
    Selector.select(apec.mesh, coil, 1, type="elem", mode="s")
    Selector.set_value_on_selection(apec.mesh, "Js", js_z, 1, type="elem")

apec.save_mesh()

apec.algorithm.input_parameters["R0"] = 6.2
apec.algorithm.input_parameters["max_picard_iterations"] = 70
apec.algorithm.input_parameters["tolerance_value"] = 1e-4
apec.algorithm.input_parameters["alpha_mp"] = 1.0001
apec.algorithm.input_parameters["alpha_np"] = 0.83197
apec.algorithm.input_parameters["AF0"] = 5.3 * 6.2
apec.algorithm.input_parameters["halo_delta1"] = 0.04
apec.algorithm.input_parameters["halo_delta2"] = 0
apec.algorithm.input_parameters["t0_halo"] = 1
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
apec.time_step["theta"] = 1.0

apec.load_mesh_into_solver()
for step_idx in range(len(time_steps)):
    apec.time_step["time"] = time_steps[step_idx]
    if step_idx == 0:
        apec.time_step["transient_effects"] = False
        apec.time_step["time"] = 1e-10
    else:
        apec.time_step["transient_effects"] = True
    apec.algorithm.input_parameters["plasma_current"] = currents[step_idx]
    apec.algorithm.input_parameters["beta_p"] = betas[step_idx]
    apec.algorithm.evolution_solve()
    apec.write_output()

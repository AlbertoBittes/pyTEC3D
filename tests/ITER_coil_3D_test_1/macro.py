from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.model_edit import ModelEdit
from pytec3d.selector import Selector
import os
import vtk.util.numpy_support as vtknp
import numpy as np

pytec3d = pyTEC3DSolver()

filename = "ITER_coil_3D_test_1"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_excel("input_template.xlsx")

Selector.select(pytec3d.model.mesh, "X_BOUND", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select(pytec3d.model.mesh, "Z_BOUND", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

Selector.select(pytec3d.model.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select(pytec3d.model.mesh, "ROD_1", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)
Selector.select(pytec3d.model.mesh, "ROD_2", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

Selector.select(pytec3d.model.mesh, "INFINITY_BOUNDARY_INNER", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.model.time_step["theta"] = 1.0

for i in range(1, 10):
    pytec3d.model.time_step["time_step_number"] = i
    if i == 1:
        pytec3d.model.time_step["time"] = i * 0.0001
        pytec3d.model.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.model.time_step["theta"] = 1.0  # autots,off (keep as before)
    else:
        pytec3d.model.time_step["time"] = i * 0.0001
        pytec3d.model.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.model.time_step["theta"] = 1.0  # autots,off (default)

    # select coil elements
    Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
    selection_cells = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetCellData().GetArray("selection_cell")
    )
    selected_ids = np.where(selection_cells == 1)[0]
    centroids_array = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetCellData().GetArray("cell_centroids")
    )
    js_array = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetCellData().GetArray("Js"))
    for l in selected_ids:
        # Get the center coordinates of the element (el_num)
        # This depends on your mesh structure; assuming pytec3d.model.mesh.get_element_center(el_num) returns (x, z)
        x_coord = centroids_array[l, 0]
        z_coord = centroids_array[l, 2]
        module_square = x_coord**2 + z_coord**2
        module = np.sqrt(module_square)
        cos_theta = z_coord / module
        sin_theta = x_coord / module
        fact = i**6
        js_x = cos_theta * fact
        js_z = -sin_theta * fact
        # Apply Js vector to the element (cell_data)
        js_array[l, 0] = js_x
        js_array[l, 2] = js_z
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

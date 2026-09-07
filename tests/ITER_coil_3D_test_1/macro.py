from apec_m_solver import APEC_M
from apec_m.model_edit import ModelEdit
from apec_m.selector import Selector
import os
import vtk.util.numpy_support as vtknp
import numpy as np

apec = APEC_M()

filename = "ITER_coil_3D_test_1"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

Selector.select(apec.model.mesh, "X_BOUND", 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select(apec.model.mesh, "Z_BOUND", 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)

Selector.select(apec.model.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select(apec.model.mesh, "ROD_1", 1)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)
Selector.select(apec.model.mesh, "ROD_2", 1)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

Selector.select(apec.model.mesh, "INFINITY_BOUNDARY_INNER", 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

apec.load_model_into_solver()

apec.model.time_step["theta"] = 1.0

for i in range(1, 10):
    apec.model.time_step["time_step_number"] = i
    if i == 1:
        apec.model.time_step["time"] = i * 0.0001
        apec.model.time_step["transient_effects"] = 0  # TIMINT,OFF
        apec.model.time_step["theta"] = 1.0  # autots,off (keep as before)
    else:
        apec.model.time_step["time"] = i * 0.0001
        apec.model.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.model.time_step["theta"] = 1.0  # autots,off (default)

    # select coil elements
    Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
    selection_cells = vtknp.vtk_to_numpy(
        apec.model.mesh.GetCellData().GetArray("selection_cell")
    )
    selected_ids = np.where(selection_cells == 1)[0]
    centroids_array = vtknp.vtk_to_numpy(
        apec.model.mesh.GetCellData().GetArray("cell_centroids")
    )
    js_array = vtknp.vtk_to_numpy(apec.model.mesh.GetCellData().GetArray("Js"))
    for l in selected_ids:
        # Get the center coordinates of the element (el_num)
        # This depends on your mesh structure; assuming apec.model.mesh.get_element_center(el_num) returns (x, z)
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
    apec.solve()
    apec.write_output()

apec.save_model()

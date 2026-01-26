from apec_m import APEC_M
from selector import Selector
import os
import vtk.util.numpy_support as vtknp
import numpy as np

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_inf_3D"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

apec.write_model_data()

Selector.select(apec.mesh, "X_BOUND", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "Z_BOUND", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)

Selector.select(apec.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "ROD_1", 1)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)
Selector.select(apec.mesh, "ROD_2", 1)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

Selector.select(apec.mesh, "INFINITY_BOUNDARY_INNER", 1)
# Selector.set_value_on_selection(apec.mesh, "boundary_id_0", 111)

apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 0.75
apec.time_parameters["write_substep_results"] = 1

for i in range(1, 10):
    apec.time_parameters["time_step_number"] = i
    if i == 1:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 0  # TIMINT,OFF
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 1  # nsub,1
        apec.time_parameters["theta"] = 0.75  # autots,off (keep as before)
    elif i == 6:
        apec.time_parameters["time"] = 1e-3
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 0  # kbc,1
        apec.time_parameters["number_of_substeps"] = 3  # nsub,3
        apec.time_parameters["theta"] = 1.0  # tintp,1
    elif i == 7:
        apec.time_parameters["time"] = 2e-3
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 2  # nsub,2
        apec.time_parameters["theta"] = 0.8  # tintp,0.8
    elif i == 8:
        apec.time_parameters["time"] = 3e-3
        apec.time_parameters["transient_effects"] = 0  # TIMINT,OFF
        apec.time_parameters["ramped_loads"] = 0  # kbc,1
        apec.time_parameters["number_of_substeps"] = 2  # nsub,2
        apec.time_parameters["theta"] = 0.8  # tintp,0.8
    elif i == 9:
        apec.time_parameters["time"] = 4e-3
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 2  # nsub,2
        apec.time_parameters["theta"] = 0.9  # tintp,0.9
    else:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 4  # nsub,4
        apec.time_parameters["theta"] = 0.75  # autots,off (default)

    # select coil elements
    Selector.select(apec.mesh, "ansys_material_type", 23, type="elem")
    selection_cells = vtknp.vtk_to_numpy(
        apec.mesh.GetCellData().GetArray("selection_cell")
    )
    selected_ids = np.where(selection_cells == 1)[0]
    centroids_array = vtknp.vtk_to_numpy(
        apec.mesh.GetCellData().GetArray("cell_centroids")
    )
    js_array = vtknp.vtk_to_numpy(apec.mesh.GetCellData().GetArray("Js"))
    for l in selected_ids:
        # Get the center coordinates of the element (el_num)
        # This depends on your mesh structure; assuming apec.mesh.get_element_center(el_num) returns (x, z)
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

    apec.write_load_step(i)

apec.save_mesh()

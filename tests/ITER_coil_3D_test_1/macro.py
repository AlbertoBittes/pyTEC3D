from apec_m import APEC_M
from mesh_edit import MeshEdit
import os
import vtk.util.numpy_support as vtknp
import numpy as np

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_3D_test_1"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

MeshEdit.reset_vtk_arrays_to_default((apec.mesh))

apec.read_excel("input_template.xlsx")

apec.write_model_data()

apec.select("X_BOUND", 1)
apec.set_value_on_selection("A", 0, 1)
apec.set_value_on_selection("A", 0, 2)

apec.select("Z_BOUND", 1)
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 1)

apec.select("AXIS_NODES", 1)
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 1)
apec.set_value_on_selection("A", 0, 2)

apec.select("ROD_1", 1)
apec.set_value_on_selection("VOLT", 0)
apec.select("ROD_2", 1)
apec.set_value_on_selection("VOLT", 0)

apec.select("INFINITY_BOUNDARY_INNER", 1)
apec.set_value_on_selection("A", 0, 0)
apec.set_value_on_selection("A", 0, 1)
apec.set_value_on_selection("A", 0, 2)

apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 1.0
apec.time_parameters["write_substep_results"] = 1

for i in range(1, 10):
    apec.time_parameters["time_step_number"] = i
    if i == 1:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 0  # TIMINT,OFF
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 1  # nsub,1
        apec.time_parameters["theta"] = 1.0  # autots,off (keep as before)
    else:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 4  # nsub,4
        apec.time_parameters["theta"] = 1.0  # autots,off (default)

    # select coil elements
    apec.select("ansys_material_type", 23, type="elem")
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

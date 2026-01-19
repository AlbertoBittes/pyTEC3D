from apec_m import APEC_M
from mesh_edit import MeshEdit
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "coil_2D_noint"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")
apec.write_model_data()

apec.select_all()
apec.nsel_ext()
apec.set_value_on_selection("A", 0, 2)

apec.add_array_to_grid("CP_0", dtype=vtk.VTK_INT)
apec.select("ansys_material_type", 23, type="elem")
apec.nsle()
apec.set_value_on_selection("CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(apec.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

apec.select_all()
apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 1
apec.time_parameters["write_substep_results"] = 1
apec.time_parameters["number_of_substeps"] = 4

for i in range(1, 3):
    apec.time_parameters["time_step_number"] = i
    apec.time_parameters["time"] = i * 0.00001
    if i == 1:
        apec.select("ID", master_node)
        apec.set_value_on_selection("AMPS", 100)
    else:
        apec.select("ID", master_node)
        apec.set_value_on_selection("AMPS", 0)
        apec.set_value_on_selection("V", 10000)
    apec.write_load_step(i)

apec.save_mesh()

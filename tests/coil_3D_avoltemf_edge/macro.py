from apec_m import APEC_M
from selector import Selector
from cyclic import Cyclic
from coordinate_system import CylindricalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import numpy as np
# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "coil_3D_avoltemf_edge"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

size = 1
coil = 0.2
extr_angle = 30

Selector.select(apec.mesh, "ansys_material_type", 23, type="elem")
apec.real_constants = Cyclic.coil_winding(apec.mesh, CylindricalCoordinateSystem(),
                                          (coil ** 2, 2000, np.pi*size*coil*coil/(360/extr_angle), 20, 1.0))
apec.write_model_data()

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "Az", 0)

Selector.select(apec.mesh, "ansys_material_type", 23, type="elem")
Selector.set_value_on_selection(apec.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(apec.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]


Selector.select_all(apec.mesh)
apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 0.5
apec.time_parameters["write_substep_results"] = 1
apec.time_parameters["number_of_substeps"] = 4

for i in range(1, 3):
    apec.time_parameters["time_step_number"] = i
    apec.time_parameters["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 100)
    else:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 0)
        Selector.set_value_on_selection(apec.mesh, "VOLT", 10000)
    apec.write_load_step(i)

apec.save_mesh()

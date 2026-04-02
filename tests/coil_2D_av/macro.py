from apec_m import APEC_M
from apec_m.selector import Selector
import os
import vtk.util.numpy_support as vtknp

apec = APEC_M()

filename = "coil_2D_av"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.select_all(apec.mesh)
Selector.nsel_ext(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(apec.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

apec.load_mesh_into_solver()

Selector.select_all(apec.mesh)
apec.time_step["theta"] = 0.5

for i in range(1, 3):
    apec.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 100)
    else:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 0)
        Selector.set_value_on_selection(apec.mesh, "_VOLT", 10000)
    apec.solve()
    apec.write_output()

apec.save_mesh()

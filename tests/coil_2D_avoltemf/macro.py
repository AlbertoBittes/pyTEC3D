from apec_m_solver import APEC_M
from apec_m.selector import Selector
import os
import vtk.util.numpy_support as vtknp

apec = APEC_M()

filename = "coil_2D_avoltemf"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

Selector.select_all(apec.model.mesh)
Selector.nsel_ext(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(apec.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

Selector.select(apec.model.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CP_0", 2)
master_node_ground = vtknp.vtk_to_numpy(apec.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]
Selector.select(apec.model.mesh, "ID", master_node_ground)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

apec.load_model_into_solver()
apec.model.time_step["theta"] = 0.5

Selector.select_all(apec.model.mesh)
for i in range(1, 3):
    apec.model.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.model.mesh, "CURR", 100)
    else:
        Selector.select(apec.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.model.mesh, "CURR", 0)
        Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 10000)
    apec.solve()
    apec.write_output()

apec.save_model()

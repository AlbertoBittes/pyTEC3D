from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.coordinate_system import YAxisCylindricalCoordinateSystem
from apec_m.tree_gauge import TreeGauge
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
import pandas as pd

apec = APEC_M()

filename = "coil_3D_avoltemf_edge_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)
Selector.select_all(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CONSTANT", 1, type="elem")
apec.model.read_excel("input_template.xlsx")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

ModelEdit.remove_array(apec.model.mesh, "CP_1")
ModelEdit.add_array(apec.model.mesh, "CP_1", dtype=vtk.VTK_INT)

Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
apec.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")


Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(apec.model.mesh)
Selector.set_value_on_selection(apec.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    apec.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

Selector.select_pos(
    apec.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
apec.model.dof_couplings = pd.concat(
    [
        apec.model.dof_couplings,
        ModelEdit.periodic_boundaries(
            apec.model.mesh,
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("Az", "V"),
        ),
    ]
)

Selector.select_all(apec.model.mesh)
Selector.select_pos(
    apec.model.mesh,
    0,
    (size-tol, size+tol),
    mode="s",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.model.mesh,
    2,
    (size-tol, size+tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.model.mesh,
    2,
    (-size-tol, -size+tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.model.mesh, "_Az", 0)
TreeGauge.tree_gauge(
    apec.model.mesh,
    periodic_point_coupling="CP_1"
)

Selector.select(apec.model.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.model.mesh)
# Take the first selected node as master node (lowest ID)
master_node_passive = vtknp.vtk_to_numpy(
    apec.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
Selector.select(apec.model.mesh, "ID", master_node_passive)
Selector.set_value_on_selection(apec.model.mesh, "_V", 0.0)

apec.load_model_into_solver()

apec.model.time_step["theta"] = 1.0
apec.save_model()
for i in range(1, 3):
    apec.model.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.model.mesh, "CURR", 100)
        Selector.set_value_on_selection(apec.model.mesh, "_VOLT", np.nan)
    else:
        Selector.select(apec.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.model.mesh, "CURR", 0)
        Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 10000)
    apec.solve()
    apec.write_output()

apec.save_model()

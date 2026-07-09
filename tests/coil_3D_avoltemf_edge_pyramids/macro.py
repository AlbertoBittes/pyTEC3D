from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.cyclic import Cyclic
from apec_m.model_edit import ModelEdit
from apec_m.coordinate_system import CylindricalCoordinateSystem
from apec_m.tree_gauge import TreeGauge
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
import pandas as pd

apec = APEC_M()

filename = "coil_3D_avoltemf_edge_pyramids"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

ModelEdit.remove_array(apec.mesh, "CP_1")
ModelEdit.add_array(apec.mesh, "CP_1", dtype=vtk.VTK_INT)

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
apec.constants = pd.concat([apec.constants, Cyclic.coil_winding(
    apec.mesh,
    CylindricalCoordinateSystem(),
    (coil**2, 2000, np.pi * size * coil * coil / (360 / extr_angle), 20, 1.0),
)], ignore_index=True)


Selector.select_pos(
    apec.mesh,
    2,
    (1-tol, 1+tol),
    mode="s",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    2,
    (-1-tol, -1+tol),
    mode="a",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    0,
    (1-tol, 1+tol),
    mode="a",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    0,
    (-tol, +tol),
    mode="a",
    coordinate_system=CylindricalCoordinateSystem(),
)
# Selector.select_all(apec.mesh)
# Selector.nsel_ext(apec.mesh)
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
Selector.set_value_on_selection(apec.mesh, "_Az", 0)

Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    apec.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh)
# Take the first selected node as master node (lowest ID)
master_node_passive = vtknp.vtk_to_numpy(
    apec.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
Selector.select(apec.mesh, "ID", master_node_passive)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0.0)

Selector.select_pos(
    apec.mesh, 1, (-tol, tol), coordinate_system=CylindricalCoordinateSystem()
)
Selector.select_pos(
    apec.mesh, 0, (-tol, tol), mode="u", coordinate_system=CylindricalCoordinateSystem()
)
Selector.select_pos(
    apec.mesh,
    0,
    (size - tol, size + tol),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    2,
    (size - tol, size + tol),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    2,
    (-size - tol, -size + tol),
    mode="u",
    coordinate_system=CylindricalCoordinateSystem(),
)
# Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
# Selector.nsle(apec.mesh, mode="u")
Selector.select(apec.mesh, "MidsideNodes", 1, mode="r")
apec.dof_couplings = pd.concat(
    [
        apec.dof_couplings,
        Cyclic.cpcyc(
            apec.mesh,
            CylindricalCoordinateSystem(),
            extr_angle,
            array_name="CP_1",
            variables=("Az"),
        ),
    ]
)

Selector.select_all(apec.mesh)
# Selector.select_pos(
#     apec.mesh, 1, (extr_angle-tol, extr_angle+tol), coordinate_system=CylindricalCoordinateSystem(), mode="u"
# )
# Selector.select_pos(
#     apec.mesh,
#     0,
#     (size - tol, size + tol),
#     mode="u",
#     coordinate_system=CylindricalCoordinateSystem(),
# )
# Selector.select_pos(
#     apec.mesh,
#     2,
#     (size - tol, size + tol),
#     mode="u",
#     coordinate_system=CylindricalCoordinateSystem(),
# )
# Selector.select_pos(
#     apec.mesh,
#     2,
#     (-size - tol, -size + tol),
#     mode="u",
#     coordinate_system=CylindricalCoordinateSystem(),
# )
Selector.nsel_ext(apec.mesh, mode="u")
TreeGauge.tree_gauge(apec.mesh)

apec.load_mesh_into_solver()

Selector.select_all(apec.mesh)
apec.time_step["theta"] = 1.0
apec.save_mesh()
for i in range(1, 3):
    apec.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 100)
        Selector.set_value_on_selection(apec.mesh, "_VOLT", np.nan)
    else:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 0)
        Selector.set_value_on_selection(apec.mesh, "_VOLT", 10000)
    apec.solve()
    apec.write_output()

apec.save_mesh()

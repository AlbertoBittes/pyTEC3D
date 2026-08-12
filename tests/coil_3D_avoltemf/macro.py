from apec_m import APEC_M
from apec_m.coordinate_system import YAxisCylindricalCoordinateSystem
from apec_m.model_edit import ModelEdit
from apec_m.selector import Selector
import os
import vtk.util.numpy_support as vtknp
import vtk
import numpy as np
import pandas as pd

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "coil_3D_avoltemf"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

ModelEdit.remove_array(apec.mesh, "CP_1")
ModelEdit.add_array(apec.mesh, "CP_1", dtype=vtk.VTK_INT)


Selector.select_pos(
    apec.mesh,
    2,
    (size - tol, size + tol),
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    2,
    (-size - tol, -size + tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    0,
    (size - tol, size + tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)


Selector.select_pos(
    apec.mesh, 0, (-tol, tol), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem()
)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)


# Coupling for stranded coils
Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
apec.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(apec.mesh, "COORDINATE_SYSTEM", 1, type="elem")
Selector.nsle(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    apec.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]


Selector.select_pos(
    apec.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh, mode="u")
apec.dof_couplings = pd.concat(
    [
        apec.dof_couplings,
        ModelEdit.periodic_boundaries(
            apec.mesh,
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A"),
        ),
    ]
)


Selector.select_pos(
    apec.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh, mode="r")
apec.dof_couplings = pd.concat(
    [
        apec.dof_couplings,
        ModelEdit.periodic_boundaries(
            apec.mesh,
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A", "VOLT"),
        ),
    ]
)

Selector.select_pos(
    apec.mesh,
    1,
    (extr_angle - tol, extr_angle + tol),
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "angles", extr_angle, 2)

Selector.select(apec.mesh, "ID", 2220)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

apec.load_mesh_into_solver()

Selector.select_all(apec.mesh)
apec.time_step["theta"] = 1.0

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

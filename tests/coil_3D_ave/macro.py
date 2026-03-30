from apec_m import APEC_M
from cyclic import Cyclic
from coordinate_system import CylindricalCoordinateSystem
from model_edit import ModelEdit
from selector import Selector
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
filename = "coil_3D_ave"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

ModelEdit.remove_array(apec.mesh, "CP_1")
ModelEdit.add_array(apec.mesh, "CP_1", dtype=vtk.VTK_INT)

# Coupling for stranded coils
Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(apec.mesh)
Selector.set_value_on_selection(apec.mesh, "CP_0", 100000)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    apec.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
apec.constants = pd.concat(
    [
        apec.constants,
        Cyclic.coil_winding(
            apec.mesh,
            CylindricalCoordinateSystem(),
            (coil**2, 2000, np.pi * size * coil * coil, 10, (360 / extr_angle)),
        ),
    ]
)

Selector.select_pos(
    apec.mesh,
    2,
    (size - tol, size + tol),
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    2,
    (-size - tol, -size + tol),
    mode="a",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh,
    0,
    (size - tol, size + tol),
    mode="a",
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select_pos(
    apec.mesh, 0, (-tol, tol), mode="s", coordinate_system=CylindricalCoordinateSystem()
)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

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
Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh, mode="u")
apec.dof_couplings = pd.concat(
    [
        apec.dof_couplings,
        Cyclic.cpcyc(
            apec.mesh,
            CylindricalCoordinateSystem(),
            extr_angle,
            array_name="CP_1",
            variables=("A"),
        ),
    ]
)

Selector.select_pos(
    apec.mesh,
    1,
    (extr_angle - tol, extr_angle + tol),
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh, 0, (-tol, tol), mode="u", coordinate_system=CylindricalCoordinateSystem()
)
Selector.select_pos(
    apec.mesh,
    0,
    (1 - tol, 1 + tol),
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
Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh, mode="u")
Selector.set_value_on_selection(apec.mesh, "angles", extr_angle, 2)

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
Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh, mode="r")
apec.dof_couplings = pd.concat(
    [
        apec.dof_couplings,
        Cyclic.cpcyc(
            apec.mesh,
            CylindricalCoordinateSystem(),
            extr_angle,
            array_name="CP_1",
            variables=("A", "V"),
        ),
    ]
)

Selector.select_pos(
    apec.mesh,
    1,
    (extr_angle - tol, extr_angle + tol),
    coordinate_system=CylindricalCoordinateSystem(),
)
Selector.select_pos(
    apec.mesh, 0, (-tol, tol), mode="u", coordinate_system=CylindricalCoordinateSystem()
)
Selector.select_pos(
    apec.mesh,
    0,
    (1 - tol, 1 + tol),
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
Selector.select(apec.mesh, "MATERIAL", 24, type="elem")
Selector.nsle(apec.mesh, mode="r")
Selector.set_value_on_selection(apec.mesh, "angles", extr_angle, 2)

Selector.select(apec.mesh, "ID", 2220)
Selector.set_value_on_selection(apec.mesh, "V", 0)

apec.write_model_data()

Selector.select_all(apec.mesh)
apec.time_parameters["theta"] = 1.0

for i in range(1, 3):
    apec.time_parameters["time_step_number"] = i
    apec.time_parameters["time"] = i * 0.00001
    if i == 1:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 100)
    else:
        Selector.select(apec.mesh, "ID", master_node)
        Selector.set_value_on_selection(apec.mesh, "AMPS", 0)
        Selector.set_value_on_selection(apec.mesh, "V", 10000)
    apec.write_load_step(i)

apec.save_mesh()

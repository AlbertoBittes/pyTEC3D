from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
from pytec3d.model_edit import ModelEdit
from pytec3d.selector import Selector
import os
import vtk.util.numpy_support as vtknp
import vtk
import pandas as pd

pytec3d = pyTEC3DSolver()

filename = "coil_3D_avoltemf_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_excel("input_template.xlsx")

size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

ModelEdit.remove_array(pytec3d.model.mesh, "CP_1")
ModelEdit.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

Selector.select_pos(
    pytec3d.model.mesh,
    2,
    (size - tol, size + tol),
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.select_pos(
    pytec3d.model.mesh,
    2,
    (-size - tol, -size + tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.select_pos(
    pytec3d.model.mesh,
    0,
    (size - tol, size + tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)


Selector.select_pos(
    pytec3d.model.mesh, 0, (-tol, tol), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem()
)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)


# Coupling for stranded coils
Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
pytec3d.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]


Selector.select_pos(
    pytec3d.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
Selector.select(pytec3d.model.mesh, "MATERIAL", 24, type="elem")
pytec3d.model.dof_couplings = pd.concat(
    [
        pytec3d.model.dof_couplings,
        ModelEdit.periodic_boundaries(
            pytec3d.model.mesh,
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A", "V"),
        ),
    ]
)

Selector.select_pos(
    pytec3d.model.mesh,
    1,
    (extr_angle - tol, extr_angle + tol),
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
Selector.set_value_on_selection(pytec3d.model.mesh, "angles", extr_angle, 2)

Selector.select(pytec3d.model.mesh, "ID", 2220)
Selector.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

pytec3d.load_model_into_solver()

Selector.select_all(pytec3d.model.mesh)
pytec3d.model.time_step["theta"] = 1.0

for i in range(1, 3):
    pytec3d.model.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(pytec3d.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        Selector.select(pytec3d.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

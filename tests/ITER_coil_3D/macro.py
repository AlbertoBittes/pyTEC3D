from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
from pytec3d.selector import Selector
from pytec3d.model_edit import ModelEdit
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

pytec3d = pyTEC3DSolver()

filename = "ITER_coil_3D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
Selector.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")

Selector.select_pos(pytec3d.model.mesh, 1, (-0.0001, 0.0001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select(pytec3d.model.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(pytec3d.model.mesh, "angles", 90.0, 2)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select_pos(pytec3d.model.mesh, 1, (89.999, 90.001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select(pytec3d.model.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(pytec3d.model.mesh, "angles", 90.0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "angles", 90.0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

Selector.select(pytec3d.model.mesh, "ROD_2", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

Selector.select(pytec3d.model.mesh, "ROD_1", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 0)

Selector.select(pytec3d.model.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="s", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.75

for i in range(1, 10):
    if i == 1:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 0.75  # autots,off (keep as before)
    elif i == 6:
        pytec3d.time_step["time"] = 1e-3
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 1.0  # tintp,1
    elif i == 7:
        pytec3d.time_step["time"] = 2e-3
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 8:
        pytec3d.time_step["time"] = 3e-3
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 9:
        pytec3d.time_step["time"] = 4e-3
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 0.9  # tintp,0.9
    else:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 0.75  # autots,off (default)

    # select coil elements
    Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
    fact = i**6
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", fact, 1, type="elem")

    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

from apec_m_solver import APEC_M
from apec_m_coordinate_system.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

apec = APEC_M()

filename = "ITER_coil_3D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
apec.model.coordinate_systems[1] = YAxisCylindricalCoordinateSystem()
Selector.set_value_on_selection(apec.model.mesh, "COORDINATE_SYSTEM", 1, type="elem")

Selector.select_pos(apec.model.mesh, 1, (-0.0001, 0.0001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select(apec.model.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(apec.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.model.mesh, "angles", 90.0, 2)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select_pos(apec.model.mesh, 1, (89.999, 90.001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select(apec.model.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(apec.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.model.mesh, "angles", 90.0, 0)
Selector.set_value_on_selection(apec.model.mesh, "angles", 90.0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)

Selector.select(apec.model.mesh, "ROD_2", 1)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

Selector.select(apec.model.mesh, "ROD_1", 1)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

Selector.select(apec.model.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select_pos(apec.model.mesh, 0, (0.999, 1.001), mode="s", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

apec.load_model_into_solver()

apec.model.time_step["theta"] = 0.75

for i in range(1, 10):
    if i == 1:
        apec.model.time_step["time"] = i * 0.0001
        apec.model.time_step["transient_effects"] = 0  # TIMINT,OFF
        apec.model.time_step["theta"] = 0.75  # autots,off (keep as before)
    elif i == 6:
        apec.model.time_step["time"] = 1e-3
        apec.model.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.model.time_step["theta"] = 1.0  # tintp,1
    elif i == 7:
        apec.model.time_step["time"] = 2e-3
        apec.model.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.model.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 8:
        apec.model.time_step["time"] = 3e-3
        apec.model.time_step["transient_effects"] = 0  # TIMINT,OFF
        apec.model.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 9:
        apec.model.time_step["time"] = 4e-3
        apec.model.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.model.time_step["theta"] = 0.9  # tintp,0.9
    else:
        apec.model.time_step["time"] = i * 0.0001
        apec.model.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.model.time_step["theta"] = 0.75  # autots,off (default)

    # select coil elements
    Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
    fact = i**6
    Selector.set_value_on_selection(apec.model.mesh, "Js", fact, 1, type="elem")

    apec.solve()
    apec.write_output()

apec.save_model()

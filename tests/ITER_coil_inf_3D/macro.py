from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import numpy as np

apec = APEC_M()

filename = "ITER_coil_inf_3D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.select_pos(apec.mesh, 1, (-0.0001, 0.0001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select(apec.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

Selector.select_pos(apec.mesh, 1, (-90.0001, -89.999), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select(apec.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 0)
Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 1)
Selector.set_value_on_selection(apec.mesh, "angles", -90.0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 1)

Selector.select(apec.mesh, "ROD_2", 1)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

Selector.select(apec.mesh, "ROD_1", 1)
Selector.set_value_on_selection(apec.mesh, "_VOLT", 0)

Selector.select(apec.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="s", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "_A", 0, 2)

apec.load_mesh_into_solver()

for i in range(1, 10):
    apec.time_step["time_step_number"] = i
    if i == 1:
        apec.time_step["time"] = i * 0.0001
        apec.time_step["transient_effects"] = 0  # TIMINT,OFF
        apec.time_step["theta"] = 0.75  # autots,off (keep as before)
    elif i == 6:
        apec.time_step["time"] = 1e-3
        apec.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.time_step["theta"] = 1.0  # tintp,1
    elif i == 7:
        apec.time_step["time"] = 2e-3
        apec.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 8:
        apec.time_step["time"] = 3e-3
        apec.time_step["transient_effects"] = 0  # TIMINT,OFF
        apec.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 9:
        apec.time_step["time"] = 4e-3
        apec.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.time_step["theta"] = 0.9  # tintp,0.9
    else:
        apec.time_step["time"] = i * 0.0001
        apec.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.time_step["theta"] = 0.75  # autots,off (default)

        # select coil elements
    Selector.select(apec.mesh, "MATERIAL", 23, type="elem")
    selection_cells = vtknp.vtk_to_numpy(
        apec.mesh.GetCellData().GetArray("selection_cell")
    )
    Selector.set_value_on_selection(apec.mesh, "Js", 0, 0, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", 0, 1, type="elem")
    Selector.set_value_on_selection(apec.mesh, "Js", 0, 2, type="elem")
    fact = i**6
    Selector.set_value_on_selection(apec.mesh, "Js", fact, 1, type="elem")
    ModelEdit.map_vector_field_to_global_components(apec.mesh, "Js",
                                                   coordinate_system=YAxisCylindricalCoordinateSystem(),
                                                   location="elem")
    apec.solve()
    apec.write_output()

apec.save_mesh()

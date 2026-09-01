from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.coordinate_system import SphericalCoordinateSystem, YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

apec = APEC_M()

filename = "ITER_coil_3D_period"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_model(filepath)

apec.model.read_excel("input_template.xlsx")

Selector.select_pos(apec.model.mesh, 1, (-90.0001, -89.999), coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select_pos(apec.model.mesh, 0, (-0.0001, 0.0001), mode="u", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select_pos(apec.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.model.mesh, "angles", -90.0, 2)

ModelEdit.remove_array(apec.model.mesh, "CP_0")
ModelEdit.add_array(apec.model.mesh, "CP_0", dtype=vtk.VTK_INT)
Selector.select_pos(apec.model.mesh, 1, 0.0, coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select_pos(apec.model.mesh,  0, (-0.0001, 0.0001), mode="u", coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.select_pos(apec.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
apec.model.dof_couplings = ModelEdit.periodic_boundaries(apec.model.mesh, YAxisCylindricalCoordinateSystem(), dy=-90.0, variables=("A", "VOLT"))

Selector.select(apec.model.mesh, "TYPE", 24, type="elem")
Selector.nsle(apec.model.mesh, "r")

Selector.select_pos(apec.model.mesh, 0, (-0.0001, 0.0001), coordinate_system=YAxisCylindricalCoordinateSystem())
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

Selector.select(apec.model.mesh, "ID", 107)
Selector.set_value_on_selection(apec.model.mesh, "_VOLT", 0)

Selector.select_pos(apec.model.mesh, 0, (0.999, 1.001), coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(apec.model.mesh, "_A", 0, 2)

apec.load_model_into_solver()

apec.model.time_step["theta"] = 0.75

for i in range(1, 3):
    if i == 1:
        apec.model.time_step["time"] = i * 0.0001
        apec.model.time_step["transient_effects"] = 0  # TIMINT,OFF
        apec.model.time_step["theta"] = 1.0  # autots,off (keep as before)
    else:
        apec.model.time_step["time"] = i * 0.0001
        apec.model.time_step["transient_effects"] = 1  # TIMINT,ON
        apec.model.time_step["theta"] = 1.0  # autots,off (default)

    # select coil elements
    Selector.select(apec.model.mesh, "MATERIAL", 23, type="elem")
    selection_cells = vtknp.vtk_to_numpy(
        apec.model.mesh.GetCellData().GetArray("selection_cell")
    )
    Selector.set_value_on_selection(apec.model.mesh, "Js", 0, 0, type="elem")
    Selector.set_value_on_selection(apec.model.mesh, "Js", 0, 1, type="elem")
    Selector.set_value_on_selection(apec.model.mesh, "Js", 0, 2, type="elem")
    fact = i**6
    Selector.set_value_on_selection(apec.model.mesh, "Js", fact, 1, type="elem")
    ModelEdit.map_vector_field_to_global_components(apec.model.mesh, "Js",
                                                   coordinate_system=YAxisCylindricalCoordinateSystem(),
                                                   location="elem")
    apec.solve()
    apec.write_output()

apec.save_model()

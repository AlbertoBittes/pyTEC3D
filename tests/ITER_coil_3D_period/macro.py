from pytec3d import mesh, model_edit, selection, Solver
from pytec3d.coordinate_system import SphericalCoordinateSystem, YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

pytec3d = Solver()

filename = "ITER_coil_3D_period"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

selection.select_pos(pytec3d.model.mesh, 1, (-90.0001, -89.999), coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 0, (-0.0001, 0.0001), mode="u", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", -90.0, 2)

mesh.remove_array(pytec3d.model.mesh, "CP_0")
mesh.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
selection.select_pos(pytec3d.model.mesh, 1, 0.0, coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh,  0, (-0.0001, 0.0001), mode="u", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(), dy=-90.0, variables=("A", "V"))

selection.select(pytec3d.model.mesh, "TYPE", 24, location="elem")
selection.nsle(pytec3d.model.mesh, "r")

selection.select_pos(pytec3d.model.mesh, 0, (-0.0001, 0.0001), coordinate_system=YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select(pytec3d.model.mesh, "ID", 107)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.75

for i in range(1, 3):
    if i == 1:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 1.0  # autots,off (keep as before)
    else:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 1.0  # autots,off (default)

    # select coil elements
    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    selection_cells = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetCellData().GetArray("selection_cell")
    )
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 0, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 1, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 2, location="elem")
    fact = i**6
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", fact, 1, location="elem")
    model_edit.map_vector_field_to_global_components(pytec3d.model.mesh, "Js",
                                                   coordinate_system=YAxisCylindricalCoordinateSystem(),
                                                   location="elem")
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

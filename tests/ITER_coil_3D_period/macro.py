from apec_m import APEC_M
from apec_m.selector import Selector
from apec_m.model_edit import ModelEdit
from apec_m.coordinate_system import SphericalCoordinateSystem, CylindricalCoordinateSystem
from apec_m.cyclic import Cyclic
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_3D_period"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

# Selector.select_pos(apec.mesh, 1, 0.0, coordinate_system=CylindricalCoordinateSystem())
# Selector.select_pos(apec.mesh, 0, 0, mode="u", coordinate_system=CylindricalCoordinateSystem())
# Selector.select_pos(apec.mesh, 0, 1, mode="u", coordinate_system=SphericalCoordinateSystem())
# Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 2)
# Selector.set_value_on_selection(apec.mesh, "A", 0, 0)

Selector.select_pos(apec.mesh, 1, (-90.0001, -89.999), coordinate_system=CylindricalCoordinateSystem())
Selector.select_pos(apec.mesh, 0, (-0.0001, 0.0001), mode="u", coordinate_system=CylindricalCoordinateSystem())
Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "angles", -90.0, 2)

Selector.select_pos(apec.mesh, 1, 0.0, coordinate_system=CylindricalCoordinateSystem())
Selector.select_pos(apec.mesh,  0, (-0.0001, 0.0001), mode="u", coordinate_system=CylindricalCoordinateSystem())
Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
apec.dof_couplings = Cyclic.cpcyc(apec.mesh, CylindricalCoordinateSystem(), -90.0, variables=("A", "VOLT"))

Selector.select(apec.mesh, "TYPE", 24, type="elem")
Selector.nsle(apec.mesh, "r")

Selector.select_pos(apec.mesh, 0, (-0.0001, 0.0001), coordinate_system=CylindricalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select(apec.mesh, "ID", 107)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)
# Selector.select(apec.mesh, "ROD_2", 1)
# Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

Selector.select_pos(apec.mesh, 0, (0.999, 1.001), coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

apec.write_model_data()

apec.time_parameters["theta"] = 0.75

for i in range(1, 3):
    apec.time_parameters["time_step_number"] = i
    if i == 1:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 0  # TIMINT,OFF
        apec.time_parameters["theta"] = 1.0  # autots,off (keep as before)
    else:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["theta"] = 1.0  # autots,off (default)

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
                                                   coordinate_system=CylindricalCoordinateSystem(),
                                                   location="elem")
    apec.write_load_step(i)

apec.save_mesh()

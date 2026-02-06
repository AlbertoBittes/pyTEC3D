from apec_m import APEC_M
from coordinate_system import CylindricalCoordinateSystem, SphericalCoordinateSystem
from selector import Selector
from model_edit import ModelEdit
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# from convert_mesh import convert_mesh_from_cdb

apec = APEC_M()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_3D"

filepath = os.path.join(os.getcwd(), "tests", filename, filename + ".vtu")

apec.read_mesh(filepath)

apec.read_excel("input_template.xlsx")

Selector.select_pos(apec.mesh, 1, (-0.0001, 0.0001), mode="s", coordinate_system=CylindricalCoordinateSystem())
Selector.select(apec.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 2)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

Selector.select_pos(apec.mesh, 1, (89.999, 90.0001), mode="s", coordinate_system=CylindricalCoordinateSystem())
Selector.select(apec.mesh, "AXIS_NODES", 1, mode="u")
Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 0)
Selector.set_value_on_selection(apec.mesh, "angles", 90.0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)

Selector.select(apec.mesh, "ROD_2", 1)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

Selector.select(apec.mesh, "ROD_1", 1)
Selector.set_value_on_selection(apec.mesh, "VOLT", 0)

Selector.select(apec.mesh, "AXIS_NODES", 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
# Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)


apec.write_model_data()

Selector.select_pos(apec.mesh, 0, (0.999, 1.001), mode="s", coordinate_system=SphericalCoordinateSystem())
Selector.set_value_on_selection(apec.mesh, "A", 0, 0)
Selector.set_value_on_selection(apec.mesh, "A", 0, 1)
Selector.set_value_on_selection(apec.mesh, "A", 0, 2)

apec.time_parameters["ramped_loads"] = 1
apec.time_parameters["theta"] = 0.75
apec.time_parameters["write_substep_results"] = 1

for i in range(1, 10):
    apec.time_parameters["time_step_number"] = i
    if i == 1:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 0  # TIMINT,OFF
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 1  # nsub,1
        apec.time_parameters["theta"] = 0.75  # autots,off (keep as before)
    elif i == 6:
        apec.time_parameters["time"] = 1e-3
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 0  # kbc,1
        apec.time_parameters["number_of_substeps"] = 3  # nsub,3
        apec.time_parameters["theta"] = 1.0  # tintp,1
    elif i == 7:
        apec.time_parameters["time"] = 2e-3
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 2  # nsub,2
        apec.time_parameters["theta"] = 0.8  # tintp,0.8
    elif i == 8:
        apec.time_parameters["time"] = 3e-3
        apec.time_parameters["transient_effects"] = 0  # TIMINT,OFF
        apec.time_parameters["ramped_loads"] = 0  # kbc,1
        apec.time_parameters["number_of_substeps"] = 2  # nsub,2
        apec.time_parameters["theta"] = 0.8  # tintp,0.8
    elif i == 9:
        apec.time_parameters["time"] = 4e-3
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 2  # nsub,2
        apec.time_parameters["theta"] = 0.9  # tintp,0.9
    else:
        apec.time_parameters["time"] = i * 0.0001
        apec.time_parameters["transient_effects"] = 1  # TIMINT,ON
        apec.time_parameters["ramped_loads"] = 1  # kbc,0
        apec.time_parameters["number_of_substeps"] = 4  # nsub,4
        apec.time_parameters["theta"] = 0.75  # autots,off (default)

    # select coil elements
    Selector.select(apec.mesh, "ansys_material_type", 23, type="elem")
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

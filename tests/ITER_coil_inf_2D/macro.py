from pytec3d import selection, Solver
import os

# from convert_mesh import convert_mesh_from_cdb

pytec3d = Solver()

# convert_mesh_from_cdb(
#     r"D:\PlasmaArc\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_inf_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh, "s")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.5

for i in range(1, 7):
    pytec3d.time_step["time_step_number"] = i
    if i == 1:
        pytec3d.time_step["transient_effects"] = 0
        pytec3d.time_step["time"] = i * 0.0001
    elif i == 5:
        pytec3d.time_step["transient_effects"] = 1
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["theta"] = 0.85
    elif i == 6:
        pytec3d.time_step["transient_effects"] = 0
        pytec3d.time_step["time"] = i * 0.001
    else:
        pytec3d.time_step["transient_effects"] = 1
        pytec3d.time_step["time"] = i * 0.0001

    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1 / i, 2, location="elem")
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

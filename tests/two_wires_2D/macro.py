from pytec3d import selection, Solver
import os

pytec3d = Solver()

filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.91

for i in range(1, 3):
    if i == 1:
        pytec3d.time_step["transient_effects"] = 0
        pytec3d.time_step["time"] = i * 0.00001
    else:
        pytec3d.time_step["transient_effects"] = 1
        pytec3d.time_step["time"] = i * 0.00002

    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", i * 1000, 2, location="elem")
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()

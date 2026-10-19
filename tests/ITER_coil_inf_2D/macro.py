# =====================================================================================================================
# ITER_coil_inf_2D - 2D axisymmetric: a source ring with imposed Js and a passive shorted ring, many step types
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A minimal axisymmetric model (one cell per ring) used to check the time integration: a sequence of static and
# transient steps with different theta values, the source current density changing at every step. It shows how the
# time_step dictionary is changed between solves: "time" (end time of the step), "transient_effects" (False = static
# step: no damping term, eddy currents switched off) and "theta". The solver reuses the factorized matrix as long as
# dt, theta, transient_effects and the set of constrained DOFs do not change.
#
# Geometry and mesh (ITER_coil_inf_2D.vtu): axisymmetric about the y axis (x = r), 0 <= r <= 1 m, -1 <= z <= 1 m,
#   330 TRI3 + 2 QUAD4 cells. Source ring TYPE / MATERIAL 23: one cell at r = 0.25 m, z = +0.25 m; passive ring
#   TYPE / MATERIAL 24: one cell at r = 0.25 m, z = -0.25 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV2D_nodal / NON_CONDUCTING            / 1 / 1   air
#   * 23: AV2D_nodal / STRANDED_CONDUCTOR        / 1 / 1   source ring: Js imposed
#   * 24: AV2D_nodal / SOLID_CONDUCTOR           / 1 / 1   passive shorted ring: eddy currents
#   Materials: 21, 23 mur = 1; 24 mur = 1, rho = 1e-9 Ohm m.  Real constants 1: thickness not
#   given: 1 (the full 360 degree ring, the default).
#
# Boundary conditions:  A_z = 0 on the external boundary (axis and far boundary).
# Loads and time history (Js_z = 1/i A/m^2 at step i):
#   * step 1, t = 1e-4 s: static;   steps 2-4, t = 2e-4 ... 4e-4 s: transient, theta = 0.5;
#   * step 5, t = 5e-4 s: transient, theta = 0.85;   step 6, t = 6e-3 s: static.
#   time_step["time_step_number"] is not a time-step setting: the solver ignores it.
# Results: output/ITER_coil_inf_2D_out_<n>.vtu, ITER_coil_inf_2D_steps.csv, ITER_coil_inf_2D.log
#
# How to run:  cd tests/ITER_coil_inf_2D && python macro.py
# Documentation: "User manual - time steps", "Theory - theta method".
# =====================================================================================================================


from pytec3d import selection, Solver
import os

# from convert_mesh import convert_mesh_from_cdb

# 1. Solver and model files. (A .cdb mesh is converted once with the command line tool pytec3d-convert-cdb.)
pytec3d = Solver()

# convert_mesh_from_cdb(
#     r"D:\pyTEC3D\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
filename = "ITER_coil_inf_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Boundary condition: A_z = 0 on the external boundary (mode "s" given positionally).
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh, "s")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Load the model into the solver session; default theta 0.5.
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.5

# 4. Time loop: static / transient steps with changing theta (see the header).
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

    # Source current density of step i: Js_z = 1/i on the source ring.
    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1 / i, 2, location="elem")
    pytec3d.solve()
    pytec3d.write_output()

# 5. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

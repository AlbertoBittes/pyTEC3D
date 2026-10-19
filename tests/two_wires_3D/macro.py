# =====================================================================================================================
# two_wires_3D - 3D: source wire with imposed Js and passive conducting wire, AV3D_nodal
# =====================================================================================================================
#
# What this example shows
# -----------------------
# Two long parallel wires along z in a box: a STRANDED CONDUCTOR with an imposed source current density (Js, element
# load) and a passive massive conductor where eddy currents are induced. One transient step from zero fields.
# Nodal vector potential, first-order Lagrange elements. It shows Dirichlet conditions set by planes (select_pos),
# the MidsideNodes point array to restrict a selection to the midside nodes of a quadratic mesh, and the
# element load Js.
# First-order mesh: no node is a midside node, so the selection of step 4 is empty and the end faces are free.
#
# Geometry and mesh (two_wires_3D.vtu): -1 <= x <= 1 m, -0.5 <= y <= 0.5 m, -0.5 <= z <= 0.5 m;
#   7701 nodes, 33571 TET4 + 1280 HEX8 + 1280 PYRAMID5 cells.
#   * Source wire (TYPE / MATERIAL 23): square cross-section about 0.1 m wide at x = -0.2 m, along z.
#   * Passive wire (TYPE / MATERIAL 24): the same at x = +0.2 m.   * Air (TYPE / MATERIAL 21).
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 1 = AV3D_nodal / NON_CONDUCTING            / 1 / 0   air
#   * 23: 1 = AV3D_nodal / STRANDED_CONDUCTOR        / 1 / 0   source wire: Js imposed
#   * 24: 1 = AV3D_nodal / SOLID_CONDUCTOR           / 1 / 0   passive wire: A + V
#   Materials: 21 mur = diag(1.3, 1, 1); 23 mur = 1; 24 mur = 1, rho = 1e-9 Ohm m.  Real constants 1: none needed
#   (empty row).
#
# Loads and time history: Js_z = 1 A/m^2 in the source wire; one step t = 1e-3 s, theta = 0.5, transient.
# Results: output/two_wires_3D_out_1.vtu (A, V, B, H, Jt), two_wires_3D_steps.csv, two_wires_3D.log
#
# How to run:  cd tests/two_wires_3D && python macro.py
# Documentation: "Element library" (AV3D_nodal), "User manual - boundary conditions", "Examples" (two wires).
# =====================================================================================================================


from pytec3d import selection, Solver
import os

# 1. Solver and model files.
pytec3d = Solver()

filename = "two_wires_3D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Planes x = +-0.5 m (axis 0): all the components of A fixed to zero.
#    Note: the x faces of this mesh are at x = +-1 m; the planes x = +-0.5 m hold (almost) no node, so this
#    selection constrains (almost) nothing.
selection.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Faces y = +-0.5 m (axis 1): A = 0.
selection.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 4. End faces z = +-0.5 m, crossed by the wires: A_x = A_y = 0 (B tangent to the faces, as for
#    infinitely long wires), restricted to the MIDSIDE nodes (MidsideNodes = 1, mode "r").
selection.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

# 5. Zero potential V on the nodes of the passive wire (an equipotential, grounded conductor).
selection.select(pytec3d.model.mesh, "TYPE", 24, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 6. Load: source current density Js_z = 1 A/m^2 (component 2) on the cells of the source wire.
selection.select(pytec3d.model.mesh, "TYPE", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1, 2, location="elem")

# 7. One transient step from zero fields: t = 1e-3 s, Crank-Nicolson (theta = 0.5). save_model() keeps the
#    prepared model before the solve.
pytec3d.time_step["time"] = 1e-3
pytec3d.time_step["theta"] = 0.5
pytec3d.time_step["transient_effects"] = 1
pytec3d.save_model()
# 8. Load the model, solve the step and write output/<mesh_name>_out_1.vtu.
pytec3d.load_model_into_solver()

pytec3d.solve()
pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

# =====================================================================================================================
# two_wires_3D_edge_el_2order - 3D edge elements: source rod with imposed Js and solid rod, second-order V
# =====================================================================================================================
#
# What this example shows
# -----------------------
# Two parallel rods along z crossing a box: a STRANDED CONDUCTOR with an imposed source current density Js and a SOLID
# conductor (non-integrated scalar potential V) in which eddy currents are induced. Edge (Nedelec) elements for A on
# a quadratic hexahedral mesh, Lagrange elements of order 2 for V. It shows:
#   * tangential A = 0 on the six faces of the box, written on the midside nodes (component 2 of _A);
#   * named point arrays of the mesh (LOWER_FACE_ROD, UPPER_FACE_ROD) used to select boundary nodes;
#   * the tree-cotree gauge: all the cells selected, grown from the external boundary nodes;
#   * a static first step followed by transient steps with different theta values and loads.
#
# Geometry and mesh (two_wires_3D_edge_el_2order.vtu): box -0.5 <= x, y, z <= 0.5 m, 1210 HEX20 cells.
#   * Source rod (TYPE / MATERIAL 23): one column of 10 cells along z at x = -0.185 m.
#   * Solid rod (TYPE / MATERIAL 24): one column of 10 cells along z at x = +0.185 m; its end faces are the point
#     arrays LOWER_FACE_ROD (z = -0.5 m) and UPPER_FACE_ROD (z = +0.5 m).
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 2 = AV3D_edge / NON_CONDUCTING     / 2 / 0   air
#   * 23: 2 = AV3D_edge / STRANDED_CONDUCTOR / 2 / 0   source rod: Js imposed
#   * 24: 2 = AV3D_edge / SOLID_CONDUCTOR    / 2 / 0   solid rod: A + V
#   Materials: 21 mur = diag(1.1, 1.2, 1.3); 23 mur = 1 (rho unused); 24 mur = 1, rho = 1e-9 Ohm m.
#   Real constants 1: none needed (empty row).
#
# Loads and time history (three steps; the factorization is rebuilt when dt or theta changes):
#   * step 1, t = 1e-3 s, STATIC (transient_effects = 0): Js_z = 1 A/m^2;
#   * step 2, t = 2e-3 s, transient, theta = 0.75: Js_z = 100 A/m^2 (sudden rise: eddy currents in the solid rod);
#   * step 3, t = 3e-3 s, transient, theta = 0.65: Js_z = 0 (switch-off: the eddy currents decay).
# Results: output/two_wires_3D_edge_el_2order_out_<n>.vtu (A, V, B, H, Jt), two_wires_3D_edge_el_2order_steps.csv, two_wires_3D_edge_el_2order.log
#
# How to run:  cd tests/two_wires_3D_edge_el_2order && python macro.py
# Documentation: "Element library" (AV3D_edge), "User manual - gauging edge elements, boundary conditions".
# =====================================================================================================================


from pytec3d import selection, Solver, tree_gauge
import os

# 1. Solver and model files.
pytec3d = Solver()

filename = "two_wires_3D_edge_el_2order"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Tangential A = 0 on the faces x = +-0.5 m: edge values (component 2 of _A) on the midside nodes.
selection.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

#    ... on the faces y = +-0.5 m,
selection.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

#    ... and on the faces z = +-0.5 m.
selection.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Zero potential V on the two end faces of the solid rod: the rod is shorted at its ends.
selection.select(pytec3d.model.mesh, "LOWER_FACE_ROD", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)
selection.select(pytec3d.model.mesh, "UPPER_FACE_ROD", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 4. Tree-cotree gauge: the tree is built over the edges of the SELECTED CELLS (here all of them) and grown from the
#    SELECTED NODES (here the external boundary, where the tangential A was fixed above). It sets _A = 0 on the
#    midside nodes of the tree edges.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
tree_gauge.tree_gauge(pytec3d.model.mesh)

# 5. Load: Js_z = 1 A/m^2 (component 2) on the cells of the source rod.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1, 2, location="elem")

# 6. Step 1: static (transient_effects = 0) with Js_z = 1 A/m^2 (theta is not used on a static step).
pytec3d.time_step["time"] = 1e-3
pytec3d.time_step["theta"] = 0.75
pytec3d.time_step["transient_effects"] = 0

# 7. Load the model and solve the static step.
pytec3d.load_model_into_solver()

pytec3d.solve()
pytec3d.write_output()

# 8. Step 2: transient, Js_z raised to 100 A/m^2, theta = 0.75. Loads can change between steps without reloading.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 100, 2, location="elem")

pytec3d.time_step["time"] = 2e-3
pytec3d.time_step["theta"] = 0.75
pytec3d.time_step["transient_effects"] = 1

pytec3d.solve()
pytec3d.write_output()

# 9. Step 3: transient, source switched off, theta = 0.65.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 2, location="elem")

pytec3d.time_step["time"] = 3e-3
pytec3d.time_step["theta"] = 0.65
pytec3d.time_step["transient_effects"] = 1

pytec3d.solve()
pytec3d.write_output()

# 10. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

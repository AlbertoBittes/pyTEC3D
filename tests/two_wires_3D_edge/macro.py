# =====================================================================================================================
# two_wires_3D_edge - 3D: source wire with imposed Js and passive conducting wire, AV3D_edge
# =====================================================================================================================
#
# What this example shows
# -----------------------
# Two long parallel wires along z in a box: a STRANDED CONDUCTOR with an imposed source current density (Js, element
# load) and a passive massive conductor where eddy currents are induced. One transient step from zero fields.
# Edge (Nedelec) vector potential, first-order Lagrange V. It shows Dirichlet conditions set by planes (select_pos),
# the MidsideNodes point array to restrict a selection to the midside nodes of a quadratic mesh, and the
# element load Js.
# Edge elements without gauge: no tree_gauge() call, so the curl-curl matrix of the air is singular (the gradient
# fields); the direct solver still solves the consistent system. two_wires_3D_edge_el shows the tree-cotree gauge.
#
# Geometry and mesh (two_wires_3D_edge.vtu): -1 <= x <= 1 m, -0.5 <= y <= 0.5 m, -0.5 <= z <= 0.5 m;
#   54156 nodes, the quadratic two_wires_3D mesh (TET10, HEX20, PYRAMID13)
#   needed by the edge elements.
#   * Source wire (TYPE / MATERIAL 23): square cross-section about 0.1 m wide at x = -0.2 m, along z.
#   * Passive wire (TYPE / MATERIAL 24): the same at x = +0.2 m.   * Air (TYPE / MATERIAL 21).
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 2 = AV3D_edge / NON_CONDUCTING            / 1 / 0   air
#   * 23: 2 = AV3D_edge / STRANDED_CONDUCTOR        / 1 / 0   source wire: Js imposed
#   * 24: 2 = AV3D_edge / SOLID_CONDUCTOR           / 1 / 0   passive wire: A + V
#   Materials: 21 mur = diag(1.3, 1, 1); 23 mur = 1; 24 mur = 1, rho = 1e-9 Ohm m.  Real constants 1: none needed
#   (empty row).
#
# Loads and time history: Js_z = 1 A/m^2 in the source wire; one step t = 1e-3 s, theta = 0.5, transient.
# Results: output/two_wires_3D_edge_out_1.vtu (A, V, B, H, Jt), two_wires_3D_edge_steps.csv, two_wires_3D_edge.log
#
# How to run:  cd tests/two_wires_3D_edge && python macro.py
# Documentation: "Element library" (AV3D_edge), "User manual - boundary conditions", "Examples" (two wires).
# =====================================================================================================================


from pytec3d import selection, Solver
import os

# 1. Solver and model files.
pytec3d = Solver()

filename = "two_wires_3D_edge"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Planes x = +-0.5 m (axis 0): tangential A = 0, component 2 of _A on the midside (edge) nodes.
#    Note: the x faces of this mesh are at x = +-1 m; the planes x = +-0.5 m hold (almost) no node, so this
#    selection constrains (almost) nothing.
selection.select_pos(pytec3d.model.mesh, 0, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 0, (-0.501, -0.49999), mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Faces y = +-0.5 m (axis 1): tangential A = 0 (edge values).
selection.select_pos(pytec3d.model.mesh, 1, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 1, (-0.501, -0.49999), mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 4. End faces z = +-0.5 m, crossed by the wires: tangential A = 0 on the midside nodes (edge values).
selection.select_pos(pytec3d.model.mesh, 2, (0.499, 0.50001))
selection.select_pos(pytec3d.model.mesh, 2, (-0.501, -0.49999), mode="a")
selection.select(pytec3d.model.mesh, "MidsideNodes", 1, mode="r")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

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

# =====================================================================================================================
# two_wires_2D - 2D planar: a source wire with imposed current density and a passive conducting wire
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The simplest planar (option3 = 0) 2D eddy-current model: two long parallel wires, normal to the xy plane. The left
# wire is a STRANDED CONDUCTOR with an imposed source current density Js (element load, out-of-plane component); the
# right wire is a massive conductor in which the changing field induces eddy currents. It shows:
#   * a planar 2D model: the out-of-plane direction is z and the real constant thickness is the depth of the model;
#   * an element load (Js, cell array, component 2 = z) set on the cells of a material;
#   * a STATIC first step (transient_effects = 0: no eddy currents, the initial state), followed by a transient step;
#   * an anisotropic material (relative permeability diag(1.3, 1, 1) in the air).
#
# Geometry and mesh (two_wires_2D.vtu): -1 <= x <= 1 m, -0.5 <= y <= 0.5 m, 1212 TRI6 + 32 QUAD8 cells.
#   * Source wire (TYPE / MATERIAL 23): 0.1 m x 0.1 m square centred at x = -0.5 m.
#   * Passive wire (TYPE / MATERIAL 24): 0.1 m x 0.1 m square centred at x = +0.5 m.
#   * Air (TYPE / MATERIAL 21).
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV2D_nodal / NON_CONDUCTING            / 2 / 0   air
#   * 23: AV2D_nodal / STRANDED_CONDUCTOR        / 2 / 0   source wire: Js imposed, no eddy currents
#   * 24: AV2D_nodal / SOLID_CONDUCTOR           / 2 / 0   passive wire: eddy currents, no applied voltage
#   Materials: 21 mur = diag(1.3, 1, 1); 23 mur = 1; 24 mur = 1, rho = 1e-8 Ohm m.
#   Real constants 1: thickness not given: 1 m of depth (the default).
#
# Boundary conditions:  A_z = 0 on the external boundary.
# Loads and time history (theta = 0.91):
#   * step 1, t = 1e-5 s, static: Js_z = 1000 A/m^2 in the source wire;
#   * step 2, t = 4e-5 s, transient (dt = 3e-5 s): Js_z = 2000 A/m^2: the passive wire reacts with eddy currents.
# Results: output/two_wires_2D_out_<n>.vtu (A, B, H, Jt), two_wires_2D_steps.csv, two_wires_2D.log
#
# How to run:  cd tests/two_wires_2D && python macro.py
# Documentation: "Element library" (AV2D_nodal, stranded conductor), "User manual - loads, static and transient steps".
# =====================================================================================================================


from pytec3d import selection, Solver
import os

# 1. Solver and model files.
pytec3d = Solver()

filename = "two_wires_2D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Boundary condition: A_z = 0 (component 2 of _A) on the external boundary: select every cell, then the nodes on
#    the external boundary of the selection.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Load the model into the solver session.
pytec3d.load_model_into_solver()

# 4. Theta method with theta = 0.91 (between Crank-Nicolson, 0.5, and backward Euler, 1).
pytec3d.time_step["theta"] = 0.91

# 5. Time loop. transient_effects = 0 solves a static step (no damping terms): the initial magnetostatic state.
#    The element load Js is a cell array: component 2 (z, out of plane) on the cells of MATERIAL 23.
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

# 6. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

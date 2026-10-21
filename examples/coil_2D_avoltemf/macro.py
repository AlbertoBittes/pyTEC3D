# =====================================================================================================================
# coil_2D_avoltemf - 2D axisymmetric stranded coil (circuit equation: VOLT + EMF) next to a passive massive ring
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A 2D axisymmetric transient analysis of a multi-turn STRANDED COIL fed through its terminals, facing a passive
# massive ring. A stranded coil is a winding of many thin turns: no eddy currents flow inside it, and its current
# density is uniform, J = (N / S) I T, with N turns, S the winding cross-section, T the winding direction and I the coil
# current. I is not imposed point by point: it follows from the circuit equation of the coil,
#       VOLT = R I + EMF,      EMF = d(Psi)/dt (the voltage induced by the flux linked with the winding),
# solved together with the field. The two circuit unknowns VOLT and EMF are DOFs of the coil nodes, coupled into one
# value each. It shows:
#   * the real constants of a stranded coil (cross_section, turns, radius, resistance, symmetry, T);
#   * a DOF coupling of two fields at once (CP 1 -> ["EMF", "VOLT"]);
#   * a current-driven step (nodal load CURR = coil current) followed by a voltage-driven step (_VOLT);
#   * a passive SOLID_CONDUCTOR: a massive ring without terminal (no electric DOF in 2D), shorted by construction.
#
# Geometry and mesh (coil_2D_avoltemf.vtu): axisymmetric about the y axis (x = r, y = z), 0 <= r <= 1 m,
#   -1 <= z <= 1 m, 800 QUAD8 cells. Lower ring TYPE / MATERIAL 23 (the coil) at r = 0.5 m, z = -0.5 m; upper ring
#   TYPE / MATERIAL 24 (passive conductor) at r = 0.5 m, z = +0.5 m; both 0.2 m x 0.2 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV2D_nodal / NON_CONDUCTING  / 2 / 1   air
#   * 23: AV2D_nodal / STRANDED_COIL   / 2 / 1   coil: A_z + VOLT (terminal voltage) + EMF (induced voltage)
#   * 24: AV2D_nodal / SOLID_CONDUCTOR / 2 / 1   passive ring: A_z only (eddy currents)
#   Materials: 21 air (mur = 1); 23 mur = 1 (its rho is not used: a stranded coil has no eddy currents);
#              24 mur = 1, rho = 1e-8 Ohm m.
#   Real constants (CONSTANT 1, used by every cell; only the coil reads the coil constants):
#     thickness = 1 (full ring), cross_section S = 0.04 m^2, turns N = 1000, radius = 0.5 m (mean radius of the coil:
#     coil volume = 2 pi radius S thickness), resistance R = 50 Ohm, symmetry = 1 (the whole coil is modelled),
#     T = (0, 0, 1) winding direction (z = out of plane = toroidal).
#   DOF couplings (set by the macro, add_dof_couplings): 1 -> ["EMF", "VOLT"] (the coil).
#
# Boundary conditions:  A_z = 0 on the external boundary;
#                       step 2: coil terminal voltage VOLT = 10 kV.
# Loads and time history (theta = 0.5, dt = 1e-5 s):
#   * step 1: coil current CURR = 100 A (the voltage needed to drive it is computed);
#   * step 2: CURR = 0 and VOLT = 10 kV (the current is computed from VOLT = R I + EMF).
# Results: output/coil_2D_avoltemf_out_<n>.vtu (A, VOLT, EMF, B, H, Jt), coil_2D_avoltemf_steps.csv, .log
#
# How to run:  cd examples/coil_2D_avoltemf && python macro.py
# Documentation: "Element library" (stranded coil), "Theory - stranded coils (circuit coupling)".
# =====================================================================================================================


from pytec3d import selection, Solver
import os
import vtk.util.numpy_support as vtknp

# 1. Solver and model files (mesh + model data).
pytec3d = Solver()

filename = "coil_2D_avoltemf"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")
# DOF couplings (set here, not in the model data): the nodes with CP_0 = 1 share one EMF and one VOLT value.
pytec3d.model.add_dof_couplings({1: ["EMF", "VOLT"]})

# 2. Boundary condition: A_z = 0 on the external boundary.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Terminals of the stranded coil: one VOLT and one EMF for the whole winding
#    CP_0 = 1 on the coil nodes; the coupling 1 -> ["EMF", "VOLT"] (add_dof_couplings) ties both fields to the master node (the lowest
#    node id), where the coil current (CURR) or the terminal voltage (_VOLT) is applied.
# ---------------------------------------------------------------------------------------------------------------------
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

# ---------------------------------------------------------------------------------------------------------------------
# 4. Passive ring: a SOLID_CONDUCTOR, a massive conductor without electric DOF in 2D (J = -sigma dA/dt): a closed,
#    shorted ring in which only the induced (eddy) current flows. Nothing to set.
# ---------------------------------------------------------------------------------------------------------------------

# 5. Load the model into the solver session; theta method with theta = 0.5.
pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 0.5

# ---------------------------------------------------------------------------------------------------------------------
# 6. Time loop
#    CURR, on the master node of the coil, is the coil current I (times the symmetry factor, 1 here).
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(pytec3d.model.mesh)
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        # Current driven: I = 100 A.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        # Voltage driven: terminal voltage 10 kV, I = (VOLT - EMF) / R.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 7. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

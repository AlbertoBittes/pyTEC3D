# =====================================================================================================================
# coil_2D_avolt_timint - coil_2D_avolt with the time-integrated electric potential (AV2D_nodal_timint)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The same analysis as coil_2D_avolt (read that macro first), solved with the element type AV2D_nodal_timint. In this
# formulation the electric DOF of the circuit-coupled conductor is the TIME INTEGRAL of the terminal voltage,
# v(t) = int_0^t V dt  [V s], so that the electric terms move from the stiffness K to the damping C and the coupled
# A-v system is symmetric. The script is identical; what changes is the meaning of the values:
#   * AMPS is still the total current through the conductor cross-section [A];
#   * _V now fixes the time-integrated voltage v [V s], not the voltage: step 2 imposes v(t = 2e-5 s) = 10000 V s.
# Compare the outputs of the two macros to see how the two formulations of the same terminal differ.
# The time-integrated DOFs have no row in K: on a static step (transient_effects = False) the solver fixes them to
# zero, so a time-integrated conductor carries no conduction current on a static step.
#
# Geometry and mesh (coil_2D_avolt_timint.vtu): as coil_2D_avolt (axisymmetric about y, 800 QUAD8 cells; lower ring
#   TYPE 23 at z = -0.5 m, upper ring TYPE 24 at z = +0.5 m, r = 0.5 m, 0.2 m x 0.2 m cross-sections).
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 3 = AV2D_nodal_timint / NON_CONDUCTING            / 2 / 1   air
#   * 23: 3 = AV2D_nodal_timint / CIRCUIT_COUPLED_CONDUCTOR / 2 / 1   active ring: A_z + time-integrated V
#   * 24: 3 = AV2D_nodal_timint / SOLID_CONDUCTOR           / 2 / 1   passive shorted ring (eddy currents)
#   Materials: 21 air, 23 rho = 1e-7, 24 rho = 1e-8 (mur = 1).  Real constants: 1: thickness = 1.
#   DOF couplings: 1 -> ["V"].
#
# Boundary conditions:  A_z = 0 on the external boundary; step 2: time-integrated voltage v = 10000 V s.
# Loads and time history (theta = 0.5, dt = 1e-5 s): step 1 AMPS = 100 A; step 2 AMPS = 0 and _V = 10000.
# Results: output/coil_2D_avolt_timint_out_<n>.vtu, coil_2D_avolt_timint_steps.csv, coil_2D_avolt_timint.log
#
# How to run:  cd tests/coil_2D_avolt_timint && python macro.py
# Documentation: "Element library" (AV2D_nodal_timint), "Theory - time-integrated potentials".
# =====================================================================================================================


from pytec3d import mesh, selection, Solver
import os
import vtk
import vtk.util.numpy_support as vtknp

# 1. Solver and model files.
pytec3d = Solver()

filename = "coil_2D_avolt_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# Make sure the COORDINATE_SYSTEM cell array exists (no-op when present; 0 = global Cartesian system).
mesh.add_array(pytec3d.model.mesh, "COORDINATE_SYSTEM", location="elem", dtype=vtk.VTK_INT)

# 2. Boundary condition: A_z = 0 (component 2 of _A) on the external boundary.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Terminal of the circuit-coupled conductor: couple the V DOFs of the ring (CP id 1) and keep the master node.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

# 4. Load the model into the solver session (the model is now locked).
pytec3d.load_model_into_solver()

# 5. Theta method, theta = 0.5.
selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 0.5

# 6. Time loop: step 1 current driven, step 2 driven by the time-integrated voltage.
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "AMPS", 100)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "AMPS", 0)
        # In this formulation _V fixes v = int V dt [V s].
        selection.set_value_on_selection(pytec3d.model.mesh, "_V", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 7. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

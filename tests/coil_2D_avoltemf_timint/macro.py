# =====================================================================================================================
# coil_2D_avoltemf_timint - coil_2D_avoltemf with time-integrated circuit unknowns (AV2D_nodal_timint)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The stranded-coil analysis of coil_2D_avoltemf (read that macro first) with the element type AV2D_nodal_timint:
# the circuit unknowns of the coil are TIME-INTEGRATED quantities,
#   VOLT -> int_0^t VOLT dt [V s],   EMF -> int_0^t EMF dt = Psi, the flux linked with the coil [Wb],
# so that every electric term is in the damping matrix C and the system is symmetric. Consequences for the input:
#   * CURR is still the coil current [A];
#   * _VOLT = 10000 fixes the time integral of the terminal voltage (10000 V s at t = 2e-5 s), not the voltage;
# Backward Euler (theta = 1) is used.
#
# Geometry and mesh (coil_2D_avoltemf_timint.vtu): as coil_2D_avoltemf (coil TYPE 23 at z = -0.5 m, passive ring
#   TYPE 24 at z = +0.5 m, r = 0.5 m, axisymmetric about y, 800 QUAD8 cells).
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV2D_nodal_timint / NON_CONDUCTING  / 2 / 1   air
#   * 23: AV2D_nodal_timint / STRANDED_COIL   / 2 / 1   coil: A_z + int VOLT + int EMF
#   * 24: AV2D_nodal_timint / SOLID_CONDUCTOR / 2 / 1   passive ring: A_z only (eddy currents)
#   Materials: 21 air; 23 coil; 24 rho = 1e-8 Ohm m (mur = 1 everywhere).
#   Real constants 1: thickness = 1, cross_section = 0.04 m^2, turns = 1000, radius = 0.5 m, resistance = 50 Ohm,
#     symmetry = 1, T = (0, 0, 1).
#   DOF couplings: 1 -> ["EMF", "VOLT"] (coil).
#
# Boundary conditions:  A_z = 0 on the external boundary; step 2: coil int VOLT = 10000 V s.
# Loads and time history (theta = 1.0, dt = 1e-5 s): step 1 CURR = 100 A; step 2 CURR = 0 and _VOLT = 10000.
# Results: output/coil_2D_avoltemf_timint_out_<n>.vtu, coil_2D_avoltemf_timint_steps.csv, .log
#
# How to run:  cd tests/coil_2D_avoltemf_timint && python macro.py
# Documentation: "Element library" (AV2D_nodal_timint), "Theory - stranded coils", "Theory - time integration".
# =====================================================================================================================


from pytec3d import selection, Solver
import os
import vtk.util.numpy_support as vtknp

# 1. Solver and model files.
pytec3d = Solver()

# convert_mesh_from_cdb(
#     r"D:\pyTEC3D\simulations\single_wire_edge_elements\single_wire_edge_elements.cdb"
# )
# (An ANSYS .cdb mesh is converted once with the command line tool:  pytec3d-convert-cdb mesh.cdb)
filename = "coil_2D_avoltemf_timint"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Boundary condition: A_z = 0 on the external boundary.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Coil terminals: couple VOLT and EMF of the winding (CP id 1); master node = lowest node id.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

# 4. Passive ring: a SOLID_CONDUCTOR (massive, no electric DOF in 2D): a shorted ring, nothing to set.

# 5. Load the model; backward Euler.
pytec3d.load_model_into_solver()
pytec3d.time_step["theta"] = 1.0
selection.select_all(pytec3d.model.mesh)

# 6. Time loop: current driven, then driven by the time-integrated terminal voltage.
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 7. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

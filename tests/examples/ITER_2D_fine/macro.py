# =====================================================================================================================
# ITER_2D_fine - free-boundary plasma EQUILIBRIUM of ITER in 2D axisymmetry (Grad-Shafranov, Picard iterations)
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The computation of a free-boundary plasma equilibrium with the Algorithm class (APEC-M). The magnetostatic field of
# the coils and of the plasma current is solved with the finite elements; around it, Algorithm.equilibrium_solve()
# iterates (Picard) on the plasma current density until the poloidal flux PSI stops changing. Each iteration:
#   1. finds the magnetic axis (PSI maximum), the X point and the limiter point, hence the last closed flux surface;
#      the cells inside it form the plasma (PLASMA cell array);
#   2. computes the plasma current density Js from the Grad-Shafranov profile
#        J_phi = lambda (1 - psi_n^alpha_mp)^alpha_np (beta_p R / R0 + (1 - beta_p) R0 / R),
#      lambda scaling the total to plasma_current;
#   3. corrects the currents of the control coils (radial and vertical correction currents, with the gains of
#      control_coils) so that the poloidal field at the target magnetic axis vanishes;
#   4. solves the field with the new Js and recomputes PSI;
#   5. stops when the largest change of PSI is below tolerance_value.
# The first solve uses the initial guess set by the macro: a uniform current density over the whole vacuum chamber.
# The result is the starting point of the evolution cases (ITER_2D_evol*, ITER_3D_evol*), which read it from
# ./output.
#
# Geometry and mesh (ITER_2D_fine.vtu): 2D axisymmetric about the y axis (x = R, y = Z), 0 <= R <= 100 m,
#   -100 <= Z <= 100 m, 17818 QUAD4 + 16 TRI3 cells (finer first-order mesh). TYPE ids used by the cells: 21 (air), 23 (coils, plasma region and structures).
#   Named cell arrays (1 on the cells of a component): CS1U, CS1L, CS2U, CS2L, CS3U, CS3L (central solenoid modules),
#   PFC_1 ... PFC_6 (poloidal field coils), COILS, VACUUM_CHAMBER (where the plasma can be), PASSIVE_STRUCTURES.
#
# Model data (model_data.json)
#   Element types (AV2D_nodal, option2 = 1, option3 = 1 axisymmetric): 21 NON_CONDUCTING, 23 STRANDED_CONDUCTOR
#   (current density Js imposed), 24 SOLID_CONDUCTOR, 25 STRANDED_COIL; the rows not used by any cell are
#   reported as warnings when the model is checked. Every step here is static, so the conductivities play no role.
#   Materials: mur = 1 everywhere (rho of the structures and coils given for the evolution cases).
#   Real constants: one row per coil (cross_section, turns, radius, resistance, symmetry, T) and rows for the other
#   regions; a stranded conductor only uses the out-of-plane thickness.
#
# Boundary conditions:  A_z = 0 on the external boundary (axis R = 0 and the far boundary).
# Loads:  Js_z = I / a in each CS module and PF coil (the coil current I in ampere-turns over its cross-section a);
#   initial plasma guess: Js_z = 15 MA / a_plasma over the vacuum chamber. One static solve per Picard iteration.
# Algorithm parameters:  target axis (6.52 m, 0.54 m), R0 = 6.2 m, plasma_current = 15 MA, beta_p = 0.682,
#   alpha_mp = 1.0001, alpha_np = 0.83197, control coils PFC_2 ... PFC_5 with gains +-0.8, field of a unit
#   correction current at the target axis (bradx/brady/bverx/bvery_value), tolerance 1e-6, up to 70 iterations.
# input_parameters["min_distance"] is not a parameter of Algorithm: it is ignored.
#
# Results: output/ITER_2D_fine_out_<n>.vtu, written once after convergence, with <n> the number of solves done (one per
#   Picard iteration): A, PSI, Js (plasma and coils), PLASMA, B, H, Jt. The evolution macros
#   read the most recent file of this series. __pytec3d__/pytec3d_iteration_diagnostics.csv: one row per Picard iteration (axis,
#   X point, limiter, flux change, correction currents, plasma current scaling). ITER_2D_fine.log: the run report.
#
# How to run:  cd tests/ITER_2D_fine && python macro.py      (a long run: Picard iterations on the full machine)
# Documentation: "User manual - plasma equilibrium and evolution (Algorithm)", "Examples" (ITER).
# =====================================================================================================================


from pytec3d import Algorithm, selection, Solver
import os

# from convert_mesh import convert_mesh_from_cdb

# 1. The solver and the plasma algorithm that drives it (Algorithm only uses the public API of the solver).
pytec3d = Solver()
algorithm = Algorithm(pytec3d)

aplasma = 0

# 2. Read the model (mesh with the named component arrays).
filename = "ITER_2D_fine"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

# 3. Parameters of the equilibrium (see Algorithm.DEFAULT_INPUT_PARAMETERS):
#    target position of the magnetic axis (x = R, y = Z), major radius R0, Picard iterations, control coils and
#    their gains, convergence tolerance on PSI, field at the target axis of a unit radial / vertical correction
#    current, Grad-Shafranov profile exponents, plasma current [A] and poloidal beta.
algorithm.input_parameters["target_central_axis_x"] = 6.52136102486433
algorithm.input_parameters["target_central_axis_y"] = 0.536537145137914
algorithm.input_parameters["target_central_axis_z"] = 0.0
algorithm.input_parameters["R0"] = 6.2
algorithm.input_parameters["max_picard_iterations"] = 70
algorithm.input_parameters["min_distance"] = 0.6
algorithm.input_parameters["control_coils"] = [
	{"component": "PFC_2", "radial_gain": 0.8, "vertical_gain": 0.0},
	{"component": "PFC_3", "radial_gain": 0.0, "vertical_gain": 0.8},
	{"component": "PFC_4", "radial_gain": 0.0, "vertical_gain": 0.8},
	{"component": "PFC_5", "radial_gain": -0.8, "vertical_gain": 0.0},
]
algorithm.input_parameters["tolerance_value"] = 1e-6
algorithm.input_parameters["bradx_value"] = -3.873415861e-08
algorithm.input_parameters["brady_value"] = 5.313857037e-09
algorithm.input_parameters["bverx_value"] = 1.822467692e-10
algorithm.input_parameters["bvery_value"] = 1.133287835e-07
algorithm.input_parameters["alpha_mp"] = 1.0001
algorithm.input_parameters["alpha_np"] = 0.83197
algorithm.input_parameters["plasma_current"] = 15e6
algorithm.input_parameters["beta_p"] = 6.8216e-001

# 4. Coil currents [A] (total current through the cross-section of each CS module and PF coil).
i_CS3U = -4.8204e06
i_CS2U = -4.0480e06
i_CS1U = -2.0804e07
i_CS1L = -2.0804e07
i_CS2L = -1.0297e07
i_CS3L = 5.3012e06

i_PF1 = 5.7325e06
i_PF2 = -2.8834e06
i_PF3 = -5.8084e06
i_PF4 = -4.7762e06
i_PF5 = -7.8133e06
i_PF6 = 1.6941e07

# Cross-section areas [m^2] of the PF coils, of a CS module and of the vacuum chamber region.
a_PF1 = 0.94375190
a_PF2 = 0.41453946
a_PF3 = 0.66413094
a_PF4 = 0.60871516
a_PF5 = 0.77496250
a_PF6 = 1.72659250
a_CS = 1.55608
a_plasma = 28.6027095

# 5. Model data: element types, materials, real constants (see the header).
pytec3d.model.read_json("model_data.json")

# 6. Boundary condition: A_z = 0 on the external boundary (axis R = 0 and far boundary). All the cells use the
#    global coordinate system.
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 7. Static analysis: every Picard iteration is a magnetostatic solve (transient_effects False).
pytec3d.time_step["time"] = 1.0
pytec3d.time_step["theta"] = 0.5
pytec3d.time_step["transient_effects"] = 0

# 8. Coil loads: Js_z = I / a (component 2 = toroidal) on the cells of each CS module ...
selection.select(pytec3d.model.mesh, "CS3L", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS3L / a_CS, 2, location="elem")
selection.select(pytec3d.model.mesh, "CS2L", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS2L / a_CS, 2, location="elem")
selection.select(pytec3d.model.mesh, "CS1L", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS1L / a_CS, 2, location="elem")
selection.select(pytec3d.model.mesh, "CS1U", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS1U / a_CS, 2, location="elem")
selection.select(pytec3d.model.mesh, "CS2U", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS2U / a_CS, 2, location="elem")
selection.select(pytec3d.model.mesh, "CS3U", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_CS3U / a_CS, 2, location="elem")

#    ... and of each PF coil.
selection.select(pytec3d.model.mesh, "PFC_1", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF1 / a_PF1, 2, location="elem")
selection.select(pytec3d.model.mesh, "PFC_2", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF2 / a_PF2, 2, location="elem")
selection.select(pytec3d.model.mesh, "PFC_3", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF3 / a_PF3, 2, location="elem")
selection.select(pytec3d.model.mesh, "PFC_4", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF4 / a_PF4, 2, location="elem")
selection.select(pytec3d.model.mesh, "PFC_5", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF5 / a_PF5, 2, location="elem")
selection.select(pytec3d.model.mesh, "PFC_6", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", i_PF6 / a_PF6, 2, location="elem")

# 9. Initial guess of the plasma current: 15 MA spread uniformly over the vacuum chamber.
selection.select(pytec3d.model.mesh, "VACUUM_CHAMBER", 1, location="elem")
selection.set_value_on_selection(pytec3d.model.mesh, "Js", 15e6 / a_plasma, 2, location="elem")

# 10. Load the model into the solver session.
pytec3d.load_model_into_solver()

# 11. Picard iterations of the free-boundary equilibrium (raises ConvergenceError if not converged).
algorithm.equilibrium_solve()

# 12. Write the converged equilibrium to ./output.
pytec3d.write_output()

# 13. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

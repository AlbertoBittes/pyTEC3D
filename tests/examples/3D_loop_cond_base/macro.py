# =====================================================================================================================
# 3D_loop_cond_base - 3D quarter model: square source loop (imposed Js) below an identical passive solid loop
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A 3D transient analysis of a passive conducting square loop placed above an identical active loop that carries an
# imposed source current density. Only one quarter of the loops and of the surrounding air is modelled, using the
# two symmetry planes x = 0 and y = 0. It shows:
#   * boundary conditions selected through named point arrays stored in the mesh (X_BOUNDS_1, Z_BOUNDS_2, ...);
#   * far-field (A = 0) and symmetry conditions (tangential A = 0) for a nodal vector potential;
#   * the electric potential V of a solid conductor fixed where the loop is cut by the symmetry planes;
#   * a vector element load Js given component by component on named groups of cells (the straight sides and the
#     45 degree corner of the loop);
#   * loads and boundary conditions set AFTER the model is loaded into the solver: once loaded, the mesh and the model
#     data are locked but the values of the arrays may change at any time;
#   * a static first step (initial magnetostatic state) followed by a transient step.
#
# Geometry and mesh (3D_loop_cond_base.vtu)
# -----------------------------------------
#   Quarter domain 0 <= x <= 0.25 m, 0 <= y <= 0.25 m, -0.125 <= z <= 0.125 m, 15625 HEX8 cells (25 x 25 x 25).
#   * Source loop (TYPE / MATERIAL 23, 9 cells): quarter of a square loop in the plane z = -0.02 m; cell groups
#     COIL_NOCOND_1 (side along y), COIL_NOCOND_2 (side along x), COIL_NOCOND_3 (corner).
#   * Passive loop (TYPE / MATERIAL 24, 9 cells): the same shape at z = +0.02 m (cell groups COIL_COND_1..3).
#   * Air (TYPE / MATERIAL 21).
#   Named point arrays (1 on the nodes of a face): X_BOUNDS_1 (x = 0), X_BOUNDS_2 (x = 0.25), Y_BOUNDS_1 (y = 0),
#   Y_BOUNDS_2 (y = 0.25), Z_BOUNDS_1 / Z_BOUNDS_2 (z = -0.125 / +0.125), VOLT_0_1 / VOLT_0_2 (the nodes of the
#   passive loop on the two symmetry planes).
#
# Model data (model_data.json)
# ----------------------------
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 1 = AV3D_nodal / 0 = NON_CONDUCTING     / 1 / 0   air: nodal A with Coulomb gauge term
#   * 23: 1 = AV3D_nodal / 3 = STRANDED_CONDUCTOR / 1 / 0   source loop: Js imposed, no eddy currents
#   * 24: 1 = AV3D_nodal / 5 = SOLID_CONDUCTOR    / 1 / 0   passive loop: A + V (non-integrated scalar potential)
#   Materials (MATERIAL id): mur relative permeability, rho resistivity [Ohm m]
#   * 21: mur = 1 (rho not used)   * 23: mur = 1 (rho not used)   * 24: mur = 1, rho = 1e-9
#   Real constants 1: none needed (empty row).   No DOF coupling.
#
# Boundary conditions
# -------------------
#   * Far field, A = 0: faces z = +-0.125 m, x = 0.25 m and y = 0.25 m.
#   * Symmetry plane x = 0: A_y = A_z = 0 (A, like the current, is normal to the plane; B is tangent to it).
#   * Symmetry plane y = 0: A_x = A_z = 0.
#   * V = 0 on the nodes of the passive loop on the symmetry planes (the potential is odd with respect to them); this
#     also fixes the free constant of V.
#
# Loads and time history (theta = 1, the default)
# ------------------------------------------------
#   * Step 1, t = 1e-4 s, static: Js of magnitude 1 A/m^2 along the source loop (counter-clockwise about z).
#   * Step 2, t = 2e-4 s, transient (dt = 1e-4 s): Js halved: eddy currents are induced in the passive loop.
#
# Results (./output): 3D_loop_cond_base_out_<n>.vtu (A, V, B, H nodal and element, Jt), 3D_loop_cond_base_steps.csv,
#   3D_loop_cond_base.log.
#
# How to run:  cd tests/3D_loop_cond_base && python macro.py
# Documentation: "Examples" (3D loops), "Element library" (AV3D_nodal), "User manual - boundary conditions, loads".
# =====================================================================================================================


# 0. Import necessary modules
from pytec3d import selection, Solver
import os

# 1. Initialize the program instance: one Solver per analysis (model, time step, solver settings, solver session).
pytec3d = Solver()

# 2. Define the file name and path for the model and read it into the program.
#    read_model() reads the VTK grid, names the model after the file and adds the arrays the element types need.
filename = "3D_loop_cond_base"
filepath = os.path.join(os.getcwd(), filename + ".vtu")
pytec3d.read_model(filepath)

# 3. Read the model data (element types, materials, real constants, DOF couplings, nonlinear laws, coordinate
#    systems) from the JSON file.
pytec3d.model.read_json("model_data.json")

# 4. Save the mesh with the model data stored in it (VTK field data) to 3D_loop_cond_base.vtu: the file can then be
#    read back alone. Run from the case folder, this overwrites the input file.
pytec3d.save_model()

# 5. Load the model into the solver module: validation, then the C++ session. From now on only array VALUES change.
pytec3d.load_model_into_solver()

# 6. Boundary conditions (Dirichlet values in the point array _A: NaN = free, a number = fixed; component 0, 1, 2 =
#    x, y, z).
# Far-field boundary conditions: set magnetic vector potential A to zero
selection.select(pytec3d.model.mesh, "Z_BOUNDS_1", 1)
selection.select(pytec3d.model.mesh, "Z_BOUNDS_2", 1, mode="a")
selection.select(pytec3d.model.mesh, "Y_BOUNDS_2", 1, mode="a")
selection.select(pytec3d.model.mesh, "X_BOUNDS_2", 1, mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# Symmetry boundary conditions: set the components of A parallel to the symmetry planes to zero
selection.select(pytec3d.model.mesh, "X_BOUNDS_1", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

selection.select(pytec3d.model.mesh, "Y_BOUNDS_1", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# Voltage boundary conditions: set the electric potential V of the passive loop on the symmetry planes
selection.select(pytec3d.model.mesh, "VOLT_0_1", 1)
selection.select(pytec3d.model.mesh, "VOLT_0_2", 1, mode="a")
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# sqrt(2)/2 factor for the corner element of the current loop (current at 45 degrees)
b = 2**0.5 / 2

# 7. Transient analysis: loop over the number of desired time-steps
for i in range(1, 3):
    # Set the current time for this time-step (the END time of the step)
    pytec3d.time_step["time"] = i * 0.0001
    # Set the transient effects flag
    # initial time-step is static to establish the initial conditions
    if i == 1:
        pytec3d.time_step["transient_effects"] = 0
    # subsequent time-steps have transient effects
    else:
        pytec3d.time_step["transient_effects"] = 1

    # 8. Loads: impose the source current density Js (cell array, one component at a time) on the source loop,
    #    magnitude 1/i A/m^2: +y on side 1, -x on side 2, (-b, +b) on the corner.
    selection.select(pytec3d.model.mesh, "COIL_NOCOND_1", 1, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 1 / i, 1, location="elem")
    selection.select(pytec3d.model.mesh, "COIL_NOCOND_2", 1, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", -1 / i, 0, location="elem")
    selection.select(pytec3d.model.mesh, "COIL_NOCOND_3", 1, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", -b / i, 0, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", b / i, 1, location="elem")

    # 9. Solve the system for this time-step (the solution is written back into the mesh arrays)
    pytec3d.solve()

    # 10. Write the output for this time-step: output/3D_loop_cond_base_out_<i>.vtu
    pytec3d.write_output()

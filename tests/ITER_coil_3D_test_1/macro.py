# =====================================================================================================================
# ITER_coil_3D_test_1 - ITER_coil_3D with named boundary arrays and Js computed cell by cell in global components
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The coil system of ITER_coil_3D (read that macro first) set up in the most explicit way:
#   * the symmetry faces are selected through point arrays stored in the mesh (X_BOUND, Z_BOUND), and their in-plane
#     components are global components, so no nodal frame rotation is needed;
#   * the toroidal source current density is computed by hand from the cell centroids (cell_centroids array) and
#     written directly into the global components of Js with NumPy: the NumPy view returned by vtk_to_numpy writes
#     into the mesh array (call Modified() on the VTK array when the solver must see a change made this way; the
#     package functions do it for you).
#
# Geometry and mesh (ITER_coil_3D_test_1.vtu): 90 degree sector of a sphere of radius 1 m about the y axis,
#   -1 <= z <= 0 m, 10 HEX8, 1275 TET4, 40 PYRAMID5 cells. Source ring TYPE / MATERIAL 23 at r = 0.25 m, y = +0.25 m;
#   solid ring TYPE / MATERIAL 24 at y = -0.25 m. Point arrays: X_BOUND (plane x = 0), Z_BOUND (plane z = 0),
#   AXIS_NODES, INFINITY_BOUNDARY_INNER (sphere r = 1 m), ROD_1, ROD_2.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV3D_nodal / NON_CONDUCTING     / 1 / 0   air
#   * 23: AV3D_nodal / STRANDED_CONDUCTOR / 1 / 0   source ring: Js imposed
#   * 24: AV3D_nodal / SOLID_CONDUCTOR    / 1 / 0   solid ring: A + V
#   Materials: mur = 1 everywhere; 24 rho = 1e-9 Ohm m.
#
# Boundary conditions:  x = 0: A_y = A_z = 0;  z = 0: A_x = A_y = 0;  axis: A = 0;  sphere r = 1 m: A = 0;
#   V = 0 on ROD_1 / ROD_2.
# Loads and time history (theta = 1): step 1, t = 1e-4 s, static; steps 2-9, t = i * 1e-4 s, transient. Toroidal
#   Js of magnitude i**6 A/m^2 at step i: J = i**6 (z / rho, 0, -x / rho), rho = sqrt(x^2 + z^2).
# Results: output/ITER_coil_3D_test_1_out_<n>.vtu, ITER_coil_3D_test_1_steps.csv, ITER_coil_3D_test_1.log
#
# How to run:  cd tests/ITER_coil_3D_test_1 && python macro.py
# Documentation: "User manual - boundary conditions, loads, working with the arrays in NumPy".
# =====================================================================================================================


from pytec3d import selection, Solver
import os
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and model files.
pytec3d = Solver()

filename = "ITER_coil_3D_test_1"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Symmetry plane x = 0: the in-plane components A_y, A_z vanish.
selection.select(pytec3d.model.mesh, "X_BOUND", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Symmetry plane z = 0: A_x, A_y vanish.
selection.select(pytec3d.model.mesh, "Z_BOUND", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

# 4. Axis: A = 0.
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 5. Zero potential V at the rod ends (ROD_1, ROD_2): fixes the potential of the solid ring.
selection.select(pytec3d.model.mesh, "ROD_1", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)
selection.select(pytec3d.model.mesh, "ROD_2", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 6. Far field: A = 0 on the sphere r = 1 m.
selection.select(pytec3d.model.mesh, "INFINITY_BOUNDARY_INNER", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 7. Load the model; backward Euler.
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 1.0

# 8. Time loop: step 1 static, then transient steps. time_step_number is ignored by the solver.
for i in range(1, 10):
    pytec3d.time_step["time_step_number"] = i
    if i == 1:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 1.0  # autots,off (keep as before)
    else:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 1.0  # autots,off (default)

    # select coil elements
    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    # Toroidal Js computed per cell from its centroid, written into the global components 0 (x) and 2 (z).
    selection_cells = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetCellData().GetArray("selection_cell")
    )
    selected_ids = np.where(selection_cells == 1)[0]
    centroids_array = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetCellData().GetArray("cell_centroids")
    )
    js_array = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetCellData().GetArray("Js"))
    for l in selected_ids:
        # Get the center coordinates of the element (el_num)
        # This depends on your mesh structure; assuming pytec3d.model.mesh.get_element_center(el_num) returns (x, z)
        x_coord = centroids_array[l, 0]
        z_coord = centroids_array[l, 2]
        module_square = x_coord**2 + z_coord**2
        module = np.sqrt(module_square)
        cos_theta = z_coord / module
        sin_theta = x_coord / module
        fact = i**6
        js_x = cos_theta * fact
        js_z = -sin_theta * fact
        # Apply Js vector to the element (cell_data)
        js_array[l, 0] = js_x
        js_array[l, 2] = js_z
    # Solve the step and write output/ITER_coil_3D_test_1_out_<i>.vtu.
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

# =====================================================================================================================
# ITER_coil_inf_3D - 3D 90 degree sector: source ring (Js in cylindrical components) and solid ring, rotated frames
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The 3D version of ITER_coil_inf_2D on a 90 degree sector, with the same sequence of static and transient steps.
# It shows:
#   * symmetry conditions on the two sector faces through ROTATED NODAL FRAMES: node_angles (thxy, thyz, thzx, in
#     degrees) rotates the frame in which the Dirichlet values _A of a node are given, so that "component 0 / 1 / 2"
#     are directions chosen by the user (here the two in-plane directions of each sector face);
#   * a source current density given in CYLINDRICAL components (Js[1] = J_theta, toroidal) and converted to global
#     components by model_edit.map_vector_field_to_global_components() before each solve (compare ITER_coil_3D,
#     which gives the same load through an element coordinate system instead);
#   * selections combined in cylindrical and spherical coordinate systems (mode "u" removes from the selection).
#
# Geometry and mesh (ITER_coil_inf_3D.vtu):
#   90 degree sector of a sphere of radius 1 m about the y axis (Y-axis cylindrical system: theta from +x towards
#   -z), 12 HEX8, 534 WEDGE6, 48 TET4, 78 PYRAMID5 cells, 0 <= z <= 1 m (theta from -90 to 0). Point arrays: AXIS_NODES (the nodes on the y axis), ROD_1, ROD_2 (end nodes of the conductors).
#   * Source ring (TYPE / MATERIAL 23): small cross-section at r = 0.25 m, y = +0.25 m.
#   * Solid ring (TYPE / MATERIAL 24): the same at y = -0.25 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV3D_nodal / NON_CONDUCTING     / 1 / 0   air
#   * 23: AV3D_nodal / STRANDED_CONDUCTOR / 1 / 0   source ring: Js imposed
#   * 24: AV3D_nodal / SOLID_CONDUCTOR    / 1 / 0   solid ring: A + V
#   Materials: 21 mur = diag(1.1, 1.2, 1.3) (anisotropic air, to test the rotations); 23 mur = 1; 24 rho = 5e-9 Ohm m.
#
# Boundary conditions
#   * Face theta = 0 (plane z = 0): nodal frames rotated by thyz = 90 degrees; components 0 and 2 (the in-plane
#     directions x and y) fixed to 0: A is normal to the face (toroidal), B tangent to it.
#   * Face theta = -90 degrees (plane x = 0): frames rotated by (thxy, thyz) = (90, -90) degrees; components 0 and 1
#     (the in-plane directions y and z) fixed to 0.
#   * Axis (AXIS_NODES): A_x = A_z = 0.   * Outer sphere r = 1 m: A = 0.
#   * V = 0 on ROD_1 / ROD_2 (ends of the solid ring).
# Loads and time history (9 steps; Js magnitude i**6 A/m^2 at step i, an ever faster rising current):
#   * step 1, t = 1e-4 s: static;   steps 2-5, t = 2e-4 ... 5e-4 s: transient, theta = 0.75;
#   * step 6, t = 1e-3 s: transient, theta = 1.0;   step 7, t = 2e-3 s: transient, theta = 0.8;
#   * step 8, t = 3e-3 s: STATIC again (eddy currents switched off, theta unused);   step 9, t = 4e-3 s: transient,
#     theta = 0.9.  The # TIMINT / tintp comments of the loop are the equivalent ANSYS commands.
#   time_step["time_step_number"] is not a time-step setting: the solver ignores it.
# Results: output/ITER_coil_inf_3D_out_<n>.vtu, ITER_coil_inf_3D_steps.csv, ITER_coil_inf_3D.log
#
# How to run:  cd tests/ITER_coil_inf_3D && python macro.py
# Documentation: "User manual - boundary conditions (rotated nodal frames), coordinate systems, loads".
# =====================================================================================================================


from pytec3d import model_edit, selection, Solver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and model files.
pytec3d = Solver()

filename = "ITER_coil_inf_3D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Face theta = 0, without the axis nodes and the outer sphere: rotate the nodal frames (thyz = 90 degrees),
#    then fix the two in-plane components (local 0 and 2).
selection.select_pos(pytec3d.model.mesh, 1, (-0.0001, 0.0001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1, mode="u")
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", 90.0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 3. Face theta = -90 degrees: frames rotated by thxy = 90, thyz = -90 (the second value overrides the first one
#    written in component 1), then the two in-plane components (local 0 and 1) fixed.
selection.select_pos(pytec3d.model.mesh, 1, (-90.0001, -89.999), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1, mode="u")
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", 90.0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", 90.0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", -90.0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

# 4. Zero potential V at the rod ends (ROD_1, ROD_2): fixes the potential of the solid ring.
selection.select(pytec3d.model.mesh, "ROD_2", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

selection.select(pytec3d.model.mesh, "ROD_1", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 5. Axis: the components normal to the axis vanish.
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 6. Far field: A = 0 on the sphere r = 1 m (index 0 of the spherical system = radius).
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="s", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 7. Load the model into the solver session.
pytec3d.load_model_into_solver()

# 8. Time loop (see the header for the step sequence).
for i in range(1, 10):
    pytec3d.time_step["time_step_number"] = i
    if i == 1:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 0.75  # autots,off (keep as before)
    elif i == 6:
        pytec3d.time_step["time"] = 1e-3
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 1.0  # tintp,1
    elif i == 7:
        pytec3d.time_step["time"] = 2e-3
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 8:
        pytec3d.time_step["time"] = 3e-3
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 0.8  # tintp,0.8
    elif i == 9:
        pytec3d.time_step["time"] = 4e-3
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 0.9  # tintp,0.9
    else:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 0.75  # autots,off (default)

        # select coil elements
    # Source current density: reset Js on the ring, write the toroidal (cylindrical component 1) value i**6, then
    # convert the cylindrical components of every cell to global x, y, z components.
    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    selection_cells = vtknp.vtk_to_numpy(
        pytec3d.model.mesh.GetCellData().GetArray("selection_cell")
    )
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 0, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 1, location="elem")
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", 0, 2, location="elem")
    fact = i**6
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", fact, 1, location="elem")
    model_edit.map_vector_field_to_global_components(pytec3d.model.mesh, "Js",
                                                   coordinate_system=YAxisCylindricalCoordinateSystem(),
                                                   location="elem")
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

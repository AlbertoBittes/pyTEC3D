# =====================================================================================================================
# ITER_coil_3D - 3D 90 degree sector: source ring with Js in an element coordinate system, and solid ring
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The same coil system as ITER_coil_inf_3D, on the other 90 degree sector (z <= 0) and with the source current
# density given through an ELEMENT COORDINATE SYSTEM: the cells of the source ring use the cylindrical system 1, so
# that Js[1] is directly the toroidal component in every cell (vector and tensor parameters are rotated from the
# element system to the global one at the element centre). No conversion of the array is needed.
# The symmetry faces use rotated nodal frames (node_angles), as in ITER_coil_inf_3D.
#
# Geometry and mesh (ITER_coil_3D.vtu):
#   90 degree sector of a sphere of radius 1 m about the y axis (Y-axis cylindrical system: theta from +x towards
#   -z), 10 HEX8, 1276 TET4, 40 PYRAMID5 cells, -1 <= z <= 0 m (theta from 0 to 90). Point arrays: AXIS_NODES (the nodes on the y axis), ROD_1, ROD_2 (end nodes of the conductors).
#   * Source ring (TYPE / MATERIAL 23) at r = 0.25 m, y = +0.25 m; solid ring (TYPE / MATERIAL 24) at y = -0.25 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV3D_nodal / NON_CONDUCTING     / 1 / 0   air
#   * 23: AV3D_nodal / STRANDED_CONDUCTOR / 1 / 0   source ring: Js imposed (element coordinate system 1)
#   * 24: AV3D_nodal / SOLID_CONDUCTOR    / 1 / 0   solid ring: A + V
#   Materials: 21, 23, 24 mur = 1.7 (a permeable medium everywhere); 24 rho = 1e-9 Ohm m.
#   Coordinate systems: 1 = Y-axis cylindrical (added by the macro).
#
# Boundary conditions
#   * Face theta = 0 (plane z = 0), without axis and sphere: frames rotated by thzx = 90 degrees; local components 1
#     and 2 (global y and x, in-plane) fixed to 0.
#   * Face theta = 90 degrees (plane x = 0): frames rotated by (thxy, thyz) = (90, 90); local components 0 and 1
#     (global y and z, in-plane) fixed to 0.
#   * Axis: A_x = A_z = 0.   * Sphere r = 1 m: A = 0.   * V = 0 on ROD_1 / ROD_2.
# Loads and time history (9 steps; Js magnitude i**6 A/m^2 at step i, an ever faster rising current):
#   * step 1, t = 1e-4 s: static;   steps 2-5, t = 2e-4 ... 5e-4 s: transient, theta = 0.75;
#   * step 6, t = 1e-3 s: transient, theta = 1.0;   step 7, t = 2e-3 s: transient, theta = 0.8;
#   * step 8, t = 3e-3 s: STATIC again (eddy currents switched off, theta unused);   step 9, t = 4e-3 s: transient,
#     theta = 0.9.  The # TIMINT / tintp comments of the loop are the equivalent ANSYS commands.
# Results: output/ITER_coil_3D_out_<n>.vtu, ITER_coil_3D_steps.csv, ITER_coil_3D.log
#
# How to run:  cd tests/ITER_coil_3D && python macro.py
# Documentation: "User manual - coordinate systems, boundary conditions", "Element library" (AV3D_nodal).
# =====================================================================================================================


from pytec3d import selection, Solver
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem, SphericalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and model files.
pytec3d = Solver()

filename = "ITER_coil_3D"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Element coordinate system of the source ring: cylindrical about y; Js[1] is then the toroidal component.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")

# 3. Face theta = 0 (without axis and sphere): nodal frames rotated by thzx = 90 degrees, in-plane components fixed.
selection.select_pos(pytec3d.model.mesh, 1, (-0.0001, 0.0001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1, mode="u")
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", 90.0, 2)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 4. Face theta = 90 degrees: frames rotated by thxy = thyz = 90 degrees, in-plane components fixed.
selection.select_pos(pytec3d.model.mesh, 1, (89.999, 90.001), mode="s", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1, mode="u")
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", 90.0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", 90.0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)

# 5. Zero potential V at the rod ends (ROD_1, ROD_2): fixes the potential of the solid ring.
selection.select(pytec3d.model.mesh, "ROD_2", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

selection.select(pytec3d.model.mesh, "ROD_1", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 6. Axis: the components normal to the axis vanish.
selection.select(pytec3d.model.mesh, "AXIS_NODES", 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 7. Far field: A = 0 on the sphere r = 1 m.
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="s", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 8. Load the model; default theta 0.75.
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.75

# 9. Time loop (see the header for the step sequence).
for i in range(1, 10):
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
    # Toroidal source current density i**6 A/m^2 (component 1 in the element cylindrical system).
    selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
    fact = i**6
    selection.set_value_on_selection(pytec3d.model.mesh, "Js", fact, 1, location="elem")

    pytec3d.solve()
    pytec3d.write_output()

# 10. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

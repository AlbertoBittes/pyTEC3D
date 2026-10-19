# =====================================================================================================================
# ITER_coil_3D_period - 3D 90 degree sector with CYCLIC (periodic) conditions instead of symmetry conditions
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The coil system of ITER_coil_inf_3D (source ring + solid ring), where the two sector faces are tied by cyclic
# symmetry instead of being symmetry planes: every field at theta = 0 equals the field at theta = -90 degrees, in
# the rotated frame. It shows:
#   * Model.periodic_boundaries(): the selected nodes of one face are turned by dy = -90 degrees about y, matched with
#     the nodes of the other face, and each pair gets a new CP id in CP_0, coupling A and V;
#   * the nodal frames of the face theta = -90 degrees rotated by thzx = -90 degrees (node_angles), so that the coupled
#     vector components are the same cylindrical components on the two faces;
#   * the grounding of the floating potential V of the solid ring (one node with _V = 0).
#
# Geometry and mesh (ITER_coil_3D_period.vtu): 90 degree sector of a sphere of radius 1 m about the y axis,
#   0 <= z <= 1 m (theta from -90 to 0), 12 HEX8, 534 WEDGE6, 48 TET4, 78 PYRAMID5 cells. Source ring TYPE / MATERIAL 23
#   at r = 0.25 m, y = +0.25 m; solid ring TYPE / MATERIAL 24 at y = -0.25 m.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: AV3D_nodal / NON_CONDUCTING     / 1 / 0   air
#   * 23: AV3D_nodal / STRANDED_CONDUCTOR / 1 / 0   source ring: Js imposed
#   * 24: AV3D_nodal / SOLID_CONDUCTOR    / 1 / 0   solid ring: A + V
#   Materials: mur = 1 everywhere; 24 rho = 1e-9 Ohm m.  DOF couplings: the periodic pairs, added by the macro.
#
# Boundary conditions:  cyclic symmetry theta = 0 <-> -90 degrees (A, V), without the axis and the outer sphere;
#   axis: A_x = A_z = 0; sphere r = 1 m: A = 0; V = 0 on node 107 (on the solid ring).
# Loads and time history (theta = 1): step 1, t = 1e-4 s, static; step 2, t = 2e-4 s, transient. Toroidal Js of
#   magnitude i**6 A/m^2 (cylindrical component 1), converted to global components before each solve.
# Results: output/ITER_coil_3D_period_out_<n>.vtu, ITER_coil_3D_period_steps.csv, ITER_coil_3D_period.log
#
# How to run:  cd tests/ITER_coil_3D_period && python macro.py
# Documentation: "User manual - DOF couplings and periodic boundaries, rotated nodal frames".
# =====================================================================================================================


from pytec3d import mesh, model_edit, selection, Solver
from pytec3d.coordinate_system import SphericalCoordinateSystem, YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and model files.
pytec3d = Solver()

filename = "ITER_coil_3D_period"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# 2. Face theta = -90 degrees, without axis (r = 0) and sphere (r = 1): rotate the nodal frames by thzx = -90 degrees.
selection.select_pos(pytec3d.model.mesh, 1, (-90.0001, -89.999), coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 0, (-0.0001, 0.0001), mode="u", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "node_angles", -90.0, 2)

# 3. Cyclic symmetry: fresh CP_0 array; select the face theta = 0 (exact value 0.0), without axis and sphere, and
#    pair its nodes with their images turned by -90 degrees, coupling A and V.
mesh.remove_array(pytec3d.model.mesh, "CP_0")
mesh.add_array(pytec3d.model.mesh, "CP_0", dtype=vtk.VTK_INT)
selection.select_pos(pytec3d.model.mesh, 1, 0.0, coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh,  0, (-0.0001, 0.0001), mode="u", coordinate_system=YAxisCylindricalCoordinateSystem())
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), mode="u", coordinate_system=SphericalCoordinateSystem())
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(), dy=-90.0, variables=("A", "V"))

# (This reselection does not affect the result: the next select_pos replaces the node selection.)
selection.select(pytec3d.model.mesh, "TYPE", 24, location="elem")
selection.nsle(pytec3d.model.mesh, "r")

# 4. Axis: the components normal to the axis vanish.
selection.select_pos(pytec3d.model.mesh, 0, (-0.0001, 0.0001), coordinate_system=YAxisCylindricalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 5. Ground the potential V of the solid ring on one of its nodes.
selection.select(pytec3d.model.mesh, "ID", 107)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# 6. Far field: A = 0 on the sphere r = 1 m.
selection.select_pos(pytec3d.model.mesh, 0, (0.999, 1.001), coordinate_system=SphericalCoordinateSystem())
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# 7. Load the model into the solver session.
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 0.75

# 8. Time loop: static step, then a transient step (theta = 1).
for i in range(1, 3):
    if i == 1:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 0  # TIMINT,OFF
        pytec3d.time_step["theta"] = 1.0  # autots,off (keep as before)
    else:
        pytec3d.time_step["time"] = i * 0.0001
        pytec3d.time_step["transient_effects"] = 1  # TIMINT,ON
        pytec3d.time_step["theta"] = 1.0  # autots,off (default)

    # select coil elements
    # Toroidal Js (cylindrical component 1) converted to global components.
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

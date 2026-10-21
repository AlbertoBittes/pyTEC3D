# =====================================================================================================================
# coil_3D_avoltemf_edge - coil_3D_avoltemf with EDGE (Nedelec) elements for A and a tree-cotree gauge
# =====================================================================================================================
#
# What this example shows
# -----------------------
# The 30 degree cyclic-sector model of coil_3D_avoltemf (stranded coil + passive ring) with the element type AV3D_edge:
# the vector potential is discretized with first-order Nedelec (edge) elements, whose DOFs are the circulations of A
# along the mesh edges, while V, VOLT and EMF stay nodal (Lagrange). Compared with the nodal macro it shows:
#   * edge DOFs live on the MIDSIDE nodes of a quadratic mesh: their Dirichlet value is component 2 of _A at those
#     nodes ("_A", 0, 2 fixes the tangential A along the selected edges);
#   * no Coulomb gauge term in the weak form: the gauge is a tree-cotree decomposition of the edges, built by
#     tree_gauge.tree_gauge(), which fixes A = 0 on the tree edges (written into _A);
#   * cyclic symmetry with edge DOFs: the coupled edge values are scalars (circulations), so no nodal frame rotation
#     (node_angles) is needed; tree_gauge() is told about the periodic pairs (periodic_point_coupling="CP_1") so that
#     it does not fix the slave edges;
#   * an explicit "free" Dirichlet value: _VOLT = NaN on step 1.
#
# Geometry and mesh (coil_3D_avoltemf_edge.vtu): 30 degree sector about y, 0 <= r <= 1 m, -1 <= y <= 1 m,
#   3800 HEX20 + 200 WEDGE15 (quadratic) cells. Coil TYPE / MATERIAL 23 at r = 0.5 m, y = -0.5 m; passive ring
#   TYPE / MATERIAL 24 at r = 0.5 m, y = +0.5 m. The point array MidsideNodes (1 on midside nodes) is computed when the
#   mesh is read.
#
# Model data (model_data.json)
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 2 = AV3D_edge / NON_CONDUCTING            / 2 / 0   air: edge A
#   * 23: 2 = AV3D_edge / STRANDED_COIL             / 2 / 0   coil: edge A + VOLT + EMF
#   * 24: 2 = AV3D_edge / SOLID_CONDUCTOR           / 2 / 0   passive ring: edge A + V
#   option2 = 2 is the order of the Lagrange fields (V, VOLT, EMF); the edge field is always first order (Nedelec).
#   Materials: 21 air; 23 coil; 24 rho = 2e-8 Ohm m (mur = 1).
#   Real constants 1 (set on every cell by the macro): cross_section = 0.04 m^2, turns = 2000, resistance = 20 Ohm,
#     volume = 0.12566 m^3 (whole coil), symmetry = 12, T = (0, 1, 0) in the cylindrical element system 1.
#   DOF couplings (set by the macro, add_dof_couplings): 1 -> ["VOLT", "EMF"]; periodic pairs (A, V) added by the macro in CP_1.
#
# Boundary conditions:  tangential A = 0 on the far faces r = 1 and y = +-1; tree gauge (tree edges A = 0); cyclic
#   symmetry theta = 0 <-> 30 deg for A and V; V = 0 on the master (lowest id) node of the passive ring.
# Loads and time history (theta = 1, dt = 1e-5 s): step 1 CURR = 100 A (VOLT free); step 2 CURR = 0, VOLT = 10 kV.
# Results: output/coil_3D_avoltemf_edge_out_<n>.vtu, coil_3D_avoltemf_edge_steps.csv, coil_3D_avoltemf_edge.log
#
# How to run:  cd examples/coil_3D_avoltemf_edge && python macro.py
# Documentation: "Element library" (AV3D_edge), "User manual - gauging edge elements", "Theory - 3D edge elements".
# =====================================================================================================================


from pytec3d import mesh, selection, Solver, tree_gauge
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem
import os
import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

# 1. Solver and mesh.
pytec3d = Solver()

filename = "coil_3D_avoltemf_edge"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)
# Every cell uses the real constants row 1 (the CONSTANT cell array selects the row of data["constants"]).
selection.select_all(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CONSTANT", 1, location="elem")
pytec3d.model.read_json("model_data.json")
# DOF couplings (set here, not in the model data): the nodes with CP_0 = 1 share one VOLT and one EMF value.
pytec3d.model.add_dof_couplings({1: ["VOLT", "EMF"]})

# Far boundary, coil size, sector angle (degrees), selection tolerance.
size = 1
coil = 0.2
extr_angle = 30
tol = 1e-5

# Empty point array for the periodic couplings.
mesh.remove_array(pytec3d.model.mesh, "CP_1")
mesh.add_array(pytec3d.model.mesh, "CP_1", dtype=vtk.VTK_INT)

# 2. Element coordinate system of the coil: cylindrical about y, so that T = (0, 1, 0) is toroidal.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
pytec3d.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem()
)
selection.set_value_on_selection(pytec3d.model.mesh, "COORDINATE_SYSTEM", 1, location="elem")


# 3. Coil terminals: VOLT and EMF of the winding coupled (CP_0 = 1); master = lowest node id.
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]

# 4. Cyclic symmetry: the nodes of the face theta = 0 (corner and midside nodes) are paired with their images on
#    theta = 30 degrees; the pairs couple A (the edge values, through the midside nodes) and V.
selection.select_pos(
    pytec3d.model.mesh, 1, (-tol, tol), coordinate_system=YAxisCylindricalCoordinateSystem()
)
pytec3d.model.periodic_boundaries(
            YAxisCylindricalCoordinateSystem(),
            dy=extr_angle,
            array_name="CP_1",
            variables=("A", "V"),
        )

# ---------------------------------------------------------------------------------------------------------------------
# 5. Far field and gauge
#    a) Tangential A = 0 on the faces r = 1, y = +1 and y = -1: component 2 of _A on the selected nodes (only the
#       midside nodes carry edge DOFs; the value on the corner nodes is not used by the edge field).
#    b) tree_gauge(): builds a spanning tree of the edges of the selected cells (all of them: select_all() selected
#       every cell), grown from the selected nodes (the far-field faces), and sets _A = 0 on the midside nodes of the
#       tree edges. The periodic slave edges of CP_1 are left to the periodic constraint.
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(pytec3d.model.mesh)
selection.select_pos(
    pytec3d.model.mesh,
    0,
    (size-tol, size+tol),
    mode="s",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.select_pos(
    pytec3d.model.mesh,
    2,
    (size-tol, size+tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.select_pos(
    pytec3d.model.mesh,
    2,
    (-size-tol, -size+tol),
    mode="a",
    coordinate_system=YAxisCylindricalCoordinateSystem(),
)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)
tree_gauge.tree_gauge(
    pytec3d.model.mesh,
    periodic_point_coupling="CP_1"
)

# 6. Ground the floating potential of the passive ring: V = 0 on its master (lowest id) node.
selection.select(pytec3d.model.mesh, "MATERIAL", 24, location="elem")
selection.nsle(pytec3d.model.mesh)
# Take the first selected node as master node (lowest ID)
master_node_passive = vtknp.vtk_to_numpy(
    pytec3d.model.mesh.GetPointData().GetArray("selection_point")
).nonzero()[0][0]
selection.select(pytec3d.model.mesh, "ID", master_node_passive)
selection.set_value_on_selection(pytec3d.model.mesh, "_V", 0.0)

# 7. Load the model; backward Euler. save_model() here keeps a copy of the prepared model (overwritten at the end).
pytec3d.load_model_into_solver()

pytec3d.time_step["theta"] = 1.0
pytec3d.save_model()

# 8. Time loop: CURR = 100 A with VOLT explicitly free (NaN), then VOLT = 10 kV.
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", np.nan)
    else:
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# 9. Save the solved model (overwrites the input .vtu when run from the case folder).
pytec3d.save_model()

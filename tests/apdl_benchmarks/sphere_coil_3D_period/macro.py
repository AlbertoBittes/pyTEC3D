# =====================================================================================================================
# sphere_coil_3D_period - toroidal coil and passive ring, quarter model, cyclic symmetry
# =====================================================================================================================
#
# Model
# -----
#   Toroidal coil (0.2 <= r <= 0.3, 0.2 <= y <= 0.3) and passive ring (0.2 <= r <= 0.3, -0.3 <= y <= -0.2)
#   around the y axis, inside the sphere r <= 1 (A = 0 on it), 90 degree sector theta in [-90, 0] (z >= 0),
#   extruded mesh (wedges).
#   Js is toroidal: Js[1] in COORDINATE_SYSTEM 1 (cylindrical about y) on the coil cells.
#   Element type: AV3D_nodal.
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * sphere_coil_3D_period.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/sphere_coil_3D_period (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/sphere_coil_3D_period && python macro.py      (outputs in ./output)
# =====================================================================================================================

import os
import sys

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
# tests/apdl_benchmarks: the shared model data, time history and benchmark_case.py (set by the test for a copied case)
BENCHMARKS_DIR = os.environ.get("PYTEC3D_BENCHMARKS_DIR", os.path.dirname(CASE_DIR))
sys.path.insert(0, BENCHMARKS_DIR)

import benchmark_case as bc  # noqa: E402
from pytec3d import selection  # noqa: E402
from pytec3d.coordinate_system import SphericalCoordinateSystem, YAxisCylindricalCoordinateSystem  # noqa: E402

CYL = YAxisCylindricalCoordinateSystem()
SPH = SphericalCoordinateSystem()
TOL = 1e-4

# ---------------------------------------------------------------------------------------------------------------------
# 1. Read the mesh, the unified model data and the time history
# ---------------------------------------------------------------------------------------------------------------------
solver, history = bc.start(CASE_DIR)
mesh = solver.model.mesh

# ---------------------------------------------------------------------------------------------------------------------
# 2. Coil frame: COORDINATE_SYSTEM 1 on the coil cells (Js[1] = toroidal current density)
# ---------------------------------------------------------------------------------------------------------------------
selection.select(mesh, "MATERIAL", bc.ACTIVE_MATERIAL, location="elem")
selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")
selection.select_all(mesh)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Cyclic symmetry between theta = 0 (z = 0) and theta = -90 (x = 0)
#    periodic_boundaries pairs the low-face nodes (axis and sphere excluded) with the high face; the high face
#    gets node_angles[2] = -90 so that the coupled components are the same physical components.
#    Axis: A along the axis; sphere r = 1: A = 0; passive ring: V = 0 on its lowest node.
# ---------------------------------------------------------------------------------------------------------------------
import vtk  # noqa: E402
from pytec3d import mesh as meshes  # noqa: E402

meshes.add_array(mesh, "CP_1", dtype=vtk.VTK_INT)
selection.select_pos(mesh, 1, (-TOL, TOL), coordinate_system=CYL)
selection.select(mesh, "AXIS_NODES", 1, mode="u")
selection.select(mesh, "SPHERE_BOUNDARY", 1, mode="u")
solver.model.periodic_boundaries(CYL, dy=-90.0, array_name="CP_1", variables=("A", "V"))
selection.select_pos(mesh, 1, (-90 - TOL, -90 + TOL), coordinate_system=CYL)
selection.select(mesh, "AXIS_NODES", 1, mode="u")
selection.select(mesh, "SPHERE_BOUNDARY", 1, mode="u")
selection.set_value_on_selection(mesh, "node_angles", -90.0, 2)
selection.select(mesh, "AXIS_NODES", 1)
selection.set_value_on_selection(mesh, "_A", 0, 0)
selection.set_value_on_selection(mesh, "_A", 0, 2)
selection.select(mesh, "SPHERE_BOUNDARY", 1)
for component in range(3):
    selection.set_value_on_selection(mesh, "_A", 0, component)
bc.select_material(solver, bc.PASSIVE_MATERIAL)
selection.select(mesh, "ID", bc.lowest_selected_node(solver))
selection.set_value_on_selection(mesh, "_V", 0)
master = None  # Js-driven: no terminal

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    Js of the time history is applied as Js[1] (toroidal) on the active cells.
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master, js_component=1)

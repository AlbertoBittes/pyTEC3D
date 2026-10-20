# =====================================================================================================================
# coil_3D_avoltemf_timint - coil_3D_avoltemf with time-integrated electric unknowns
# =====================================================================================================================
#
# Model
# -----
#   The coil_2D cross-section revolved about the y axis by 30 degrees in 5 divisions (cyclic symmetry).
#   * Active coil (TYPE 1206, MATERIAL 2) at r = 0.5, y = -0.5: stranded coil (CONSTANT 3: the coil as a 30 degree
#     sector, symmetry 12), current along theta (T = [0, 1, 0] in COORDINATE_SYSTEM 1, cylindrical about y).
#   * Passive ring (TYPE 1205, MATERIAL 3) at r = 0.5, y = +0.5: solid conductor, V = 0 on its lowest node.
#   Element type: AV3D_nodal_timint (Coulomb gauge penalty). V, VOLT, EMF are time integrals [V s].
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * coil_3D_avoltemf_timint.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/coil_3D_avoltemf_timint (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/coil_3D_avoltemf_timint && python macro.py      (outputs in ./output)
# =====================================================================================================================

import os
import sys

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
# tests/apdl_benchmarks: the shared model data, time history and benchmark_case.py (set by the test for a copied case)
BENCHMARKS_DIR = os.environ.get("PYTEC3D_BENCHMARKS_DIR", os.path.dirname(CASE_DIR))
sys.path.insert(0, BENCHMARKS_DIR)

import benchmark_case as bc  # noqa: E402
import vtk  # noqa: E402
from pytec3d import mesh as meshes, selection, tree_gauge  # noqa: E402
from pytec3d.coordinate_system import YAxisCylindricalCoordinateSystem  # noqa: E402

CYL = YAxisCylindricalCoordinateSystem()
SIZE = 1.0                  # cross-section: 0 <= r <= SIZE, -SIZE <= y <= SIZE
SECTOR = 30.0               # sector angle [deg]
TOL = 1e-5

# ---------------------------------------------------------------------------------------------------------------------
# 1. Read the mesh, the unified model data and the time history
# ---------------------------------------------------------------------------------------------------------------------
solver, history = bc.start(CASE_DIR)
mesh = solver.model.mesh

# ---------------------------------------------------------------------------------------------------------------------
# 2. Coil frame: COORDINATE_SYSTEM 1 (cylindrical about y, defined in the model data) on the coil cells
# ---------------------------------------------------------------------------------------------------------------------
selection.select(mesh, "MATERIAL", bc.ACTIVE_MATERIAL, location="elem")
selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")
selection.select_all(mesh)   # the next selections are of nodes; tree_gauge works on the selected cells

# ---------------------------------------------------------------------------------------------------------------------
# 3. Boundary conditions
#    A = 0 on r = SIZE and y = +-SIZE; on the axis A is along the axis (A_x = A_z = 0).
# ---------------------------------------------------------------------------------------------------------------------
selection.select_pos(mesh, 0, (SIZE - TOL, SIZE + TOL), coordinate_system=CYL)
selection.select_pos(mesh, 2, (SIZE - TOL, SIZE + TOL), mode="a", coordinate_system=CYL)
selection.select_pos(mesh, 2, (-SIZE - TOL, -SIZE + TOL), mode="a", coordinate_system=CYL)
for component in range(3):
    selection.set_value_on_selection(mesh, "_A", 0, component)
selection.select_pos(mesh, 0, (-TOL, TOL), coordinate_system=CYL)
selection.set_value_on_selection(mesh, "_A", 0, 0)
selection.set_value_on_selection(mesh, "_A", 0, 2)

# ---------------------------------------------------------------------------------------------------------------------
# 4. Cyclic symmetry: couple A and V of theta = 0 to theta = 30
#    The nodes of the high face get node_angles[2] = 30 (rotated nodal frame), so that the coupled components
#    are the same physical components on both faces.
# ---------------------------------------------------------------------------------------------------------------------
meshes.add_array(mesh, "CP_1", dtype=vtk.VTK_INT)
selection.select_pos(mesh, 1, (-TOL, TOL), coordinate_system=CYL)
solver.model.periodic_boundaries(CYL, dy=SECTOR, array_name="CP_1", variables=("A", "V"))
selection.select_pos(mesh, 1, (SECTOR - TOL, SECTOR + TOL), coordinate_system=CYL)
selection.set_value_on_selection(mesh, "node_angles", SECTOR, 2)

# ---------------------------------------------------------------------------------------------------------------------
# Passive ring: V = 0 on its lowest node (V of a solid conductor is defined up to a constant)
# ---------------------------------------------------------------------------------------------------------------------
bc.select_material(solver, bc.PASSIVE_MATERIAL)
selection.select(mesh, "ID", bc.lowest_selected_node(solver))
selection.set_value_on_selection(mesh, "_V", 0.0)

# ---------------------------------------------------------------------------------------------------------------------
# Coil terminal: VOLT and EMF coupled over the coil (DOF coupling, free CP id); run_history loads its master node
# ---------------------------------------------------------------------------------------------------------------------
master = bc.couple_active(solver, ["VOLT", "EMF"])

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master)

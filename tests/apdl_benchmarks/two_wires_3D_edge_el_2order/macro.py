# =====================================================================================================================
# two_wires_3D_edge_el_2order - 3D: a source bar crossing a conducting bar, edge elements, V of order 2
# =====================================================================================================================
#
# Model
# -----
#   Cube -0.5 <= x, y, z <= 0.5 m, 1210 HEX20. Active bar (x = -0.185, y = 0): Js along z. Passive bar
#   (x = +0.185, y = 0): V = 0 on its two end faces (point arrays LOWER_FACE_ROD, UPPER_FACE_ROD).
#   Element type: AV3D_edge, V of order 2 (TYPE 220x).
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * two_wires_3D_edge_el_2order.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/two_wires_3D_edge_el_2order (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/two_wires_3D_edge_el_2order && python macro.py      (outputs in ./output)
# =====================================================================================================================

import os
import sys

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
# tests/apdl_benchmarks: the shared model data, time history and benchmark_case.py (set by the test for a copied case)
BENCHMARKS_DIR = os.environ.get("PYTEC3D_BENCHMARKS_DIR", os.path.dirname(CASE_DIR))
sys.path.insert(0, BENCHMARKS_DIR)

import benchmark_case as bc  # noqa: E402
from pytec3d import selection, tree_gauge  # noqa: E402

TOL = 1e-5

# ---------------------------------------------------------------------------------------------------------------------
# 1. Read the mesh, the unified model data and the time history
# ---------------------------------------------------------------------------------------------------------------------
solver, history = bc.start(CASE_DIR)
mesh = solver.model.mesh

# ---------------------------------------------------------------------------------------------------------------------
# 2. Boundary condition: A_t = 0 on the six faces (component 2 of _A on the midside nodes)
# ---------------------------------------------------------------------------------------------------------------------
for axis in range(3):
    selection.select_pos(mesh, axis, (0.5 - TOL, 0.5 + TOL))
    selection.select_pos(mesh, axis, (-0.5 - TOL, -0.5 + TOL), mode="a")
    selection.select(mesh, "MidsideNodes", 1, mode="r")
    selection.set_value_on_selection(mesh, "_A", 0, 2)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Passive bar: V = 0 on its end faces
# ---------------------------------------------------------------------------------------------------------------------
for face in ("LOWER_FACE_ROD", "UPPER_FACE_ROD"):
    selection.select(mesh, face, 1)
    selection.set_value_on_selection(mesh, "_V", 0)

# ---------------------------------------------------------------------------------------------------------------------
# 4. Gauge: tree-cotree decomposition (the tree starts on the boundary)
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(mesh)
selection.nsel_ext(mesh)
tree_gauge.tree_gauge(mesh)
master = None  # Js-driven: no terminal

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    Js of the time history is applied along z (component 2) on the active cells.
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master)

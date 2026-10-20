# =====================================================================================================================
# sphere_coil_2D - 2D axisymmetric source coil and passive ring in a half disc, linear elements
# =====================================================================================================================
#
# Model
# -----
#   Axisymmetric about y, half disc r <= 1 m; active coil 0.2 <= x <= 0.3, 0.2 <= y <= 0.3 (Js toroidal, one
#   cell), passive ring 0.2 <= x <= 0.3, -0.3 <= y <= -0.2 (one cell). Linear cells (TRI3, QUAD4).
#   Element type: AV2D_nodal, linear, axisymmetric.
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * sphere_coil_2D.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/sphere_coil_2D (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/sphere_coil_2D && python macro.py      (outputs in ./output)
# =====================================================================================================================

import os
import sys

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
# tests/apdl_benchmarks: the shared model data, time history and benchmark_case.py (set by the test for a copied case)
BENCHMARKS_DIR = os.environ.get("PYTEC3D_BENCHMARKS_DIR", os.path.dirname(CASE_DIR))
sys.path.insert(0, BENCHMARKS_DIR)

import benchmark_case as bc  # noqa: E402

# ---------------------------------------------------------------------------------------------------------------------
# 1. Read the mesh, the unified model data and the time history
# ---------------------------------------------------------------------------------------------------------------------
solver, history = bc.start(CASE_DIR)
mesh = solver.model.mesh

# ---------------------------------------------------------------------------------------------------------------------
# 2. Boundary condition: A_z = 0 on the axis and on r = 1
# ---------------------------------------------------------------------------------------------------------------------
bc.selection.select_all(mesh)
bc.selection.nsel_ext(mesh)
bc.selection.set_value_on_selection(mesh, "_A", 0, 2)
master = None  # Js-driven: no terminal

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    Js of the time history is applied along z (toroidal, component 2) on the active cell.
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master)

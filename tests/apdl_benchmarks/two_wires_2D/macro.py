# =====================================================================================================================
# two_wires_2D - 2D planar: a wire with imposed current density next to a conducting wire
# =====================================================================================================================
#
# Model
# -----
#   Planar, -1 <= x <= 1, -0.5 <= y <= 0.5 m; QUAD8 in the wires, TRI6 in the air.
#   * Active wire (TYPE 203, MATERIAL 2) at x = -0.5: STRANDED_CONDUCTOR, imposed Js along z.
#   * Passive wire (TYPE 205, MATERIAL 3) at x = +0.5: SOLID_CONDUCTOR, Jz = -sigma dA/dt (an infinitely long
#     bar shorted at infinity).
#   Element type: AV2D_nodal, planar.
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * two_wires_2D.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/two_wires_2D (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/two_wires_2D && python macro.py      (outputs in ./output)
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
# 2. Boundary condition: A_z = 0 on the external boundary
# ---------------------------------------------------------------------------------------------------------------------
bc.selection.select_all(mesh)
bc.selection.nsel_ext(mesh)
bc.selection.set_value_on_selection(mesh, "_A", 0, 2)
master = None  # Js-driven: no terminal

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    Js of the time history is applied along z (component 2) on the active cells.
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master)

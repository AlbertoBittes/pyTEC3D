# =====================================================================================================================
# coil_2D_avoltemf_timint - coil_2D_avoltemf with time-integrated coil voltages
# =====================================================================================================================
#
# Model
# -----
#   2D axisymmetric about the y axis (x = r, y = z), 0 <= r <= 1, -1 <= z <= 1 m, 800 QUAD8 of 0.05 m.
#   * Active ring (TYPE 406, MATERIAL 2) at r = 0.5, z = -0.5, 0.2 x 0.2 m: a stranded coil (STRANDED_COIL, CONSTANT 2), driven by its current CURR or its
#     terminal voltage VOLT.
#   * Passive ring (TYPE 405, MATERIAL 3) at r = 0.5, z = +0.5: a shorted ring (SOLID_CONDUCTOR) carrying the
#     induced currents.
#   Element type: AV2D_nodal_timint  (V, VOLT, EMF are time integrals [V s], also V / VOLT of the time history).
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * coil_2D_avoltemf_timint.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/coil_2D_avoltemf_timint (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/coil_2D_avoltemf_timint && python macro.py      (outputs in ./output)
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
# 2. Boundary condition: A_z = 0 on the external boundary (axis and far boundary)
# ---------------------------------------------------------------------------------------------------------------------
bc.selection.select_all(mesh)
bc.selection.nsel_ext(mesh)
bc.selection.set_value_on_selection(mesh, "_A", 0, 2)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Terminal of the active ring
#    couple_active() couples its VOLT and EMF over all its nodes (a DOF coupling with a free CP id): one terminal, loaded on the
#    master node (the lowest node id) by run_history. The passive ring needs nothing (a 2D solid conductor is a
#    closed, shorted ring).
# ---------------------------------------------------------------------------------------------------------------------
master = bc.couple_active(solver, ['VOLT', 'EMF'])

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master)

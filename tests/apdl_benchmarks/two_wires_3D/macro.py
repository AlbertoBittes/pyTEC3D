# =====================================================================================================================
# two_wires_3D - 3D: a wire with imposed current density next to a conducting wire, nodal elements
# =====================================================================================================================
#
# Model
# -----
#   Box -1 <= x <= 1, -0.5 <= y, z <= 0.5 m, two 0.1 x 0.1 m wires along z at x = -0.2 (active, Js along z)
#   and x = +0.2 (passive). Hexahedra in the wires, free tetrahedra and pyramids in the air.
#   Element type: AV3D_nodal.
#
# Data
# ----
#   * tests/apdl_benchmarks/model_data.json: the model data shared by all the benchmark cases (element types, materials, real
#     constants, B-H curve, coordinate system), see tests/apdl_benchmarks/benchmark_case.py for its ids:
#     TYPE = 100 * formulation + role, MATERIAL 1 air, 2 active, 3 passive, 4 iron.
#   * tests/apdl_benchmarks/time_history.json: the steps (time, theta, transient_effects) and the loads of the active coil (Js,
#     AMPS / V, CURR / VOLT; the formulation of the coil decides which are used), the same for every benchmark case.
#   * two_wires_3D.vtu: the mesh, written from the MAPDL mesh of pyTEC3D_APDL_benchmarks/cases/two_wires_3D (same nodes).
#
# Benchmark: tests/test_apdl_benchmarks.py compares the outputs with MAPDL (pyTEC3D_APDL_benchmarks).
# How to run:  cd tests/apdl_benchmarks/two_wires_3D && python macro.py      (outputs in ./output)
# =====================================================================================================================

import os
import sys

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
# tests/apdl_benchmarks: the shared model data, time history and benchmark_case.py (set by the test for a copied case)
BENCHMARKS_DIR = os.environ.get("PYTEC3D_BENCHMARKS_DIR", os.path.dirname(CASE_DIR))
sys.path.insert(0, BENCHMARKS_DIR)

import benchmark_case as bc  # noqa: E402
from pytec3d import selection  # noqa: E402

TOL = 1e-5

# ---------------------------------------------------------------------------------------------------------------------
# 1. Read the mesh, the unified model data and the time history
# ---------------------------------------------------------------------------------------------------------------------
solver, history = bc.start(CASE_DIR)
mesh = solver.model.mesh

# ---------------------------------------------------------------------------------------------------------------------
# 2. Boundary conditions
#    A = 0 on x = +-1 and y = +-0.5. The wires cross the faces z = +-0.5: the current is normal to them and B
#    tangent, so the tangential components of A vanish there (A_x = A_y = 0).
# ---------------------------------------------------------------------------------------------------------------------
selection.select_pos(mesh, 0, (1 - TOL, 1 + TOL))
selection.select_pos(mesh, 0, (-1 - TOL, -1 + TOL), mode="a")
selection.select_pos(mesh, 1, (0.5 - TOL, 0.5 + TOL), mode="a")
selection.select_pos(mesh, 1, (-0.5 - TOL, -0.5 + TOL), mode="a")
for component in range(3):
    selection.set_value_on_selection(mesh, "_A", 0, component)
selection.select_pos(mesh, 2, (0.5 - TOL, 0.5 + TOL))
selection.select_pos(mesh, 2, (-0.5 - TOL, -0.5 + TOL), mode="a")
selection.set_value_on_selection(mesh, "_A", 0, 0)
selection.set_value_on_selection(mesh, "_A", 0, 1)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Passive wire: V = 0 on its two end faces (z = +-0.5)
# ---------------------------------------------------------------------------------------------------------------------
bc.select_material(solver, bc.PASSIVE_MATERIAL)
selection.select_pos(mesh, 2, (0.5 - TOL, 0.5 + TOL), mode="r")
selection.set_value_on_selection(mesh, "_V", 0)
bc.select_material(solver, bc.PASSIVE_MATERIAL)
selection.select_pos(mesh, 2, (-0.5 - TOL, -0.5 + TOL), mode="r")
selection.set_value_on_selection(mesh, "_V", 0)
master = None  # Js-driven: no terminal

# ---------------------------------------------------------------------------------------------------------------------
# Load the model and run the time history
#    Js of the time history is applied along z (component 2) on the active cells.
#    run_history: for every step, time / theta / transient_effects and the loads of the active coil, then
#    solve() and write_output() (output/<case>_out_<n>.vtu).
# ---------------------------------------------------------------------------------------------------------------------
solver.load_model_into_solver()
bc.run_history(solver, history, master)

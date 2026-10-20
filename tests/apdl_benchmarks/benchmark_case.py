"""Shared setup of the benchmark cases (the pyTEC3D twins of pyTEC3D_APDL_benchmarks).

Every benchmark case has the same structure: air, an ACTIVE coil or conductor (MATERIAL 2) that carries the excitation,
a PASSIVE conductor (MATERIAL 3) in which eddy currents are induced and, in the *_ferr cases, iron (MATERIAL 4). The
data shared by all of them is in two files:

* ``tests/apdl_benchmarks/model_data.json``: one model data for every case. TYPE id = 100 * formulation + role:

  ======  =================================================  ======  =====================================
  code    formulation                                        role    domain (option1)
  ======  =================================================  ======  =====================================
  1       AV2D_nodal, quadratic, axisymmetric                1       air: NON_CONDUCTING
  2       AV2D_nodal, quadratic, planar                      2       iron: SOFT_MAGNETIC
  3       AV2D_nodal, linear, axisymmetric                   3       active: STRANDED_CONDUCTOR (Js)
  4       AV2D_nodal_timint, quadratic, axisymmetric         4       active: CIRCUIT_COUPLED_CONDUCTOR (2D)
  11      AV3D_nodal, linear                                 5       passive: SOLID_CONDUCTOR
  12      AV3D_nodal_timint, linear                          6       active: STRANDED_COIL
  21      AV3D_edge, V of order 1
  22      AV3D_edge, V of order 2
  23      AV3D_edge_timint, V of order 2
  ======  =================================================  ======  =====================================

  MATERIAL 1 air, 2 active, 3 passive, 4 iron (B-H curve): the same values in every case. CONSTANT 1 generic
  (thickness 1); 2 and 3 the same coil (cross-section, turns, resistance, volume), whole (symmetry 1, current out of
  plane: 2D) or as a 30 degree sector (symmetry 12, current along theta of COORDINATE_SYSTEM 1: 3D).
  COORDINATE_SYSTEM 1: cylindrical about Y. No DOF couplings: each case couples what it needs (couple_active).

* ``tests/apdl_benchmarks/time_history.json``: the steps of every case, each with ``time``, ``theta``, ``transient_effects``
  and the loads of the active coil: ``Js`` [A/m2], one of ``AMPS`` / ``V`` and one of ``CURR`` / ``VOLT`` (null = not given). The
  domain of the active cells decides which are used: STRANDED_CONDUCTOR -> Js; CIRCUIT_COUPLED_CONDUCTOR or
  SOLID_CONDUCTOR -> AMPS (current) or V (voltage); STRANDED_COIL -> CURR or VOLT. In the ``_timint`` element types
  V and VOLT are time integrals [V s].

The MAPDL macros of pyTEC3D_APDL_benchmarks read the same two files, so changing a value there changes both models.
Other files can replace them: ``$PYTEC3D_BENCHMARK_MODEL_DATA`` and ``$PYTEC3D_BENCHMARK_TIME_HISTORY`` (used by
``bench.py assess`` of pyTEC3D_APDL_benchmarks).

Usage in a case ``macro.py``::

    solver, history = benchmark_case.start(CASE_DIR)
    ...  # boundary conditions of the case
    master = benchmark_case.couple_active(solver, ["VOLT", "EMF"])   # or ["V"]; None for a Js-driven coil
    solver.load_model_into_solver()
    benchmark_case.run_history(solver, history, master)
"""

from __future__ import annotations

import json
import os

import numpy as np
import vtk.util.numpy_support as vtknp

from pytec3d import Solver, selection

from benchmark_inputs import BENCHMARKS_DIR, MODEL_DATA, TIME_HISTORY  # noqa: E402,F401  (the inputs of the run)

ACTIVE_MATERIAL = 2
PASSIVE_MATERIAL = 3
IRON_MATERIAL = 4

# domain (option1) -> load family of the active coil
_KIND = {3: "Js", 4: "AMPS/V", 5: "AMPS/V", 6: "CURR/VOLT"}


def start(case_dir: str) -> tuple[Solver, list[dict]]:
    """Read ``<case_dir>/<case>.vtu``, the unified model data and the unified time history."""
    name = os.path.basename(os.path.normpath(case_dir))
    solver = Solver()
    solver.read_model(os.path.join(case_dir, name + ".vtu"))
    solver.model.read_json(MODEL_DATA)
    with open(TIME_HISTORY) as f:
        history = json.load(f)["steps"]
    return solver, history


def _cell_ids(mesh, material: int) -> np.ndarray:
    values = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("MATERIAL"))
    return np.nonzero(values == material)[0]


def active_kind(solver: Solver) -> str:
    """'Js', 'AMPS/V' or 'CURR/VOLT', from the domain of the TYPE of the active cells."""
    mesh = solver.model.mesh
    types = vtknp.vtk_to_numpy(mesh.GetCellData().GetArray("TYPE"))[_cell_ids(mesh, ACTIVE_MATERIAL)]
    kinds = {_KIND[solver.model.data["element_types"][str(int(t))]["option1"]] for t in np.unique(types)}
    if len(kinds) != 1:
        raise ValueError(f"the active cells (MATERIAL {ACTIVE_MATERIAL}) must have one load family, got {kinds}")
    return kinds.pop()


def select_material(solver: Solver, material: int) -> None:
    """Select the cells of a material and their nodes."""
    selection.select(solver.model.mesh, "MATERIAL", material, location="elem")
    selection.nsle(solver.model.mesh)


def lowest_selected_node(solver: Solver) -> int:
    flags = vtknp.vtk_to_numpy(solver.model.mesh.GetPointData().GetArray("selection_point"))
    return int(flags.nonzero()[0][0])


def couple_active(solver: Solver, fields: list[str]) -> int:
    """Couple ``fields`` over all the nodes of the active coil (one terminal) and return its master node.

    The coupling gets the next free CP id (``periodic_boundaries`` numbers its own pairs from 1) and is added to the
    model data here, by the case: the unified model data holds no DOF couplings, since a coupling of a field the model
    does not have is an error. Circuit-coupled conductor: ``["V"]``; stranded coil: ``["VOLT", "EMF"]``. The master
    node (lowest node id) carries the loads of the terminal.
    """
    cp = max(solver.model._used_coupling_ids(), default=0) + 1
    solver.model.add_dof_couplings({cp: list(fields)})
    select_material(solver, ACTIVE_MATERIAL)
    selection.set_value_on_selection(solver.model.mesh, "CP_0", cp)
    return lowest_selected_node(solver)


def _set_node(solver: Solver, node: int, array: str, value: float) -> None:
    selection.select(solver.model.mesh, "ID", node)
    selection.set_value_on_selection(solver.model.mesh, array, value)


def apply_step_loads(solver: Solver, step: dict, master: int | None, js_component: int = 2) -> None:
    """Loads of one step on the active coil, as its formulation needs them."""
    kind = active_kind(solver)
    mesh = solver.model.mesh
    if kind == "Js":
        if step.get("Js") is None:
            raise ValueError(f"step at t = {step['time']}: Js is needed by a stranded conductor")
        selection.select(mesh, "MATERIAL", ACTIVE_MATERIAL, location="elem")
        for component in range(3):
            value = step["Js"] if component == js_component else 0.0
            selection.set_value_on_selection(mesh, "Js", value, component, location="elem")
        return
    current, voltage = ("AMPS", "V") if kind == "AMPS/V" else ("CURR", "VOLT")
    if step.get(current) is not None:
        _set_node(solver, master, current, step[current])
        _set_node(solver, master, "_" + voltage, np.nan)
    elif step.get(voltage) is not None:
        _set_node(solver, master, current, 0.0)
        _set_node(solver, master, "_" + voltage, step[voltage])
    else:
        raise ValueError(f"step at t = {step['time']}: one of {current} / {voltage} is needed by the active coil")


def run_history(solver: Solver, history: list[dict], master: int | None, js_component: int = 2) -> None:
    """Solve the steps of the time history, writing output/<case>_out_<n>.vtu after each."""
    for step in history:
        solver.time_step["time"] = step["time"]
        solver.time_step["theta"] = step["theta"]
        solver.time_step["transient_effects"] = bool(step["transient_effects"])
        apply_step_loads(solver, step, master, js_component)
        selection.select_all(solver.model.mesh)
        solver.solve()
        solver.write_output()

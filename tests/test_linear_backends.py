"""The linear solvers of the MFEM backend (``solver_settings["linear"]["backend"] = "mfem"``) against PETSc.

Patch tests (the meshes and exact fields of the formulations suite, pyTEC3D_APDL_benchmarks/suites/formulations/
common/patch.py): every solver must reproduce the uniform field of each formulation, the direct ones to round-off, the
iterative ones to their tolerance. The stranded ball of the suite with first-order edges: CG with AMS without any gauge
(the singular curl-curl system, its right-hand side made consistent) against MUMPS with the tree gauge.

Run: ``pytest tests/test_linear_backends.py`` (the benchmark repository as for test_benchmarks.py).
"""

import json
import os
import sys

import numpy as np
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCHMARKS = os.path.abspath(os.environ.get("PYTEC3D_BENCHMARKS")
                             or os.path.join(PROJECT_ROOT, "pyTEC3D_APDL_benchmarks"))
SUITE = os.path.join(BENCHMARKS, "suites", "formulations")
if not os.path.isdir(SUITE):
    pytest.skip(f"benchmark repository not found at {BENCHMARKS}", allow_module_level=True)
sys.path.insert(0, os.path.join(SUITE, "common"))

import patch  # noqa: E402
import physics  # noqa: E402
import pytec3d_driver  # noqa: E402

from pytec3d.solver import DEFAULT_SOLVER_SETTINGS  # noqa: E402
from pytec3d.validation import ValidationError, check_solver_settings  # noqa: E402

ITERATIVE = {"relative_tolerance": 1e-12, "absolute_tolerance": 1e-30, "max_iterations": 1000}
SOLVERS = {
    "petsc_lu": {},
    "mfem_lu": {"backend": "mfem"},
    "mfem_cholesky": {"backend": "mfem", "pc_type": "cholesky"},
    "mfem_cg_boomeramg": {"backend": "mfem", "ksp_type": "cg", "pc_type": "boomeramg", **ITERATIVE},
    "mfem_gmres_jacobi": {"backend": "mfem", "ksp_type": "gmres", "pc_type": "jacobi", **ITERATIVE},
    "mfem_cg_ams": {"backend": "mfem", "ksp_type": "cg", "pc_type": "ams", **ITERATIVE},
}
#: The patch variants of each solver (AMS: one first-order Nedelec field only).
VARIANTS = {
    "petsc_lu": ["PHI3D_nodal_tet_perturbed", "AVA3D_edge_tet_o2_perturbed", "AVA2D_nodal_tri_o2_quadratic"],
    "mfem_lu": ["PHI3D_nodal_tet_perturbed", "TPHI3D_edge_pyramid", "AVA3D_edge_tet_o2_perturbed",
                "AVA3D_nodal_wedge_o2_perturbed", "AVA2D_nodal_tri_o2_quadratic", "TPHI2D_edge_quad_axisym"],
    "mfem_cholesky": ["PHI3D_nodal_hex", "AVA3D_nodal_tet_o2_quadratic"],
    "mfem_cg_boomeramg": ["PHI3D_nodal_tet_perturbed", "AVA3D_nodal_hex_o2", "AVA2D_nodal_quad_o2_axisym"],
    "mfem_gmres_jacobi": ["TPHI3D_edge_tet", "AVA2D_nodal_tri_o1_perturbed"],
    "mfem_cg_ams": ["AVA3D_edge_tet_o1_perturbed", "AVA3D_edge_hex_o1", "AVA3D_edge_wedge_o1", "AVA3D_edge_pyramid_o1"],
}


@pytest.fixture
def linear(monkeypatch, tmp_path):
    """Sets the linear solver settings of the drivers; the runs in a temporary folder."""
    monkeypatch.chdir(tmp_path)

    def use(settings):
        monkeypatch.setattr(pytec3d_driver, "QUIET", {**pytec3d_driver.QUIET, "linear": settings})
    return use


@pytest.mark.parametrize("solver,variant", [(s, v) for s, variants in VARIANTS.items() for v in variants])
def test_patch(solver, variant, linear):
    linear(SOLVERS[solver])
    out = pytec3d_driver.patch_test(variant)
    exact = patch.exact_cells(variant)
    tol = 1e-10 if solver.endswith(("_lu", "_cholesky")) else 1e-8
    for key in ("B_cells", "H_cells"):
        error = np.linalg.norm(out[key] - exact[key]) / np.linalg.norm(exact[key])
        assert error < tol, f"{solver} {variant} {key}: {error:.2e}"


def test_ams_without_gauge(linear, monkeypatch):
    """The stranded ball, first-order AVA3D_edge: MUMPS with the tree gauge against CG + AMS without gauge."""
    with open(os.path.join(SUITE, "inputs", "standard", "settings.json")) as f:
        settings = json.load(f)
    monkeypatch.setattr(pytec3d_driver, "FIRST_ORDER", pytec3d_driver.FIRST_ORDER + ("AVA3D_edge",))
    linear(SOLVERS["mfem_lu"])
    gauged = pytec3d_driver.stranded_sphere(settings, "AVA3D_edge")
    linear({**SOLVERS["mfem_cg_ams"], "relative_tolerance": 1e-10})
    monkeypatch.setattr(pytec3d_driver.tree_gauge, "tree_gauge", lambda *args, **kwargs: None)
    free = pytec3d_driver.stranded_sphere(settings, "AVA3D_edge")
    exact = physics.rotating_b_mean(settings["stranded_sphere"]["j0"], settings["mesh3d"]["radius"])
    # The same B (the gauge only changes A); the ungauged source has no gradient part, a small change at this mesh.
    assert np.linalg.norm(free["B_cells"] - gauged["B_cells"]) / np.linalg.norm(gauged["B_cells"]) < 0.02
    assert abs(free["B_mean"][2] - exact) / exact < 0.02  # first order, the coarse mesh of the suite: 1.3 %


@pytest.mark.parametrize("settings,message", [
    ({"backend": "trilinos"}, "linear.backend"),
    ({"backend": "mfem", "pc_type": "ilu"}, "pc_type"),
    ({"backend": "mfem", "ksp_type": "lgmres"}, "ksp_type"),
])
def test_invalid_settings(settings, message):
    merged = {**DEFAULT_SOLVER_SETTINGS, "linear": {**DEFAULT_SOLVER_SETTINGS["linear"], **settings}}
    with pytest.raises(ValidationError, match=message):
        check_solver_settings(merged)

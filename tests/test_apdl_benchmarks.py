"""Benchmark of the pyTEC3D test cases against ANSYS MAPDL.

pyTEC3D does not need ANSYS: the MAPDL reference results live in the repository pyTEC3D_APDL_benchmarks (APDL macros,
the tools that run them and convert the results, and reference sets of ``<case>/<case>_ref_<set>.vtu``). A reference
set holds the results of one pair of inputs (model data, time history), identified by their fingerprint
(``reference_set.json``, see ``tests/apdl_benchmarks/benchmark_inputs.py``).

Reference set used (the first found):

* ``$PYTEC3D_APDL_REFERENCES``: a reference set folder, as is (e.g. one made by ``bench.py run`` with other inputs);
* the set ``reference/`` of the revision pinned in ``tests/apdl_benchmarks/references.json``, fetched into a temporary
  folder from ``$PYTEC3D_APDL_BENCHMARKS``, else from the repository named there (a path relative to the pyTEC3D root,
  or a URL). Updating the references of pyTEC3D = committing them there and pinning the new revision here.

Inputs: ``tests/apdl_benchmarks/model_data.json`` and ``time_history.json``, or ``$PYTEC3D_BENCHMARK_MODEL_DATA`` and
``$PYTEC3D_BENCHMARK_TIME_HISTORY``. The tests are skipped, with the reason, when no reference set is found or when it
was computed with other inputs (different fingerprint).

For every ``<case>/manifest.json`` of the set:

1. the pyTEC3D twin ``tests/apdl_benchmarks/<pytec3d_case>`` is copied to a temporary folder and its ``macro.py`` is
   run there (the macros overwrite their input mesh with ``save_model()``, the copy keeps ``tests/`` clean);
2. each output ``output/<case>_out_<n>.vtu`` is paired with the reference set ``n`` (one MAPDL substep per pyTEC3D
   step);
3. nodes are matched by coordinates and cells by their corner nodes (the meshes are the same, MAPDL numbering
   included; matching geometrically keeps the test independent of the numbering and of degenerate MAPDL shapes);
4. for each field of the manifest the relative L2 error ``||pyTEC3D - MAPDL|| / ||MAPDL||`` over the matched nodes or
   cells is computed for every set and compared with the field tolerance ``rtol``. Fields with ``"assert": false`` are
   reported only (known formulation differences, explained in their ``why``).

Run: ``pytest tests/test_apdl_benchmarks.py -q`` (summary table at the end, ``-s`` also prints every field and set;
``PYTEC3D_APDL_REPORT=<file.json>`` writes the table as JSON and Markdown, see ``tests/conftest.py``).
"""

import glob
import json
import os
import atexit
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest
import vtk
from vtk.util.numpy_support import vtk_to_numpy

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(TESTS_DIR, "apdl_benchmarks")   # the pyTEC3D twins, model_data.json, time_history.json
PROJECT_ROOT = os.path.dirname(TESTS_DIR)

sys.path.insert(0, CASES_DIR)
import benchmark_inputs  # noqa: E402

# Number of corner nodes of the VTK cell types (linear and quadratic)
CORNERS = {5: 3, 9: 4, 10: 4, 12: 8, 13: 6, 14: 5, 22: 3, 23: 4, 24: 4, 25: 8, 26: 6, 27: 5, 28: 4}


def fetch_pinned_references():
    """(reference set folder, description): the set of the revision pinned in references.json, in a temporary folder."""
    with open(os.path.join(CASES_DIR, "references.json")) as f:
        pin = json.load(f)
    repo = os.environ.get("PYTEC3D_APDL_BENCHMARKS") or pin["repository"]
    if "://" not in repo and not repo.startswith("git@"):
        repo = os.path.abspath(os.path.join(PROJECT_ROOT, repo))
        if not os.path.isdir(repo):
            return None, f"repository {repo} not found (set $PYTEC3D_APDL_BENCHMARKS or $PYTEC3D_APDL_REFERENCES)"
    tmp = tempfile.mkdtemp(prefix="pytec3d_apdl_refs_")
    atexit.register(shutil.rmtree, tmp, ignore_errors=True)
    try:
        subprocess.run(["git", "clone", "--quiet", "--no-checkout", repo, tmp], check=True, capture_output=True)
        subprocess.run(["git", "-C", tmp, "checkout", "--quiet", pin["revision"], "--", pin["reference_set"]],
                       check=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = (getattr(exc, "stderr", None) or b"").decode(errors="replace").strip()
        return None, f"cannot fetch {repo} at {pin['revision']}: {detail or exc}"
    return os.path.join(tmp, pin["reference_set"]), f"{repo} @ {pin['revision'][:12]} ({pin['reference_set']})"


def find_references():
    """(reference set folder, its description, reason to skip: None when usable)."""
    folder = os.environ.get("PYTEC3D_APDL_REFERENCES")
    source = os.path.abspath(folder) if folder else None
    if not folder:
        folder, source = fetch_pinned_references()
        if folder is None:
            return None, None, source
    info_path = os.path.join(folder, "reference_set.json")
    if not os.path.isfile(info_path):
        return None, source, f"{info_path} not found: not a reference set"
    with open(info_path) as f:
        info = json.load(f)
    fingerprint = benchmark_inputs.file_fingerprint()
    if info["inputs_fingerprint"] != fingerprint:
        return None, source, (
            f"the references in {source} were computed with inputs {info['inputs_fingerprint']}, the tests use "
            f"{fingerprint} ({benchmark_inputs.MODEL_DATA}, {benchmark_inputs.TIME_HISTORY}): compute references for "
            "these inputs with 'bench.py run --model-data ... --time-history ...' in pyTEC3D_APDL_benchmarks")
    return folder, source, None


def load_manifests(folder):
    if folder is None:
        return []
    with open(os.path.join(folder, "reference_set.json")) as f:
        cases = json.load(f)["cases"]
    return [json.load(open(os.path.join(folder, c, "manifest.json"))) for c in cases]


def read_vtu(path):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    return reader.GetOutput()


def array(grid, name, location):
    data = grid.GetPointData() if location == "point" else grid.GetCellData()
    arr = data.GetArray(name)
    if arr is None:
        raise KeyError(f"array '{name}' ({location}) not found")
    values = vtk_to_numpy(arr).astype(float)
    return values.reshape(len(values), -1)


def match_points(ref, out, rtol=1e-9):
    """Index of the output point at the location of every reference point (closest point within rtol * size)."""
    ref_pts = vtk_to_numpy(ref.GetPoints().GetData())
    out_pts = vtk_to_numpy(out.GetPoints().GetData())
    tol = rtol * max(np.ptp(out_pts, axis=0).max(), 1e-30)
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(out)
    locator.BuildLocator()
    index = np.array([locator.FindClosestPoint(p) for p in ref_pts])
    far = np.linalg.norm(out_pts[index] - ref_pts, axis=1) > tol
    if np.any(far):
        raise AssertionError(f"{np.sum(far)} reference nodes have no pyTEC3D node at the same location")
    return index


def cell_corner_keys(grid, point_map=None):
    keys = []
    for c in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(c)
        ids = [cell.GetPointId(i) for i in range(CORNERS[grid.GetCellType(c)])]
        if point_map is not None:
            ids = [int(point_map[i]) for i in ids]
        keys.append(frozenset(ids))
    return keys


def match_cells(ref, out, point_map):
    lookup = {k: i for i, k in enumerate(cell_corner_keys(out))}
    index = np.array([lookup.get(k, -1) for k in cell_corner_keys(ref, point_map)])
    if np.any(index < 0):
        raise AssertionError(f"{np.sum(index < 0)} reference cells have no pyTEC3D cell with the same corners")
    return index


def material_nodes(grid, material):
    """Mask of the points of the cells of a MATERIAL."""
    mat = vtk_to_numpy(grid.GetCellData().GetArray("MATERIAL"))
    keep = np.zeros(grid.GetNumberOfPoints(), dtype=bool)
    for c in np.nonzero(mat == material)[0]:
        cell = grid.GetCell(int(c))
        keep[[cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]] = True
    return keep


def relative_error(p, r):
    norm = np.linalg.norm(r)
    if norm == 0.0:
        return 0.0 if np.linalg.norm(p) == 0.0 else np.inf
    return float(np.linalg.norm(p - r) / norm)


def run_pytec3d(case_dir, workdir):
    shutil.copytree(case_dir, workdir, ignore=shutil.ignore_patterns("output"))
    # the copied macro finds benchmark_case.py there, and uses the inputs of the test
    env = dict(os.environ, PYTEC3D_BENCHMARKS_DIR=CASES_DIR, PYTEC3D_BENCHMARK_MODEL_DATA=benchmark_inputs.MODEL_DATA,
               PYTEC3D_BENCHMARK_TIME_HISTORY=benchmark_inputs.TIME_HISTORY)
    proc = subprocess.run([sys.executable, "macro.py"], cwd=workdir, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-40:])
        raise RuntimeError(f"pyTEC3D macro failed (exit {proc.returncode}):\n{tail}")


REFERENCES, REFERENCES_SOURCE, SKIP_REASON = find_references()
manifests = load_manifests(REFERENCES)

#: One entry per (case, field), filled by the tests and printed by conftest.pytest_terminal_summary.
SUMMARY: list[dict] = []
#: The reference set and the inputs of the run, for the report.
RUN_INFO = {"references": REFERENCES_SOURCE, "inputs_fingerprint": benchmark_inputs.file_fingerprint(),
            "model_data": benchmark_inputs.MODEL_DATA, "time_history": benchmark_inputs.TIME_HISTORY,
            "solver": manifests[0].get("solver") if manifests else None, "skipped": SKIP_REASON}


@pytest.mark.skipif(not manifests, reason=f"APDL benchmarks skipped: {SKIP_REASON or 'empty reference set'}")
@pytest.mark.parametrize("manifest", manifests, ids=lambda m: m["case"])
def test_apdl_benchmark(manifest, tmp_path):
    case = manifest["case"]
    pytec3d_case = manifest["pytec3d_case"]
    workdir = str(tmp_path / pytec3d_case)
    try:
        run_pytec3d(os.path.join(CASES_DIR, pytec3d_case), workdir)
    except Exception as exc:   # counted as failing in the summary, with the last line of the error
        last = (str(exc).strip().splitlines() or ["?"])[-1]
        SUMMARY.append({"case": case, "field": "(pyTEC3D run)", "rtol": float("nan"), "err": float("nan"), "set": None,
                        "info_err": None, "passed": False, "error": last})
        raise

    rows, failures = [], []
    worst = {}   # field name -> summary entry
    for s in manifest["sets"]:
        ref = read_vtu(os.path.join(REFERENCES, case, s["file"]))
        out_path = os.path.join(workdir, "output", f"{pytec3d_case}_out_{s['set']}.vtu")
        assert os.path.isfile(out_path), f"pyTEC3D output {out_path} not found"
        out = read_vtu(out_path)
        point_map = match_points(ref, out)
        cell_map = None
        for f in manifest["fields"]:
            if f.get("sets") and s["set"] not in f["sets"]:
                continue
            if f["location"] == "cell" and cell_map is None:
                cell_map = match_cells(ref, out, point_map)
            index = point_map if f["location"] == "point" else cell_map
            r = array(ref, f["reference"], f["location"])[:, f["components"]]
            p = array(out, f["pytec3d"], f["location"])[index][:, f["components"]]
            mask = f.get("mask")
            if mask == "reference_nonzero":
                keep = np.any(r != 0.0, axis=1)
                r, p = r[keep], p[keep]
            elif isinstance(mask, dict) and "material" in mask:   # nodes of the cells of a material
                keep = material_nodes(ref, mask["material"])
                r, p = r[keep], p[keep]
            p = np.nan_to_num(p)
            err = relative_error(p, r)
            checked = f.get("assert", True)
            ok = err <= f["rtol"]
            entry = worst.setdefault(f["name"], {"case": case, "field": f["name"], "rtol": f["rtol"], "err": None,
                                                 "set": None, "info_err": None, "passed": True})
            if checked:
                if entry["err"] is None or err >= entry["err"]:
                    entry["err"], entry["set"] = err, s["set"]
                entry["passed"] = entry["passed"] and ok
            else:
                entry["info_err"] = max(err, entry["info_err"] or 0.0)
            rows.append(f"  set {s['set']:2d}  t={s['time']:<10.4g} {f['name']:14s} rel.L2 err {err:10.3e}  "
                        f"(rtol {f['rtol']:.0e}) {'ok' if ok else ('FAIL' if checked else 'info')}")
            if checked and not ok:
                failures.append(f"set {s['set']} {f['name']}: {err:.3e} > {f['rtol']:.0e}")

    SUMMARY.extend(worst.values())
    print(f"\n[APDL benchmark] {case} vs MAPDL ({manifest.get('solver', '?')})\n" + "\n".join(rows))
    assert not failures, f"{case}: " + "; ".join(failures)

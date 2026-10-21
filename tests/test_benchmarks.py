"""Validation of the installed pyTEC3D: runs the benchmark repository (pyTEC3D_APDL_benchmarks).

The benchmark cases, their inputs, the pyTEC3D drivers and the stored ANSYS MAPDL / APEC reference results all live
in the benchmark repository (the single source of truth); this test only invokes its validation,
``bench.py validate``, with the Python of this process (so the pytec3d it imports is the one installed here). No
ANSYS installation is needed.

Benchmark repository: ``$PYTEC3D_BENCHMARKS`` (a path), else ``pyTEC3D_APDL_benchmarks`` in the pyTEC3D tree (the
git submodule). Options, passed on to the validation:

* ``PYTEC3D_BENCHMARKS_LONG=1``: also the long cases (3D ITER cases, about an hour more);
* ``PYTEC3D_BENCHMARKS_REPORT=<folder>``: write the reports (simple.json/.md, iter.json/.md) there;
* ``PYTEC3D_BENCHMARKS_ARGS``: more arguments of ``bench.py validate`` (e.g. ``"--suite simple -k coil_2D"``).

Run: ``pytest tests/test_benchmarks.py -s`` (``-s`` shows the summary tables of the validation as they are printed).
"""

import os
import shlex
import subprocess
import sys

import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BENCHMARKS = os.path.abspath(os.environ.get("PYTEC3D_BENCHMARKS")
                             or os.path.join(PROJECT_ROOT, "pyTEC3D_APDL_benchmarks"))


def test_benchmarks():
    bench = os.path.join(BENCHMARKS, "bench.py")
    if not os.path.isfile(bench):
        pytest.fail(f"benchmark repository not found at {BENCHMARKS}: run 'git submodule update --init' or set "
                    "$PYTEC3D_BENCHMARKS to a clone of pyTEC3D_APDL_benchmarks")
    command = [sys.executable, bench, "validate"]
    if os.environ.get("PYTEC3D_BENCHMARKS_LONG", "") not in ("", "0"):
        command.append("--long")
    if os.environ.get("PYTEC3D_BENCHMARKS_REPORT"):
        command += ["--report", os.path.abspath(os.environ["PYTEC3D_BENCHMARKS_REPORT"])]
    command += shlex.split(os.environ.get("PYTEC3D_BENCHMARKS_ARGS", ""))
    print("$ " + " ".join(command), flush=True)
    code = subprocess.call(command, cwd=BENCHMARKS)
    assert code == 0, f"validation against the benchmark references failed (exit {code}); see the summary above"

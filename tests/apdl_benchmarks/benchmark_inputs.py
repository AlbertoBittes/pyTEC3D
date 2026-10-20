"""Inputs of the benchmark cases: which model data and time history they use, and their fingerprint.

No pyTEC3D import: tests/test_apdl_benchmarks.py uses it in the pytest process, benchmark_case.py in the case macros.

The defaults are ``model_data.json`` and ``time_history.json`` of this folder; ``$PYTEC3D_BENCHMARK_MODEL_DATA`` and
``$PYTEC3D_BENCHMARK_TIME_HISTORY`` replace them (assessment of other inputs against MAPDL results computed with the
same files: ``bench.py assess`` of pyTEC3D_APDL_benchmarks).
"""

import hashlib
import json
import os

BENCHMARKS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MODEL_DATA = os.path.join(BENCHMARKS_DIR, "model_data.json")
DEFAULT_TIME_HISTORY = os.path.join(BENCHMARKS_DIR, "time_history.json")
MODEL_DATA = os.path.abspath(os.environ.get("PYTEC3D_BENCHMARK_MODEL_DATA") or DEFAULT_MODEL_DATA)
TIME_HISTORY = os.path.abspath(os.environ.get("PYTEC3D_BENCHMARK_TIME_HISTORY") or DEFAULT_TIME_HISTORY)


def inputs_fingerprint(model_data, time_history):
    """Fingerprint of a pair of inputs (dicts): sha256 of the model data and of the steps of the time history.

    Canonical JSON (sorted keys, no spaces), so formatting and descriptions do not count. The same function is in
    pyTEC3D_APDL_benchmarks/tools/inputs.py: every MAPDL reference set records the fingerprint of its inputs.
    """
    payload = json.dumps({"model_data": model_data, "time_history": time_history["steps"]}, sort_keys=True,
                         separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def file_fingerprint(model_data_path=MODEL_DATA, time_history_path=TIME_HISTORY):
    """Fingerprint of the inputs in two files (by default the inputs of this run)."""
    with open(model_data_path) as f:
        model_data = json.load(f)
    with open(time_history_path) as f:
        time_history = json.load(f)
    return inputs_fingerprint(model_data, time_history)

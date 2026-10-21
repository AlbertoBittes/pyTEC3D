"""Pytest configuration: run the tests without writing bytecode into the source tree.

The benchmark tables and reports are produced by the benchmark repository (pyTEC3D_APDL_benchmarks, validation/), which
tests/test_benchmarks.py invokes.
"""

import sys

sys.dont_write_bytecode = True

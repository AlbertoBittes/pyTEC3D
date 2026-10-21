Testing and Benchmarks
======================

pyTEC3D is validated against ANSYS Mechanical APDL. Its field formulations are compared with MAPDL on 23 analyses:
2D and 3D, nodal and edge elements on hexahedra, prisms and pyramids, stranded coils, massive and circuit-coupled
conductors, time-integrated formulations, iron with a B-H curve. Its plasma algorithm (:mod:`pytec3d.algorithm`) is
compared with APEC, the APDL implementation of the same algorithm, on ITER equilibria and vertical displacements in
2D and 3D, with halo currents and stranded coils.

Everything about the benchmarks lives in a separate repository, **pyTEC3D_APDL_benchmarks** (a git submodule of
pyTEC3D, folder ``pyTEC3D_APDL_benchmarks``):

* the cases, their inputs and the pyTEC3D drivers;
* the MAPDL / APEC models and macros;
* the stored reference results;
* the comparison, the parametric studies and the reports.

pyTEC3D only invokes it. **pyTEC3D does not need ANSYS to be tested**: the MAPDL and APEC results of the default
inputs are stored in that repository.

Quick start
-----------

.. code-block:: bash

   git submodule update --init                 # the benchmark repository (once)
   pytest tests/test_benchmarks.py -s          # simple suite + 2D ITER cases, about 10 minutes

``build.sh`` runs the same validation at the end of a build (``--validate=yes``, the default; ``--validate=long`` also
runs the 3D ITER cases, ``--validate=no`` skips it). The test calls ``bench.py validate`` of the benchmark repository
with the Python of the test, so the pytec3d it imports is the one installed in this environment. Options:

.. code-block:: bash

   PYTEC3D_BENCHMARKS_LONG=1 pytest tests/test_benchmarks.py -s            # also the 3D ITER cases (about 1 h more)
   PYTEC3D_BENCHMARKS_REPORT=out pytest tests/test_benchmarks.py -s        # JSON / Markdown reports in out/
   PYTEC3D_BENCHMARKS_ARGS="--suite simple -k coil_2D" pytest tests/test_benchmarks.py -s
   PYTEC3D_BENCHMARKS=/path/to/pyTEC3D_APDL_benchmarks pytest tests/test_benchmarks.py -s   # another clone

or directly, in the benchmark repository: ``python bench.py validate [--long] [--suite simple|iter] [-k EXPR]``.

The summary
-----------

The validation prints one table per suite.

* **Simple suite:** per case and field, the worst relative L2 error ``||pyTEC3D - MAPDL|| / ||MAPDL||`` over the
  result sets, its tolerance and the status.
* **ITER suite:** per case and quantity, the worst difference over the steps, for example plasma current, area,
  current centroid, li, axis and boundary flux, passive, coil and halo currents, chamber flux, fields A and B. Each has
  its tolerance, and is either asserted or only reported.

The validation **fails** when the stored references do not match the current default inputs of the benchmark
repository, rather than skipping: a validation that cannot run is not a passed validation.

Parametric studies
------------------

The same repository runs the benchmarks with other inputs and writes the data of a sensitivity study:

.. code-block:: bash

   cd pyTEC3D_APDL_benchmarks
   python bench.py simple study --param model_data.materials.3.rho --values 1e-7,2e-7,4e-7 --no-mapdl
   python bench.py iter study ITER_2D_evol --param algorithm.halo_delta1 --values 0.02,0.04,0.08 --no-apec

With ``--no-mapdl`` / ``--no-apec`` only pyTEC3D runs (no ANSYS); without them MAPDL / APEC run too and the two codes
are compared for every value.

Each study writes:

* the input files of every value, so every point can be reproduced;
* ``study.md``, a summary;
* ``study_long.csv``, one row per value, case, code, step and quantity, ready for pandas.

The README of the benchmark repository describes the inputs (one model data and one time history for the simple suite,
one ``parameters.json`` for the ITER suite), the reference sets and how the default references are updated.

Assumptions
-----------

The comparison rests on choices documented in ``ASSUMPTIONS.md`` of the benchmark repository. The main ones:

* **Same formulations.** For example, the stranded coil of :ref:`theory-stranded-coil` [MAPDLTheory2024]_: real coil
  volume ``volume``, with the unmodelled part accounted for by ``symmetry``.
* **Element values.** pyTEC3D element values are cell means (weighted by the quadrature JxW); the reference element
  values are centroid values. They are close where the fields are resolved on the mesh.
* **Quadrature rules differ.** Where the fields are under-resolved on the mesh, the solutions differ by more than the
  default tolerances.
* **Constant theta.** Theta must be constant over the time history: MAPDL does not apply a change of theta between
  load steps correctly. For theta < 1 the two codes evaluate the eddy current density Jt differently: pyTEC3D uses the
  velocity of the theta scheme, MAPDL the mean rate over the step.
* **Edge elements.** Their gauge-dependent quantities (A, V) are not compared.
* **ITER suite.** The one intended difference is the identification of the plasma boundary: APEC uses the poloidal
  field B, pyTEC3D the flux PSI.

What is still open is in ``OPEN_POINTS.md`` of the same repository.

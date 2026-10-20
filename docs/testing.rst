Testing and Benchmarks
======================

pyTEC3D is tested by benchmarking it against ANSYS Mechanical APDL: 23 analyses (2D and 3D, nodal and edge elements
on hexahedra, prisms and pyramids, stranded coils, massive and circuit-coupled conductors, time-integrated
formulations, iron with a B-H curve) are run by pyTEC3D and compared, field by field and time step by time step,
with the same analyses solved by MAPDL.

**pyTEC3D does not need ANSYS to run the tests.** The MAPDL results are stored in a separate repository,
``pyTEC3D_APDL_benchmarks``, together with the APDL macros and the tools that produce them; the tests fetch them.

Quick start
-----------

After the build (:doc:`installation`), in the pyTEC3D environment:

.. code-block:: bash

   python -m pip install -e ".[dev]"              # pytest (once)
   pytest tests/test_apdl_benchmarks.py -q         # about 2 minutes; summary table at the end

The test needs ``git`` and the benchmarks repository: by default ``pyTEC3D_APDL_benchmarks`` inside the pyTEC3D tree
(clone it there), or the path or URL in ``$PYTEC3D_APDL_BENCHMARKS``.

Other ways to run it:

.. code-block:: bash

   pytest tests/test_apdl_benchmarks.py -q -s                 # also the error of every field and result set
   pytest tests/test_apdl_benchmarks.py -q -k coil_2D         # selected cases (substring of the case names)
   pytest tests/test_apdl_benchmarks.py -q -x                 # stop at the first failing case
   PYTEC3D_APDL_REPORT=report.json pytest tests/test_apdl_benchmarks.py -q   # also report.json and report.md

The summary
-----------

At the end of the run pytest prints one line per case and compared field:

.. code-block:: text

   ========= APDL benchmark summary (relative L2 error pyTEC3D vs MAPDL) ==========
   references .../pyTEC3D_APDL_benchmarks @ 88dfc2a1428c (reference), inputs 275d7a804f6a07f6
   case                           field                 err     rtol  set       info  status
   coil_2D_avoltemf               A_z              6.30e-05    1e-03    4             ok
                                  B                4.10e-03    2e-02    1             ok
                                  Jt               8.14e-06    2e-02    4             ok
                                  VOLT terminal    8.22e-06    1e-03    1             ok
                                  EMF terminal     1.23e-05    1e-03    1             ok
   ...
   21 / 21 cases within tolerance

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Column
     - Meaning
   * - ``err``
     - Worst relative L2 error :math:`\|u_\text{pyTEC3D}-u_\text{MAPDL}\|/\|u_\text{MAPDL}\|` over the result sets
       (time steps) where the field is checked, over the nodes (``A``, terminal ``VOLT`` / ``EMF``) or the cells
       (``B``, ``Jt``) of the model.
   * - ``rtol``
     - Tolerance of the field for that case (``case.json`` in the benchmarks repository).
   * - ``set``
     - Result set (time step) of the worst error.
   * - ``info``
     - Worst error on sets where the field is reported but not checked (a field known to differ by formulation).
   * - ``status``
     - ``ok``; ``FAIL`` (an error above ``rtol``, or the pyTEC3D run itself failed: the last line of its error is
       shown); ``info`` for a field that is reported but never checked (``err`` is ``-``, the error is in ``info``),
       e.g. when the two models are not like-for-like (the reason is in the case's ``case.json``).

The report written with ``PYTEC3D_APDL_REPORT`` holds the same table (Markdown and JSON) with the reference set, the
inputs and the pyTEC3D revision, for a record or a merge request.

What a case does
----------------

The cases are in ``tests/apdl_benchmarks/<case>/`` (mesh, ``macro.py``). For each of them the test:

1. copies the case to a temporary folder and runs its ``macro.py`` (``tests/`` stays clean);
2. pairs each output ``output/<case>_out_<n>.vtu`` with the MAPDL result set ``n`` (one MAPDL substep per pyTEC3D
   step);
3. matches nodes by position and cells by their corner nodes (the meshes are the MAPDL ones);
4. computes the relative L2 error of each field and compares it with its tolerance.

All the cases share **one model data and one time history**, ``tests/apdl_benchmarks/model_data.json`` and
``time_history.json``, read by the macros through ``benchmark_case.py``. Every case has air (MATERIAL 1), an active
coil or conductor (2), a passive conductor (3) and, in the ``_ferr`` cases, iron (4); the TYPE id of a cell is
``100 * formulation + role``. The time history has 4 steps: current drive, voltage drive, zero. See
:ref:`examples-benchmarks` for the families of cases.

Reference results and inputs
----------------------------

A **reference set** is a folder of MAPDL results for one pair of inputs (model data, time history), identified by a
*fingerprint*, a hash of the content of the two files (descriptions and formatting do not count).

* The set of the default inputs is ``reference/`` in the benchmarks repository. pyTEC3D pins the revision it is
  tested against in ``tests/apdl_benchmarks/references.json``; the test clones that revision into a temporary folder
  and deletes it afterwards.
* The test is **skipped, with the reason**, when the repository cannot be found or fetched, or when the reference
  set was computed with inputs other than those of the run (fingerprints differ). A skipped run is not a pass: read
  the line printed under ``APDL benchmarks skipped``.

.. list-table:: Environment variables
   :header-rows: 1
   :widths: 40 60

   * - Variable
     - Effect
   * - ``PYTEC3D_APDL_BENCHMARKS``
     - Path or URL of pyTEC3D_APDL_benchmarks (default: as in ``references.json``, ``./pyTEC3D_APDL_benchmarks``).
   * - ``PYTEC3D_APDL_REFERENCES``
     - A reference set folder, used as is instead of the pinned revision (e.g. references not yet committed).
   * - ``PYTEC3D_BENCHMARK_MODEL_DATA``, ``PYTEC3D_BENCHMARK_TIME_HISTORY``
     - Other inputs instead of ``model_data.json``, ``time_history.json`` (they must match the reference set).
   * - ``PYTEC3D_APDL_REPORT``
     - Writes the summary to ``<file>.json`` and ``<file>.md``.

Benchmarking other inputs (needs MAPDL)
---------------------------------------

To see how pyTEC3D compares with MAPDL for another model data or time history (other materials, a shorter time step,
...), edit copies of the two files and, in the benchmarks repository, with the Python of the pyTEC3D environment:

.. code-block:: bash

   cd pyTEC3D_APDL_benchmarks
   python bench.py run --model-data my_model_data.json --time-history my_time_history.json   # MAPDL, ~2 minutes
   python bench.py assess results/<fingerprint>      # runs the pyTEC3D tests against it
   python bench.py status                            # all reference sets and their last assessment

``run`` writes the MAPDL results to ``results/<fingerprint>/`` (the default inputs go to ``reference/``); ``assess``
runs ``tests/test_apdl_benchmarks.py`` with those inputs and references and writes ``assessment.md`` and
``assessment.json`` in the same folder. Neither pyTEC3D nor the committed references change.

Updating the references
-----------------------

When the default inputs, a case macro or a mesh change, the MAPDL results are recomputed and pinned again:

.. code-block:: bash

   cd pyTEC3D_APDL_benchmarks
   python bench.py run [cases] [--meshes]    # --meshes: also rewrite the pyTEC3D meshes from MAPDL
   python bench.py assess                    # check pyTEC3D against them
   git commit -am "..."                      # commit reference/ (and the cases)
   python bench.py pin                       # write the revision into pyTEC3D's references.json

then commit ``tests/apdl_benchmarks/references.json`` in pyTEC3D together with the inputs or meshes that changed.

Assumptions
-----------

The comparison rests on choices documented in ``ASSUMPTIONS.md`` of the benchmarks repository; the main ones:

* The formulations are the same in both codes: e.g. the stranded coil of :ref:`theory-stranded-coil`
  [MAPDLTheory2024]_ (real coil volume ``volume``, the unmodelled part accounted for by ``symmetry``).
* pyTEC3D element values are cell means (weighted by the quadrature JxW), the reference element values are centroid
  values: close where the fields are resolved on the mesh, as with the default inputs.
* Quadrature rules differ: with fields under-resolved on the mesh (short time steps, coarse cells) the solutions
  differ by more than the default tolerances.
* Theta is constant over the time history (MAPDL does not apply a change of theta between load steps correctly). The
  default history uses theta = 0.5, which checks every term of the theta scheme; backward Euler (theta = 1) is the
  variant ``variants/time_history_theta1.json`` of the benchmarks repository.
* Gauge-dependent quantities of edge elements (A, V) are not compared.

What is still open is in ``OPEN_POINTS.md`` of the same repository.

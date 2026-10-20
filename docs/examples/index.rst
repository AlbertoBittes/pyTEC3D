Examples
========

The folder ``tests/`` holds two kinds of complete analyses:

* ``tests/examples/``: the examples of this section, each in its own folder with the mesh (``<case>.vtu``), the model
  data (``model_data.json``) and the script (``macro.py``). Every macro starts with a summary (physics, geometry,
  element types, materials, real constants, boundary conditions, loads, results) and comments each block of commands;
  the pages of this section explain the families of cases and show the macros in full.
* ``tests/apdl_benchmarks/``: the cases benchmarked against ANSYS MAPDL (:ref:`examples-benchmarks`).

The :doc:`tutorial notebooks <../tutorials/index>` build and run analyses step by step; their own inputs are in
``docs/tutorials/data``.

Running an example
------------------

The macros use paths relative to the working directory: run them from their folder.

.. code-block:: bash

   cd tests/examples/ITER_2D
   python macro.py                      # or: mpirun -np 4 python macro.py

They write ``<case>.log``, ``output/`` and ``emat/`` in the folder, and most of them end with ``save_model()``, which
overwrites the input ``<case>.vtu`` with the solved model. To keep the repository clean, copy the folder first (with
``tests/examples/files`` and, for an evolution, its equilibrium case, as below).

The ITER evolution cases read the output of an equilibrium case (``../ITER_2D/output/...``) and the waveforms of
``tests/examples/files``: run the equilibrium first, and keep the folder structure.

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Case
     - Element type
     - What it shows
   * - :doc:`conductors`: ``two_wires_3D_2order``, ``two_wires_3D_edge``
     - ``AV3D_nodal``, ``AV3D_edge``
     - Source wire with imposed Js and a passive wire in a box, second-order nodal and edge elements.
   * - :doc:`conductors`: ``3D_loop_cond_base``
     - ``AV3D_nodal``
     - Quarter model of two square loops: named boundary arrays, symmetry planes, potential on symmetry cuts.
   * - :doc:`create`: ``CREATE_B_measurement(_nodal)``, ``CREATE_force_measurement(_nodal)``
     - ``AV3D_edge``, ``AV3D_nodal``
     - A laboratory electromagnet: racetrack stranded coil with element coordinate systems, saturating conducting
       yoke, voltage waveform.
   * - :doc:`iter`: ``ITER_2D*``, ``ITER_2D_evol*``, ``ITER_3D_evol*``
     - all
     - Plasma equilibrium and vertical displacement of ITER with the plasma algorithm, in 2D and on 3D toroidal sectors.

Note: the evolution cases read the most recent ``<equilibrium case>_out_*.vtu`` of the equilibrium case: the
equilibrium writes one output, numbered by the solves of its Picard loop.

.. _examples-benchmarks:

Benchmarks against ANSYS MAPDL
------------------------------

The cases of ``tests/apdl_benchmarks/`` (2D and 3D coils, conductors and rings: nodal, edge on hexahedra, prisms
and pyramids, time-integrated, with iron) share one model data and one time history (``model_data.json``, ``time_history.json`` in that folder): TYPE ids
select formulation and role, materials and real constants are the same for all, and the time history drives the active
coil with Js, a current or a voltage, as its formulation requires (``benchmark_case.py``). Their macros differ from the
examples above: they read the shared data through ``benchmark_case.py`` and their meshes are written from the MAPDL
models.

``pytest tests/test_apdl_benchmarks.py`` runs them and compares their fields with MAPDL results: how to run the tests,
read the summary, benchmark other inputs and update the references is in :doc:`../testing`.

.. toctree::
   :maxdepth: 2

   conductors
   create
   iter

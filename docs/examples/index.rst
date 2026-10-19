Examples
========

The folder ``tests/`` holds complete analyses, each in its own folder with the mesh (``<case>.vtu``), the model data
(``model_data.json``) and the script (``macro.py``). Every macro starts with a summary (physics, geometry, element
types, materials, real constants, boundary conditions, loads, results) and comments each block of commands; the
pages of this section explain the families of cases and show the macros in full.

Running an example
------------------

The macros use paths relative to the working directory: run them from their folder.

.. code-block:: bash

   cd tests/coil_2D_avolt
   python macro.py                      # or: mpirun -np 4 python macro.py

They write ``<case>.log``, ``output/`` and ``emat/`` in the folder, and most of them end with ``save_model()``, which
overwrites the input ``<case>.vtu`` with the solved model. To keep the repository clean, copy the folder first:

.. code-block:: bash

   cp -r tests/coil_2D_avolt /tmp/ && cd /tmp/coil_2D_avolt && python macro.py

The ITER evolution cases read the output of an equilibrium case (``../ITER_2D/output/...``) and the waveforms of
``tests/files``: run the equilibrium first, and keep the folder structure.

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - Case
     - Element type
     - What it shows
   * - :doc:`coils_2d`: ``coil_2D_avolt``, ``_ferr``, ``_timint``
     - ``AV2D_nodal(_timint)``
     - 2D axisymmetric circuit-coupled conductor driven by current then voltage, passive solid ring; iron with B-H
       curve; time-integrated potential.
   * - :doc:`coils_2d`: ``coil_2D_avoltemf``, ``_ferr``, ``_timint``
     - ``AV2D_nodal(_timint)``
     - Stranded coil with circuit equation (VOLT, EMF), passive solid ring.
   * - :doc:`coils_3d`: ``coil_3D_avoltemf``, ``_ferr``, ``_timint``
     - ``AV3D_nodal(_timint)``
     - The stranded coil in 3D on a 30 degree cyclic sector: periodic couplings, rotated nodal frames, element
       coordinate system for the winding direction.
   * - :doc:`coils_3d`: ``coil_3D_avoltemf_edge``, ``_ferr``, ``_timint``, ``_prisms``, ``_pyramids``
     - ``AV3D_edge``
     - The same with edge elements and a tree-cotree gauge, on hexahedra, prisms and pyramids.
   * - :doc:`conductors`: ``two_wires_2D``, ``two_wires_3D``, ``_2order``, ``_edge``, ``_edge_el``, ``_edge_el_2order``
     - all
     - Source wire with imposed Js and a passive wire, 2D and 3D, nodal and edge, first and second order.
   * - :doc:`conductors`: ``3D_loop_cond_base``
     - ``AV3D_nodal``
     - Quarter model of two square loops: named boundary arrays, symmetry planes, potential on symmetry cuts.
   * - :doc:`conductors`: ``ITER_coil_inf_2D``, ``ITER_coil_inf_3D``, ``ITER_coil_3D``, ``_test_1``, ``_period``
     - ``AV2D_nodal``, ``AV3D_nodal``
     - Rings on 90 degree sectors: static / transient step sequences, rotated frames, Js in cylindrical components,
       cyclic symmetry.
   * - :doc:`create`: ``CREATE_B_measurement(_nodal)``, ``CREATE_force_measurement(_nodal)``
     - ``AV3D_edge``, ``AV3D_nodal``
     - A laboratory electromagnet: racetrack stranded coil with element coordinate systems, saturating conducting
       yoke, voltage waveform.
   * - :doc:`iter`: ``ITER_2D*``, ``ITER_2D_evol*``, ``ITER_3D_evol*``
     - all
     - Plasma equilibrium and vertical displacement of ITER with the plasma algorithm, in 2D and on 3D toroidal sectors.

Notes on the test inputs
------------------------

* ``two_wires_3D``, ``_2order``, ``_edge``: the planes ``x = +-0.5`` selected for the far field are interior to the
  mesh (its x faces are at +-1 m), so that selection constrains (almost) nothing; the macros say so in a comment.
* The evolution cases read the most recent ``<equilibrium case>_out_*.vtu`` of the equilibrium case: the equilibrium
  writes one output, numbered by the solves of its Picard loop.

.. toctree::
   :maxdepth: 2

   coils_2d
   coils_3d
   conductors
   create
   iter

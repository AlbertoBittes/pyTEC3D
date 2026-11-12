Examples
========

The folder ``examples/`` of the repository holds complete analyses, each in its own folder:

* the mesh (``<case>.vtu``);
* the model data (``model_data.json``);
* the script (``macro.py``).

Every macro starts with a summary (physics, geometry, element types, materials, real constants, boundary conditions,
loads, results) and comments each block of commands. The pages of this section explain the families of cases and show
the macros in full. ``examples/files`` holds the waveforms they read.

The :doc:`tutorial notebooks <../tutorials/index>` run some of these examples step by step, with pre- and
post-processing.

The benchmarks against ANSYS MAPDL and APEC are not examples: they live in the benchmark repository
(``pyTEC3D_APDL_benchmarks``, :doc:`../testing`).

Running an example
------------------

The macros use paths relative to the working directory: run them from their folder.

.. code-block:: bash

   cd examples/coil_2D_avolt
   python macro.py                      # or: mpirun -np 4 python macro.py

They write ``<case>.log``, ``output/`` and ``emat/`` in the folder, and most of them end with ``save_model()``, which
overwrites the input ``<case>.vtu`` with the solved model. To keep the repository clean, copy the folder first, together
with ``examples/files`` and, for an evolution, its equilibrium case, keeping the folder structure. The ITER evolution
cases read the output of their equilibrium case (``../ITER_2D_fine2/output/...``) and the waveforms of
``../files``: run the equilibrium first.

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 32 20 48

   * - Case
     - Element type
     - What it shows
   * - :doc:`coils_2d`: ``coil_2D_avolt``, ``coil_2D_avoltemf``, ``coil_2D_avoltemf_ferr``, ``coil_2D_avolt_timint``
     - ``AVA2D_nodal``, ``AVA2D_nodal_timint``
     - Axisymmetric coils and rings: circuit-coupled conductor, stranded coil with its circuit equation, saturating
       iron, time-integrated potential; current- and voltage-driven steps.
   * - :doc:`coils_3d`: ``coil_3D_avoltemf``, ``coil_3D_avoltemf_edge``, ``coil_3D_avoltemf_edge_ferr``,
       ``coil_3D_avoltemf_timint``
     - ``AVA3D_nodal``, ``AVA3D_edge``, ``AVA3D_nodal_timint``
     - The same rings on a 30 degree cyclic sector: coordinate systems, cyclic symmetry of nodal and edge fields, tree
       gauge, iron with edge elements, time-integrated unknowns.
   * - :doc:`iter`: ``ITER_2D_carma``, ``ITER_2D_fine2``, ``ITER_2D_evol_carma``, ``ITER_2D_evol_fine2``,
       ``ITER_3D_evol_fine2``, ``ITER_3D_evol_coils_edge``
     - all
     - Plasma equilibria and vertical displacement of ITER with the plasma algorithm: linear and quadratic meshes, 3D
       sectors with nodal and edge elements, halo currents, voltage-driven coils.
   * - :doc:`formulations`: ``core_3D_phi``, ``ring_3D_tphi``, ``coil_core_3D_tphi``, ``wire_2D_tphi``,
       ``bulk_2D_tphi_axisym``
     - ``PHI3D_nodal``, ``TPHI3D_edge``, ``TPHI2D_edge``
     - The scalar potentials: magnetostatics with the source field in the system (iron loop, B-H curve), T-phi eddy
       currents with cuts and fixed-T coils, superconductors with power laws (AC transport loss, a field-cooled bulk
       with Jc(B)), 3D, 2D and axisymmetric.

.. toctree::
   :maxdepth: 2

   coils_2d
   coils_3d
   iter
   formulations

2D Axisymmetric Coils and Rings
===============================

Four transient eddy-current analyses of the same geometry: two rings of square cross-section, 0.2 m x 0.2 m, at
r = 0.5 m, z = -0.5 m (active) and z = +0.5 m (passive), in a 1 m x 2 m axisymmetric domain of second-order cells
(``AV2D_nodal``, x = r, y = z). The active ring is driven through its terminal, by a current or by a voltage. The passive
ring carries the induced eddy currents. Each variant adds one feature:

.. list-table::
   :header-rows: 1
   :widths: 26 22 52

   * - Case
     - Active ring
     - What it shows
   * - ``coil_2D_avolt``
     - massive conductor with a terminal (``CIRCUIT_COUPLED_CONDUCTOR``)
     - Reading a mesh and its model data, a boundary condition from a selection, a DOF coupling (one terminal V), a
       current-driven step (``AMPS``) then a voltage-driven step (``_V``), theta method, output files.
   * - ``coil_2D_avoltemf``
     - multi-turn ``STRANDED_COIL``
     - The circuit equation of a coil, VOLT = R I + EMF (VOLT and EMF coupled on the coil nodes), the real constants of
       a coil (cross section, turns, radius, resistance, winding direction), current (``CURR``) and voltage drive.
   * - ``coil_2D_avoltemf_ferr``
     - stranded coil, with an iron frame
     - A saturating B-H curve (``bh_curve`` law, Newton-Raphson at every step), backward Euler, the operating range of
       the law.
   * - ``coil_2D_avolt_timint``
     - massive conductor, ``AV2D_nodal_timint``
     - The time-integrated electric potential (symmetric A-v system): the terminal DOF is the time integral of the
       voltage, AMPS is still the current.

Each folder holds the mesh (``<case>.vtu``), its model data (``model_data.json``) and the commented script
(``macro.py``). The same geometries, with one shared model data, are benchmarked against ANSYS MAPDL in the benchmark
repository (:doc:`../testing`). :doc:`../tutorials/tutorial_coil_2D` and
:doc:`../tutorials/tutorial_stranded_coil_iron` run the first and third case step by step, with plots.

coil_2D_avolt
-------------

.. literalinclude:: ../../examples/coil_2D_avolt/macro.py
   :language: python

coil_2D_avoltemf
----------------

.. literalinclude:: ../../examples/coil_2D_avoltemf/macro.py
   :language: python

coil_2D_avoltemf_ferr
---------------------

.. literalinclude:: ../../examples/coil_2D_avoltemf_ferr/macro.py
   :language: python

coil_2D_avolt_timint
--------------------

.. literalinclude:: ../../examples/coil_2D_avolt_timint/macro.py
   :language: python

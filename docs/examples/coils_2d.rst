2D Axisymmetric Coils
=====================

Six cases on the same mesh: two rings of 0.2 m x 0.2 m cross-section at r = 0.5 m, z = -0.5 m (lower) and
z = +0.5 m (upper), in an air domain 0 <= r <= 1 m, -1 <= z <= 1 m (800 QUAD8 cells, second-order elements,
axisymmetric about the y axis). They are the reference examples of the terminal formulations of pyTEC3D.

.. list-table::
   :header-rows: 1
   :widths: 26 20 54

   * - Case
     - Lower ring (active)
     - Variation
   * - ``coil_2D_avolt``
     - ``CIRCUIT_COUPLED_CONDUCTOR``
     - Reference: current-driven step (``AMPS``), then voltage-driven step (``_V``); upper ring
       ``SOLID_CONDUCTOR`` (a closed ring: induced current only).
   * - ``coil_2D_avolt_ferr``
     - ``CIRCUIT_COUPLED_CONDUCTOR``
     - Iron frame (``HARD_MAGNETIC`` with M = 0) with a B-H curve; 30 kA; Newton iterations.
   * - ``coil_2D_avolt_timint``
     - ``CIRCUIT_COUPLED_CONDUCTOR``
     - ``AV2D_nodal_timint``: V is the time integral of the voltage.
   * - ``coil_2D_avoltemf``
     - ``STRANDED_COIL``
     - 1000-turn coil with its circuit equation (VOLT, EMF); upper ring ``SOLID_CONDUCTOR``.
   * - ``coil_2D_avoltemf_ferr``
     - ``STRANDED_COIL``
     - With the iron frame; backward Euler.
   * - ``coil_2D_avoltemf_timint``
     - ``STRANDED_COIL``
     - Time-integrated VOLT and EMF.

Physics
-------

The lower ring is driven; the changing flux induces in the upper ring a current opposing the change (Lenz). With a
**circuit-coupled conductor** the current distribution in the driven ring itself is computed (skin effect); its
terminal, V coupled over the cross-section, imposes either the total current (``AMPS``) or the voltage across the
circumference (:ref:`theory-solid-2d`). The passive ring is a **solid conductor**: no terminal, only the induced
current of a closed ring. With a **stranded coil** the current density is uniform, and the coil current follows the
circuit equation ``VOLT = R i + dPsi/dt`` (:ref:`theory-stranded-coil`). The current-driven first step at 100 A is a
step from zero at t = 0, so the eddy currents of step 1 are the response to that step.

Setting up the coil
-------------------

The key lines of the stranded-coil version:

.. code-block:: python

   selection.select(mesh, "MATERIAL", 23, location="elem")      # the coil cells
   selection.nsle(mesh)
   selection.set_value_on_selection(mesh, "CP_0", 1)             # dof_couplings["1"] = ["EMF", "VOLT"]
   master_node = selection.selected_ids(mesh, "node")[0]

   # real constants of the coil (model_data.json, CONSTANT 1):
   #   cross_section 0.04 m^2, turns 1000, radius 0.5 m, resistance 50 Ohm, symmetry 1, T = (0, 0, 1)

   selection.select(mesh, "ID", master_node)
   selection.set_value_on_selection(mesh, "CURR", 100)          # step 1: coil current 100 A
   ...
   selection.set_value_on_selection(mesh, "CURR", 0)            # step 2: terminal voltage 10 kV
   selection.set_value_on_selection(mesh, "_VOLT", 10000)

What to look at in the results
------------------------------

* ``Jt`` (component 2) in the two rings: current density; in the stranded coil it is uniform.
* ``VOLT`` on the master node at step 1 (the voltage needed to drive 100 A) and ``EMF`` (the induced voltage).
* ``B``: the flux lines are the contours of :math:`r A_\theta` (``model_edit.poloidal_flux``).
* With iron: the *law operating range* of each step in the log (largest :math:`|B|` in the frame).

The macros
----------

coil_2D_avolt
~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_2D_avolt/macro.py
   :language: python

coil_2D_avolt_ferr
~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_2D_avolt_ferr/macro.py
   :language: python

coil_2D_avolt_timint
~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_2D_avolt_timint/macro.py
   :language: python

coil_2D_avoltemf
~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_2D_avoltemf/macro.py
   :language: python

coil_2D_avoltemf_ferr
~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_2D_avoltemf_ferr/macro.py
   :language: python

coil_2D_avoltemf_timint
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_2D_avoltemf_timint/macro.py
   :language: python

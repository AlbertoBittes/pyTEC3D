3D Coils on a Cyclic Sector
===========================

The 2D stranded-coil model of :doc:`coils_2d`, revolved about the y axis and modelled on a **30 degree sector** with
cyclic symmetry: a stranded coil (``TYPE 23``, lower ring, 2000 turns, 20 Ohm) and a passive massive ring (``TYPE 24``,
``SOLID_CONDUCTOR``, electric potential V). Sector 0 <= r <= 1 m, -1 <= y <= 1 m,
0 <= theta <= 30 degrees (Y-axis cylindrical system: theta from +x towards -z).

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Case
     - Element type
     - Mesh / variation
   * - ``coil_3D_avoltemf``
     - ``AV3D_nodal``
     - HEX8 + WEDGE6, nodal A with rotated frames on the 30 degree face.
   * - ``coil_3D_avoltemf_ferr``
     - ``AV3D_nodal``
     - Iron frame with B-H curve, dt = 1e5 s (DC steady states).
   * - ``coil_3D_avoltemf_timint``
     - ``AV3D_nodal_timint``
     - Time-integrated VOLT and EMF.
   * - ``coil_3D_avoltemf_edge``
     - ``AV3D_edge``
     - HEX20 + WEDGE15, edge A, tree-cotree gauge.
   * - ``coil_3D_avoltemf_edge_ferr``
     - ``AV3D_edge``
     - Edge elements with the iron frame.
   * - ``coil_3D_avoltemf_edge_timint``
     - ``AV3D_edge``
     - Same as ``_edge`` (its model data selects ``AV3D_edge``).
   * - ``coil_3D_avoltemf_edge_prisms``, ``_pyramids``
     - ``AV3D_edge``
     - All-prism and pyramid meshes (Nedelec on wedges and pyramids).

What these cases teach
----------------------

**Sector models and the symmetry real constant.** The coil real constants describe the **whole** coil (``volume =
2 pi 0.5 x 0.04`` m^3) and ``symmetry = 12`` (360 / 30): the circuit equation is then written for the whole coil, and
``CURR`` / ``VOLT`` are the current and voltage of the whole coil (:ref:`theory-stranded-coil`).

**Winding direction in an element coordinate system.** The coil cells use ``COORDINATE_SYSTEM = 1``, a Y-axis
cylindrical system, so that ``T = (0, 1, 0)`` is the toroidal direction in every cell:

.. code-block:: python

   selection.select(mesh, "MATERIAL", 23, location="elem")
   solver.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
   selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")

**Cyclic symmetry.** The nodes of the face theta = 0 are paired with their images at 30 degrees, coupling A and V; for
nodal elements the frames of the 30 degree face are rotated so that the same cylindrical components are coupled:

.. code-block:: python

   selection.select_pos(mesh, 1, (-tol, tol), coordinate_system=cyl)
   solver.model.periodic_boundaries(cyl, dy=30, array_name="CP_1", variables=("A", "V"))
   selection.select_pos(mesh, 1, (30 - tol, 30 + tol), coordinate_system=cyl)
   selection.set_value_on_selection(mesh, "node_angles", 30, 2)

**Nodal versus edge boundary conditions.** Nodal: A = 0 (all components) on the far faces and the components normal
to the axis on the axis. Edge: component 2 of ``_A`` on the far faces (tangential A), no axis condition, then
``tree_gauge(mesh, periodic_point_coupling="CP_1")``.

**Grounding.** The potential V of the passive ring floats: one node is fixed to 0.

The macros
----------

coil_3D_avoltemf
~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf/macro.py
   :language: python

coil_3D_avoltemf_ferr
~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_ferr/macro.py
   :language: python

coil_3D_avoltemf_timint
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_timint/macro.py
   :language: python

coil_3D_avoltemf_edge
~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_edge/macro.py
   :language: python

coil_3D_avoltemf_edge_ferr
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_edge_ferr/macro.py
   :language: python

coil_3D_avoltemf_edge_timint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_edge_timint/macro.py
   :language: python

coil_3D_avoltemf_edge_prisms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_edge_prisms/macro.py
   :language: python

coil_3D_avoltemf_edge_pyramids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/coil_3D_avoltemf_edge_pyramids/macro.py
   :language: python

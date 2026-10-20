CREATE Electromagnet
====================

Four models of a laboratory electromagnet (half models, ``symmetry = 2``): a racetrack stranded coil
(``STRANDED_COIL``, 156 or 127 turns) around a conducting iron yoke (``SOLID_CONDUCTOR`` with a B-H curve),
with a support plate and air.

.. list-table::
   :header-rows: 1
   :widths: 34 18 48

   * - Case
     - Element type
     - Analysis
   * - ``CREATE_B_measurement``
     - ``AV3D_edge``
     - Transient, coil driven by a measured voltage waveform (0.2 s, 2000 points), field measurement.
   * - ``CREATE_B_measurement_nodal``
     - ``AV3D_nodal``
     - The same with nodal elements.
   * - ``CREATE_force_measurement``
     - ``AV3D_edge``
     - Static force setup: coil current 5 A.
   * - ``CREATE_force_measurement_nodal``
     - ``AV3D_nodal``
     - The same with nodal elements.

What these cases teach
----------------------

**A coil made of pieces.** The racetrack is split into five real constants rows (``CONSTANT 1-5``): the straight parts
with their global winding direction ``T``, and the two arcs with ``T = (0, 1, 0)`` in their own cylindrical element
coordinate systems, centred on the arc centres:

.. code-block:: python

   solver.model.set_coordinate_system(91, XAxisCylindricalCoordinateSystem(origin=[xc0, yc0, zc0]))
   solver.model.set_coordinate_system(92, XAxisCylindricalCoordinateSystem(origin=[xc1, yc1, zc1]))
   selection.select(mesh, "CONSTANT", 2, location="elem")
   selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 91, location="elem")

All the coil nodes share one terminal (``CP_0 = 1``, VOLT and EMF), driven by ``_VOLT`` from the waveform.

**Conducting, saturating iron.** The yoke has eddy currents and a B-H curve at the same time: every step is a Newton
solve (:doc:`../user_manual/nonlinear`). The curve ends with a far point (B = 5 T) that sets the slope in deep
saturation.

**Edge or nodal.** The edge versions fix the tangential A on the whole boundary and gauge with a tree grown from it;
the nodal versions write symmetry-plane and far-field conditions component by component.

These are large models (170 000 - 250 000 cells, quadratic for the edge versions): run them in parallel.

The macros
----------

CREATE_B_measurement
~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/CREATE_B_measurement/macro.py
   :language: python

CREATE_B_measurement_nodal
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/CREATE_B_measurement_nodal/macro.py
   :language: python

CREATE_force_measurement
~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/CREATE_force_measurement/macro.py
   :language: python

CREATE_force_measurement_nodal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/CREATE_force_measurement_nodal/macro.py
   :language: python

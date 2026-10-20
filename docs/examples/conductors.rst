Wires and Loops
===============

Small models of a source conductor with an imposed current density (``STRANDED_CONDUCTOR``, element load ``Js``) next
to a passive conductor where eddy currents are induced.

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Case
     - Element type
     - Notes
   * - ``two_wires_3D_2order``
     - ``AV3D_nodal``, second order
     - Two parallel wires along z in a box (source at x < 0, passive at x > 0), quadratic nodal elements.
   * - ``two_wires_3D_edge``
     - ``AV3D_edge``
     - The same box with edge elements.
   * - ``3D_loop_cond_base``
     - ``AV3D_nodal``
     - Quarter model of a square source loop below an identical passive loop.

The first-order nodal and the HEX20 edge versions of the two wires, and the rings on 90 degree sectors, are MAPDL
benchmarks (``tests/apdl_benchmarks/two_wires_*``, ``sphere_coil_*``, see :ref:`examples-benchmarks`).

Two wires
---------

Boundary conditions of the box: A = 0 (tangential A for edge elements) on the faces x = +-1 m and y = +-0.5 m; the
faces z = +-0.5 m are crossed by the wires, so the current is normal to them and only the tangential components of A
vanish there. The passive wire has V = 0 on all its nodes.

3D loop
-------

``3D_loop_cond_base``: a quarter model of a square source loop below an identical passive solid loop. It shows named
boundary arrays (``X_BOUNDS_1``...), symmetry planes (tangential A = 0 where the current crosses the plane), the
potential fixed where the passive loop is cut by the symmetry planes, a vector Js given component by component on the
sides and the corner of the loop, and loads set after the model is loaded. The tutorial notebook
:doc:`../tutorials/tutorial_loops_3d` runs it step by step.

The macros
----------

two_wires_3D_2order
~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/two_wires_3D_2order/macro.py
   :language: python

two_wires_3D_edge
~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/two_wires_3D_edge/macro.py
   :language: python

3D_loop_cond_base
~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/examples/3D_loop_cond_base/macro.py
   :language: python

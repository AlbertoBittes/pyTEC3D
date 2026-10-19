Wires, Loops and Rings
======================

Small models of a source conductor with an imposed current density (``STRANDED_CONDUCTOR``, element load ``Js``) next
to a passive conductor where eddy currents are induced. They exercise every element type, first and second order,
static and transient steps, and the various ways of writing boundary conditions.

Two wires
---------

Two long parallel conductors in a box, the source at x < 0 and the passive one at x > 0.

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Case
     - Element type
     - Notes
   * - ``two_wires_2D``
     - ``AV2D_nodal`` planar
     - Static then transient step, anisotropic air.
   * - ``two_wires_3D``, ``_2order``
     - ``AV3D_nodal``
     - First / second-order nodal elements.
   * - ``two_wires_3D_edge``
     - ``AV3D_edge``
     - Edge elements without gauge.
   * - ``two_wires_3D_edge_el``, ``_2order``
     - ``AV3D_edge``
     - Box of HEX20 cells, solid rod (``SOLID_CONDUCTOR``), tree gauge call, several steps.

3D loop
-------

``3D_loop_cond_base``: a quarter model of a square source loop below an identical passive solid loop. It shows named
boundary arrays (``X_BOUNDS_1``...), symmetry planes (tangential A = 0 where the current crosses the plane), the
potential fixed where the passive loop is cut by the symmetry planes, a vector Js given component by component on the
sides and the corner of the loop, and loads set after the model is loaded.

Rings on 90 degree sectors
--------------------------

A source ring and a solid ring of small cross-section (r = 0.25 m, y = +-0.25 m) inside a sphere of radius 1 m:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Case
     - Notes
   * - ``ITER_coil_inf_2D``
     - The axisymmetric version; six steps alternating static and transient with different theta.
   * - ``ITER_coil_inf_3D``
     - Symmetry faces with rotated nodal frames; Js given in cylindrical components and converted with
       ``model_edit.map_vector_field_to_global_components``.
   * - ``ITER_coil_3D``
     - The other sector; Js in a cylindrical element coordinate system.
   * - ``ITER_coil_3D_test_1``
     - Symmetry faces from named arrays; Js computed per cell with NumPy.
   * - ``ITER_coil_3D_period``
     - Cyclic symmetry between the two faces (periodic_boundaries over -90 degrees), grounding of V.

Three ways of imposing a toroidal current density, from these cases:

.. code-block:: python

   # 1. element coordinate system (ITER_coil_3D): Js[1] is toroidal in every cell
   solver.model.set_coordinate_system(1, YAxisCylindricalCoordinateSystem())
   selection.select(mesh, "MATERIAL", 23, location="elem")
   selection.set_value_on_selection(mesh, "COORDINATE_SYSTEM", 1, location="elem")
   selection.set_value_on_selection(mesh, "Js", J, 1, location="elem")

   # 2. cylindrical components converted to global ones (ITER_coil_inf_3D)
   selection.set_value_on_selection(mesh, "Js", J, 1, location="elem")
   model_edit.map_vector_field_to_global_components(mesh, "Js", coordinate_system=YAxisCylindricalCoordinateSystem(),
                                                    location="elem")

   # 3. global components computed from the centroids (ITER_coil_3D_test_1)
   c = meshes.array_values(mesh, "cell_centroids", "elem"); js = meshes.array_values(mesh, "Js", "elem")
   rho = np.hypot(c[ids, 0], c[ids, 2]); js[ids, 0] = J * c[ids, 2] / rho; js[ids, 2] = -J * c[ids, 0] / rho

The macros
----------

two_wires_2D
~~~~~~~~~~~~

.. literalinclude:: ../../tests/two_wires_2D/macro.py
   :language: python

two_wires_3D
~~~~~~~~~~~~

.. literalinclude:: ../../tests/two_wires_3D/macro.py
   :language: python

two_wires_3D_2order
~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/two_wires_3D_2order/macro.py
   :language: python

two_wires_3D_edge
~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/two_wires_3D_edge/macro.py
   :language: python

two_wires_3D_edge_el
~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/two_wires_3D_edge_el/macro.py
   :language: python

two_wires_3D_edge_el_2order
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/two_wires_3D_edge_el_2order/macro.py
   :language: python

3D_loop_cond_base
~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/3D_loop_cond_base/macro.py
   :language: python

ITER_coil_inf_2D
~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_coil_inf_2D/macro.py
   :language: python

ITER_coil_inf_3D
~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_coil_inf_3D/macro.py
   :language: python

ITER_coil_3D
~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_coil_3D/macro.py
   :language: python

ITER_coil_3D_test_1
~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_coil_3D_test_1/macro.py
   :language: python

ITER_coil_3D_period
~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_coil_3D_period/macro.py
   :language: python

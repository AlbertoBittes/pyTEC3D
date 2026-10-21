3D Coil Sectors
===============

The rings of :doc:`coils_2d`, revolved about the y axis and modelled on a 30 degree cyclic sector. The lower ring is a
multi-turn stranded coil fed through its terminals (VOLT, EMF); the upper ring is a passive massive conductor (A and
the electric scalar potential V). Four variants cover the 3D element families:

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Case
     - Element type
     - What it shows
   * - ``coil_3D_avoltemf``
     - ``AV3D_nodal``
     - Selections in a cylindrical coordinate system, an element coordinate system for the winding direction
       (T = toroidal), far-field, axis and cyclic-symmetry conditions on a nodal vector potential (nodal frames of the
       second face rotated), the real constants of a coil sector (whole-coil volume, symmetry).
   * - ``coil_3D_avoltemf_edge``
     - ``AV3D_edge``
     - Edge (Nedelec) elements for A: edge DOFs on the midside nodes, the tree-cotree gauge
       (:func:`pytec3d.tree_gauge.tree_gauge`), cyclic symmetry of edge values (no frame rotation).
   * - ``coil_3D_avoltemf_edge_ferr``
     - ``AV3D_edge`` + B-H law
     - Edge elements with a saturating iron frame: the gauge with a nonlinear law, Newton-Raphson, regions defined by
       editing TYPE and MATERIAL before the model is loaded.
   * - ``coil_3D_avoltemf_timint``
     - ``AV3D_nodal_timint``
     - Time-integrated electric unknowns (VOLT, EMF, V as time integrals; the EMF integral is the linked flux): a
       symmetric system, the meaning of the inputs of a voltage-driven step.

Each folder holds the mesh, its model data and the commented script; the benchmark repository compares the same
geometries with ANSYS MAPDL. :doc:`../tutorials/tutorial_coil_3d_sector` and :doc:`../tutorials/tutorial_edge_elements`
run the nodal and edge sectors step by step.

coil_3D_avoltemf
----------------

.. literalinclude:: ../../examples/coil_3D_avoltemf/macro.py
   :language: python

coil_3D_avoltemf_edge
---------------------

.. literalinclude:: ../../examples/coil_3D_avoltemf_edge/macro.py
   :language: python

coil_3D_avoltemf_edge_ferr
--------------------------

.. literalinclude:: ../../examples/coil_3D_avoltemf_edge_ferr/macro.py
   :language: python

coil_3D_avoltemf_timint
-----------------------

.. literalinclude:: ../../examples/coil_3D_avoltemf_timint/macro.py
   :language: python

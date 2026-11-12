Scalar Potentials and Superconductors
=====================================

Five analyses with the scalar-potential element types (:ref:`el-phi3d`, :ref:`el-tphi3d`, :ref:`el-tphi2d`):
magnetostatics with the reduced potential and its source field in one system, the T-phi (H) formulation of eddy
currents with cuts and fixed-T coils, and superconductors with power laws, in 3D, 2D and axisymmetric.

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - Case
     - Element type
     - What it shows
   * - ``core_3D_phi``
     - ``PHI3D_nodal``
     - A coil (``Js``) on a closed iron frame, static: the source field Hs solved with PHI in one system, its tree
       gauge, a B-H curve (``bh_curve_h``), Newton on static load steps; no cut for an iron loop around a current.
   * - ``ring_3D_tphi``
     - ``TPHI3D_edge``
     - A conducting ring in a field switched on: the gauges of T-phi, the applied field as a boundary potential, a cut
       layer spanning the hole with its crossing edges coupled on one DOF (the ring current).
   * - ``coil_core_3D_tphi``
     - ``TPHI3D_edge``
     - A current-driven coil as T fixed on its edges (T0 over the winding and "the disk filling the hole"), a
       conducting core in the hole (T0 moved out of it by a gradient), :func:`pytec3d.tree_gauge.edges` and
       :func:`pytec3d.tree_gauge.edge_circulations`, the flux of a turn.
   * - ``wire_2D_tphi``
     - ``TPHI2D_edge``
     - The AC transport loss of a superconducting round wire against Norris's formula: a power law as a table with a
       normal-state end, the nonlinear settings for high n, the cut of a 2D conductor and the jump of PHI on the box,
       the E output.
   * - ``bulk_2D_tphi_axisym``
     - ``TPHI2D_edge`` (axisymmetric)
     - The field-cooled magnetization of a REBCO bulk: ``power_law_b`` (Kim model Jc(B)), a static step (the cooled
       state) followed by transient steps, the trapped field and its creep; no cut (the section touches the axis).

The meshes are structured grids built by ``examples/files/formulation_meshes.py`` (only numpy and VTK); each folder
holds its mesh, its model data and the commented script. The tutorials :doc:`../tutorials/tutorial_scalar_potential`,
:doc:`../tutorials/tutorial_tphi_cuts` and :doc:`../tutorials/tutorial_superconductors` run them step by step, with
variants (drives of a cut, Hs given on the edges, exponent sweeps). The benchmarks of these formulations (analytic and
against A, V-A or MAPDL) are the ``formulations`` suite of the benchmark repository.

core_3D_phi
-----------

.. literalinclude:: ../../examples/core_3D_phi/macro.py
   :language: python

ring_3D_tphi
------------

.. literalinclude:: ../../examples/ring_3D_tphi/macro.py
   :language: python

coil_core_3D_tphi
-----------------

.. literalinclude:: ../../examples/coil_core_3D_tphi/macro.py
   :language: python

wire_2D_tphi
------------

.. literalinclude:: ../../examples/wire_2D_tphi/macro.py
   :language: python

bulk_2D_tphi_axisym
-------------------

.. literalinclude:: ../../examples/bulk_2D_tphi_axisym/macro.py
   :language: python

ITER Plasma Equilibrium and Vertical Displacement
=================================================

The full workflow of the plasma algorithm (:doc:`../user_manual/plasma_algorithm`) on ITER: a free-boundary
equilibrium in 2D axisymmetry, then a vertical displacement event (the poloidal beta drops from 0.68 to 0.24 in
0.636 s at about 15 MA: the plasma moves up and becomes limited around 0.36 s) with the eddy currents of the vacuum
vessel and of the passive structures, in 2D and on thin 3D toroidal sectors. The waveforms end with the thermal
quench at 0.6356 s.

Two notebooks run the 2D cases step by step with plots and analysis of the results:
:doc:`../tutorials/tutorial_iter_equilibrium` and :doc:`../tutorials/tutorial_iter_evolution`; a third one,
:doc:`../tutorials/tutorial_iter_3d_halo`, runs ``ITER_3D_evol_coils`` (3D sector, voltage-driven coils, halo
currents) with 30 ms steps.

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Case
     - Element type
     - Content
   * - ``ITER_2D``, ``_fine``, ``_fine2``, ``_coils``, ``_carma``
     - ``AV2D_nodal``
     - Equilibria on meshes of different refinement and order (``_coils``, ``_carma`` second order).
   * - ``ITER_2D_evol``, ``_fine``, ``_fine2``, ``_carma``
     - ``AV2D_nodal``
     - Vertical displacement from the matching equilibrium; coils with fixed current densities.
   * - ``ITER_2D_evol_coils``
     - ``AV2D_nodal``
     - Coils as voltage-driven stranded coils (circuit equations).
   * - ``ITER_3D_evol``, ``_fine``, ``_fine2``
     - ``AV3D_nodal_timint``
     - Vertical displacement on 2, 20 and 1 degree sectors, with halo currents from the limiter phase.
   * - ``ITER_3D_evol_edge``, ``_fine2_edge``
     - ``AV3D_edge_timint``
     - The same with edge elements.
   * - ``ITER_3D_evol_coils``, ``_coils_edge``
     - ``AV3D_nodal``, ``AV3D_edge``
     - 3D sector with voltage-driven stranded coils.

Workflow
--------

#. **Equilibrium** (``ITER_2D``). Coil current densities ``Js = I / a`` in the CS modules and PF coils, a uniform
   15 MA in the vacuum chamber as initial guess, A = 0 on the boundary, static. ``Algorithm.equilibrium_solve()``
   iterates: plasma boundary from the flux, Grad-Shafranov current density scaled to 15 MA, control-coil corrections to
   keep the magnetic axis at the target, field solve, until the flux converges. ``write_output()`` stores it.
#. **Evolution** (``ITER_2D_evol``). A mesh with the same cells, where the vessel and structures are
   ``SOLID_CONDUCTOR`` cells (closed toroidal shells, in 3D with the potential V). The equilibrium output becomes the
   vacuum-chamber mesh of the algorithm and its ``Jt`` (equal to ``Js`` on the stranded conductors) the source ``Js``
   of the evolution. A static first step at ``t = 1e-10`` starts from the equilibrium; then,
   for each point of the waveforms, a transient step with the plasma current and beta of that time and one
   ``evolution_solve()``.
#. **3D** (``ITER_3D_evol``). The 2D mesh revolved over 2 degrees with cyclic symmetry; the coil current densities of
   the equilibrium are written as the toroidal component of ``Js`` in a cylindrical element coordinate system; the
   algorithm maps the 2D vacuum chamber to the toroidal stations (``extrusion_angle``, ``n_toroidal_divisions``) and,
   from ``t0_halo``, injects halo currents into the wall (``AMPS``).
#. **Voltage-driven coils** (``*_coils``). Each coil a stranded coil with its terminal (VOLT, EMF coupled), driven by the
   steady-state voltage of its equilibrium current: the coil currents then respond to the displacement.

These runs are long (hundreds of steps, several Picard iterations each).

The macros
----------

ITER_2D
~~~~~~~

.. literalinclude:: ../../tests/ITER_2D/macro.py
   :language: python

ITER_2D_evol
~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_2D_evol/macro.py
   :language: python

ITER_2D_evol_coils
~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_2D_evol_coils/macro.py
   :language: python

ITER_3D_evol
~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_3D_evol/macro.py
   :language: python

ITER_3D_evol_edge
~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_3D_evol_edge/macro.py
   :language: python

ITER_3D_evol_coils
~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_3D_evol_coils/macro.py
   :language: python

ITER_3D_evol_coils_edge
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../tests/ITER_3D_evol_coils_edge/macro.py
   :language: python

The other variants (``ITER_2D_fine``, ``_fine2``, ``_coils``, ``_carma``, ``ITER_2D_evol_fine``, ``_fine2``, ``_carma``,
``ITER_3D_evol_fine``, ``_fine2``, ``_fine2_edge``) follow the same scripts; their headers list the differences.

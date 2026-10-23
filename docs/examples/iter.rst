ITER Plasma Equilibrium and Vertical Displacement
=================================================

The full workflow of the plasma algorithm (:doc:`../user_manual/plasma_algorithm`) on ITER: a free-boundary
equilibrium in 2D axisymmetry, then a vertical displacement event, in 2D and on thin 3D toroidal sectors. During the
event the poloidal beta drops from 0.68 to 0.24 in 0.636 s at about 15 MA; the plasma moves up and becomes limited
around 0.36 s. The eddy currents of the vacuum vessel and of the passive structures react to it. The waveforms end
with the thermal quench at 0.6356 s.

The examples cover different meshes and element families. The standard cases (``ITER_2D_evol``, ``ITER_2D_coils``,
``ITER_3D_evol``, ...) are part of the benchmark of the plasma algorithm against APEC in the benchmark repository
(:doc:`../testing`).

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Case
     - Element type
     - Content
   * - ``ITER_2D_carma``
     - ``AV2D_nodal`` (quadratic)
     - Equilibrium on a second-order mesh (QUAD8 / TRI6): the starting point of ``ITER_2D_evol_carma`` and
       ``ITER_3D_evol_coils_edge``.
   * - ``ITER_2D_fine2``
     - ``AV2D_nodal`` (linear)
     - Equilibrium on a refined linear mesh: the starting point of ``ITER_2D_evol_fine2`` and ``ITER_3D_evol_fine2``.
   * - ``ITER_2D_evol_carma``
     - ``AV2D_nodal`` (quadratic)
     - Vertical displacement from ``ITER_2D_carma``; coils with fixed current densities.
   * - ``ITER_2D_evol_fine2``
     - ``AV2D_nodal`` (linear)
     - Vertical displacement from ``ITER_2D_fine2``.
   * - ``ITER_3D_evol_fine2``
     - ``AV3D_nodal_timint``
     - Vertical displacement on a 1 degree sector of the fine2 mesh (time-integrated V), with halo currents
       injected into the vessel wall from the limiter phase.
   * - ``ITER_3D_evol_coils_edge``
     - ``AV3D_edge``
     - Vertical displacement on the ``ITER_2D_carma`` mesh revolved over 2 degrees in 4 cells (quadratic), with edge elements and a tree gauge, the coils as voltage-driven stranded coils (their
       currents respond to the displacement), halo currents.

Three notebooks run these analyses step by step, with plots of the results:

* :doc:`../tutorials/tutorial_iter_equilibrium` (``ITER_2D_fine2``);
* :doc:`../tutorials/tutorial_iter_evolution` (``ITER_2D_evol_fine2``);
* :doc:`../tutorials/tutorial_iter_3d_halo` (``ITER_3D_evol_coils_edge``, coarser steps).

Workflow
--------

#. **Equilibrium** (``ITER_2D_carma``, ``ITER_2D_fine2``). The macro sets:

   * coil current densities ``Js = I / a`` in the CS modules and PF coils;
   * a uniform 15 MA in the vacuum chamber as initial guess;
   * A = 0 on the boundary, and a static solve.

   ``Algorithm.equilibrium_solve()`` iterates until the flux converges. Each iteration finds the plasma boundary from
   the flux, scales the Grad-Shafranov current density to 15 MA, corrects the control coils so that the magnetic axis
   stays at the target, and solves the field. ``write_output()`` stores the result.
#. **2D evolution** (``ITER_2D_evol_*``). The evolution mesh has the same cells as the equilibrium mesh, but the
   vessel and structures are ``SOLID_CONDUCTOR`` cells (closed toroidal shells). The equilibrium output becomes the
   plasma-region mesh of the algorithm, and its ``Jt`` (equal to ``Js`` on the stranded conductors) becomes the source
   ``Js`` of the evolution. A static first step at ``t = 1e-10`` starts from the equilibrium. Then each point of the
   waveforms gives one transient step, with the plasma current and beta of that time and one ``evolution_solve()``.
#. **3D evolution** (``ITER_3D_evol_*``). The setup:

   * the 2D mesh is revolved over a small sector with cyclic symmetry;
   * the coil current densities of the equilibrium are written as the toroidal component of ``Js``, in a cylindrical
     element coordinate system;
   * the algorithm maps the 2D vacuum chamber to the toroidal stations (``extrusion_angle``, ``n_toroidal_divisions``)
     and averages the flux over them;
   * from ``t0_halo``, it injects halo currents into the wall (``AMPS``).

   With edge elements (``ITER_3D_evol_coils_edge``) a tree-cotree gauge replaces the Coulomb gauge.
#. **Voltage-driven coils** (``ITER_3D_evol_coils_edge``). Each coil is a stranded coil with its terminal (VOLT, EMF
   coupled), driven by the steady-state voltage of its equilibrium current. The coil currents then respond to the
   displacement.

Run the equilibrium first: an evolution reads the most recent ``<equilibrium>_out_*.vtu`` of its equilibrium folder
(``../ITER_2D_fine2/output``, ...). These runs are long: hundreds of steps, with several Picard iterations each.

The macros
----------

ITER_2D_fine2
~~~~~~~~~~~~~

.. literalinclude:: ../../examples/ITER_2D_fine2/macro.py
   :language: python

ITER_2D_evol_fine2
~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../examples/ITER_2D_evol_fine2/macro.py
   :language: python

ITER_3D_evol_fine2
~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../examples/ITER_3D_evol_fine2/macro.py
   :language: python

ITER_3D_evol_coils_edge
~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../examples/ITER_3D_evol_coils_edge/macro.py
   :language: python

``ITER_2D_carma`` and ``ITER_2D_evol_carma`` follow the scripts of ``ITER_2D_fine2`` and ``ITER_2D_evol_fine2`` on
the quadratic mesh; their headers list the differences (mesh, order, files).

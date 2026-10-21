Plasma Equilibrium and Evolution
================================

:class:`pytec3d.algorithm.Algorithm` (APEC-M) computes free-boundary plasma equilibria and follows their evolution in time,
coupled with the finite-element field solution of the coils and of the conducting structures. It only uses the public
API of the solver: it writes the plasma current density into ``Js``, calls ``solver.solve()`` and reads the result.
The notebooks :doc:`../tutorials/tutorial_iter_equilibrium` and :doc:`../tutorials/tutorial_iter_evolution` use it
step by step on ITER.

The vacuum-chamber mesh
-----------------------

The plasma lives on a 2D (R, Z) mesh of the region where it can be, the cells with ``VACUUM_CHAMBER = 1``. The solver
mesh is the same 2D axisymmetric mesh, or a 3D toroidal sector made by revolving it about the y axis. In 3D the
algorithm builds toroidal maps between the two meshes: the plasma current density is projected on the cells of every
toroidal station, the flux is averaged back over the stations.

* For an equilibrium, the vacuum-chamber mesh is taken from the solver mesh (the VACUUM_CHAMBER cells).
* For an evolution, set it from the output of the equilibrium, which carries the flux ``PSI``, ``Js`` and the plasma
  cells: ``algorithm.vacuum_chamber_mesh = Model.read_model(path).mesh`` (or
  ``mesh.read_mesh(path, reset_arrays=False)``), with ``path`` the last ``output/<name>_out_<n>.vtu`` of the
  equilibrium (``_1``: the equilibrium writes one output, after its Picard loop).
* The test macros select the vacuum-chamber cells and their nodes on the solver mesh before the evolution solves
  (``selection.select(mesh, "VACUUM_CHAMBER", 1, location="elem"); selection.nsle(mesh)``). This is not required: the
  algorithm makes its own selections on the solver mesh and on its copy of the vacuum chamber.
* The current densities of the equilibrium (coils and plasma) are in its ``Js`` and ``Jt`` (equal on the stranded
  conductors): either can be copied into the ``Js`` of the evolution model.
* The ``PLASMA`` cell array marks the cells that carry plasma current: the cells connected to the axis inside the
  X-point flux surface, and inside the boundary of the current (the limiter flux surface of a limited plasma, the
  halo boundary when a halo is active).
* In the evolution models the vacuum vessel and the passive structures are ``SOLID_CONDUCTOR`` cells: in 2D closed
  toroidal shells (no electric DOF), in 3D conductors with the potential V (fixed on one node, coupled between the
  faces of the sector).
* The position control is proportional: the correction currents are recomputed at every iteration from the residual
  field at the target and replace the previous correction. The axis converges near the target, not on it
  (:ref:`theory-plasma` *Position control*).

Required arrays: cell arrays ``VACUUM_CHAMBER``, the component arrays of the control coils (e.g. ``PFC_2``),
``cell_centroids``, ``cell_volume``, ``Js``; point arrays ``A``, ``PSI``. ``AMPS`` (halo currents) is used when the model
has 3D conductors.

The Picard iteration
--------------------

Each iteration (``equilibrium_solve()`` and ``evolution_solve()``):

#. finds the magnetic axis (maximum of PSI), the X point (saddle of PSI) and the limiter point, hence the last closed
   flux surface and the separatrix configuration; the cells inside it, connected to the axis, form the plasma
   (``PLASMA`` cell array);
#. computes the plasma current density from the Grad-Shafranov profile

   .. math::

      J_\varphi = \lambda\,\big(1-\psi_n^{\alpha_m}\big)^{\alpha_n}
                  \Big(\beta_p\frac{R}{R_0}+(1-\beta_p)\frac{R_0}{R}\Big),
      \qquad \psi_n=\frac{\psi-\psi_\text{axis}}{\psi_\text{boundary}-\psi_\text{axis}},

   with :math:`\lambda` scaling the total current to ``plasma_current``, and the poloidal current density
   :math:`F'(\psi)\mathbf{B}_p` from :math:`F^2=F_0^2+\dots` with :math:`F_0` = ``AF0`` (:math:`R_0B_{\varphi0}`);
#. **equilibrium only**: corrects the currents of the control coils, :math:`I = g_r I_r + g_v I_v` for each coil of
   ``control_coils``, with the radial and vertical correction currents :math:`I_r, I_v` that cancel the poloidal
   field at the target axis (from the field per unit current ``bradx/brady/bverx/bvery_value``);
#. projects Js (and the halo currents) on the solver mesh, solves, computes and averages the new flux PSI;
#. stops when the largest change of PSI is below ``tolerance_value``; raises
   :class:`~pytec3d.errors.ConvergenceError` after ``max_picard_iterations``.

The first evolution call only initializes the flux (one solve). Each Picard iteration is one ``solve()`` at the
same time: the time step is solved again with the updated plasma current.

Input parameters
----------------

``algorithm.input_parameters`` (defaults in :data:`pytec3d.algorithm.DEFAULT_INPUT_PARAMETERS`):

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Parameter
     - Meaning
   * - ``target_central_axis_x``, ``_y``, ``_z``
     - Target position of the magnetic axis (x = R, y = Z) for the equilibrium control.
   * - ``R0``
     - Major radius [m] of the profile.
   * - ``alpha_mp``, ``alpha_np``
     - Exponents :math:`\alpha_m, \alpha_n` of the profile.
   * - ``beta_p``
     - Poloidal beta.
   * - ``AF0``
     - :math:`F_0 = R_0 B_{\varphi 0}` [T m], the vacuum toroidal field function.
   * - ``plasma_current``
     - Total plasma current [A].
   * - ``max_picard_iterations``, ``tolerance_value``
     - Iteration limit and tolerance on the change of PSI.
   * - ``bradx_value``, ``brady_value``, ``bverx_value``, ``bvery_value``
     - Field (R and Z components) at the target axis of a unit radial / vertical correction current.
   * - ``control_coils``
     - ``[{"component": "PFC_2", "radial_gain": 1.0, "vertical_gain": 0.0}, ...]``: the coils (cell arrays) that
       receive the correction currents, as a current density ``I / cross_section`` added to their Js.
   * - ``extrusion_angle``, ``n_toroidal_divisions``
     - 3D solver mesh: toroidal extent [deg] and number of cells along the torus.
   * - ``halo_delta1``, ``halo_delta2``, ``t0_halo``
     - Halo region from ``t0_halo`` on: :math:`\psi_{2,\min}=\psi_{\min}-\delta_1(1+(t-t_0)\delta_2)(\psi_{\max}-\psi_{\min})`;
       the halo currents are injected as ``AMPS`` along the wall (3D).
   * - ``equilibrium_mesh``, ``static_mesh``, ``evolution_mesh``
     - File names, for reference.

Keys not in the defaults are ignored (e.g. ``min_distance``, ``n_slices`` in some test macros).

Results
-------

``algorithm.status`` holds the results of the last iteration: axis, X point and limiter positions and fluxes, the
separatrix configuration, the correction currents, the plasma current scaling and cell count, convergence. Every
iteration appends a row to ``output/<model>_iteration_diagnostics.csv``.

A complete workflow
-------------------

#. **Equilibrium** (``examples/ITER_2D_fine2``): 2D axisymmetric, static; coil current densities in the coils, a uniform initial
   plasma current density in the vacuum chamber; ``algorithm.equilibrium_solve()``; ``solver.write_output()``.
#. **Evolution** (``examples/ITER_2D_evol_fine2``, ``examples/ITER_3D_evol_fine2``): read the equilibrium output as the vacuum-chamber
   mesh, copy its current densities (coils, plasma) into ``Js`` of the evolution model (same cells), set the
   conductors (vessel, structures) as conducting domains; then a static first step and transient steps following the
   plasma current and beta waveforms, one ``evolution_solve()`` per step. Coils can be fixed current densities
   (stranded conductors) or voltage-driven stranded coils (``examples/ITER_3D_evol_coils_edge``).

See :doc:`../examples/iter`.

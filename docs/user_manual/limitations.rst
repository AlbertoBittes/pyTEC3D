Limitations and Open Points
===========================

What pyTEC3D does not do, does only approximately, or leaves to the user. Each point says what to do about it.

General
-------

* **One element type per model**, any of its domains; one FE order for all the Lagrange fields of a model; edge
  (Nedelec) fields are first order only (:doc:`model_data`).
* **No restart from a saved state**: a new session starts from zero fields. Run an analysis in one session, or rebuild
  the initial state with a static step.
* **Time schemes**: theta method and Newmark, fixed steps set by the script; no adaptive time stepping.
* **Linear solver**: MUMPS (direct) is the practical choice for every formulation; iterative solvers are available
  through PETSc options but not tuned (:doc:`solver_settings`).
* **Element results are cell means** (:doc:`outputs`). Integrals of products computed from them (Joule power, AC loss
  :math:`\sum \mathbf{E}\cdot\mathbf{J}\,V`) are slightly low where the fields vary inside a cell (2 % on the ring of
  :doc:`../tutorials/tutorial_tphi_cuts`, 4 cells across): refine where such integrals matter.

Nonlinear laws
--------------

* **Tabulated, isotropic, stateless**: the curves are piecewise linear tables; no hysteresis (iron, magnets) and no
  history-dependent law; a law depends on the magnitude of its input (:doc:`nonlinear`).
* **Law inputs** are fields at the quadrature points (``bh_curve``: B; ``power_law_b``: J and H); no element means,
  circuit or global quantities, no temperature.
* **Jacobian without cross-coupling**: where a law's output depends on another unknown than its own input
  (``power_law_b``: :math:`\rho(J, B)` through B), the Jacobian holds only the derivative with respect to the own
  input. Newton converges to the right solution, more slowly (linearly where that coupling matters).
* **Superconductors**: ``power_law_b`` takes :math:`|B|=\mu_0|H|` (a non-magnetic conductor); no anisotropic
  :math:`J_c(B_\parallel, B_\perp)`, no :math:`n(B)`, no thermal model or quench. The nonlinear settings for high n
  (an absolute residual tolerance scaled to the problem, ``require_residual_reduction`` on or off) are the user's
  choice: :ref:`superconductor-settings`. A non-converged step raises ``ConvergenceError``.
* The 3D A, V-A types with ``power_law_e`` are not reliable at real superconductor parameters (no convergence or
  oscillations): use the T-phi types.

Scalar-potential element types
------------------------------

* **PHI3D_nodal is static.** Its source field Hs adds one unknown per edge of the whole mesh (the alternative: Hs given
  on every edge). Hs carries the solenoidal part of ``Js``: a source current that is not divergence-free is used
  without a warning.
* **T-phi has no source current and no circuits.** A coil is T fixed on its edges (current drive; its voltage comes
  afterwards from the flux of a turn), T0 is built by the user for each winding. A conductor's net current is its cut
  DOF: series connections (shared coupling) and parallel connections at the same voltage (same ``VLTG``) are
  possible; a total current shared between parallel conductors, a resistor or supply impedance in series, an
  external network are not. Coils with circuits are an A, V-A model.
* **Cuts are made by hand** (:ref:`el-cuts`): the layer, the classification of its edges, their orientation and the
  couplings. **A missing or broken cut is not detected**: the conductor is then modelled as cut open (no net
  current). One free cut DOF per independent hole; a second free cut on the same loop makes the system singular.
* **Static steps are managed by the user**: ``_T`` = 0 in the conductors and ``_PHI`` free inside on a static step,
  the transient gauges restored afterwards (:doc:`boundary_conditions`). A static step wipes the currents of the
  conductors.
* **Accuracy is first order** (no second-order Nedelec): the error falls linearly with the cell size. Measured on the
  benchmarks: 0.3-3 % on fields of magnetic and conducting balls, 3 % (B) and 6 % (core current) on a coil around a
  conducting core with 4-cell regions, 8-9 % on the current profile of a superconducting slab with 40 cells; in
  conducting iron with :math:`\mu_r` = 100 the coarse-mesh difference to the A formulation reaches 20 % (10 % with
  half the cell size: the two formulations bracket the solution). Refine, and check with a finer mesh.
* **No thin-conductor formulations**: no T-A (tapes as lines in a 2D mesh, shells in 3D), no element types of a lower
  dimension than the mesh.

Validation coverage
-------------------

The ``formulations`` suite of the benchmark repository (:doc:`../testing`) checks the scalar-potential types against
analytical solutions, against the A, V-A types on the same mesh and, where MAPDL has the formulation, against MAPDL
(the reduced scalar potential of a magnetic ball, the eddy currents of a conducting ball). Not covered yet: a
comparison with an independent H-formulation code on the published superconductor benchmarks, a transient T-phi
model with a B-H curve against a reference, several conductors with cuts in series or in parallel, ``TPHI2D_edge``
with a B-H curve in the suite (checked by hand against the exact solution of an iron cylinder).

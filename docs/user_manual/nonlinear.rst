Nonlinear Materials
===================

Nonlinear laws
--------------

A nonlinear law replaces coefficients of the element types by values computed, at each quadrature point, from a field
of the solution and from the data of the law. The laws are defined in ``element_library/pytec3d_element_library/data/nonlinear_laws.json``
and compiled into the solver. The library has these laws:

.. list-table::
   :header-rows: 1
   :widths: 14 10 10 16 50

   * - Law
     - Source
     - Input
     - Data
     - Replaces
   * - ``bh_curve``
     - ``mat``
     - B
     - ``B``, ``H`` arrays
     - ``nu`` (reluctivity tensor) and ``nu_gauge`` (Coulomb gauge coefficient of the 3D nodal types): secant
       :math:`\nu_s=H(|B|)/|B|` in the residual, tangent
       :math:`\nu_s\mathbf{I}+(H'(|B|)-\nu_s)\hat{\mathbf{b}}\otimes\hat{\mathbf{b}}` in the Jacobian.
   * - ``bh_curve_h``
     - ``mat``
     - H
     - curve ``H``, ``B``
     - ``mu`` of ``TPHI3D_edge`` and ``PHI3D_nodal``: :math:`\mu_s=B(|H|)/|H|` in the residual,
       tangent :math:`\mu_s\mathbf{I}+(B'(|H|)-\mu_s)\hat{\mathbf{h}}\otimes\hat{\mathbf{h}}` in the Jacobian. A law of its
       own, not the inverse of ``nu``: ``set_bh_curve`` selects it for these types, from the same table. In a magnet
       :math:`\mathbf{B}=\mu_s\mathbf{H}+\mathbf{B}_r`, the remanence a constant load. In the
       :math:`\nabla\cdot\mathbf{B}=0` equation it is an ordinary secant and tangent (static steps too); in a conducting
       iron the Faraday term is the rate of B by the time scheme (backward Euler :math:`(B(H)-B(H_n))/\Delta t`; any
       ``theta``: the rates of B are stored per element, :ref:`theory-nonlinear-elements`).
   * - ``power_law_e``
     - ``mat``
     - E
     - curve ``J``, ``E``
     - ``sigma`` of the 3D A, V-A types: :math:`\sigma_s=J(|E|)/|E|` in the residual, tangent
       :math:`\sigma_s\mathbf{I}+(J'(|E|)-\sigma_s)\hat{\mathbf{e}}\otimes\hat{\mathbf{e}}` in the Jacobian.
   * - ``power_law``
     - ``mat``
     - Jt
     - curve ``J``, ``E``
     - ``rho`` of ``TPHI3D_edge`` and ``TPHI2D_edge``: :math:`\rho_s=E(|J|)/|J|` in the residual, tangent
       :math:`\rho_s\mathbf{I}+(E'(|J|)-\rho_s)\hat{\mathbf{j}}\otimes\hat{\mathbf{j}}` in the Jacobian.
   * - ``power_law_b``
     - ``mat``
     - Jt, H
     - curves ``J``, ``E`` and ``Bjc``, ``jc``
     - ``rho`` of the T-phi types with a critical current depending on the field:
       :math:`E(J,B)=E_0(J/s(|B|))`, :math:`E_0` the E-J curve at B = 0, :math:`s=J_c(B)/J_c(0)` the second curve,
       :math:`|B|=\mu_0|H|` (a non-magnetic conductor), e.g. the Kim model :math:`s=1/(1+|B|/B_0)`. The tangent is
       at fixed B: the Jacobian has no term of the variation of B with the unknowns (no cross-coupling), so Newton
       needs more iterations where s varies fast (``kim_slab``: up to 18 per step). Its floor (the first segment of
       :math:`E_0`) is scaled by 1/s.

Both conduction laws read the same E-J table, rows ``[J, E]`` from J = 0, both increasing:
``solver.model.set_power_law(material_id, table)`` selects ``power_law`` for the T-phi types and ``power_law_e`` for
the 3D A, V-A types; ``set_power_law(material_id, table, jc_b=rows)`` selects ``power_law_b`` (rows
``[|B|, Jc(B) / Jc(0)]``, ``|B|`` increasing). The law is the table, nothing else: a superconductor's power law with a normal-state floor is
``E = rho0 * J + Ec * (J / Jc)**n``, sampled densely (the curve is piecewise linear: with log-spaced points of relative
spacing :math:`\delta` its error is about :math:`n^2\delta^2/8`), and its first segment sets the resistivity at zero
current.

The **source** tells which id selects the law: ``mat`` a MATERIAL id, ``rc`` a CONSTANT id, ``me`` a TYPE id. A law
applies on the cells with that id, in every element type and domain that has a coefficient it replaces, restricted to
the ``element_types`` of the law when it lists them (``bh_curve``: all of them; the two power laws: the types above).
Where a law is selected the coefficient is its output, elsewhere the linear value. A model with at least one such cell is solved with Newton's method at every step
(:ref:`theory-newton`).

.. warning::

   On the slab of the ``formulations`` suite at real superconductor parameters (:math:`n = 10`--:math:`20`, no resistivity
   floor) only ``TPHI3D_edge`` is accurate: the A, V-A types with ``power_law_e`` fail to converge or oscillate there.

A law may read more than one field: ``"inputs"`` in ``nonlinear_laws.json`` lists extra fields read at the
quadrature points, values only (no Jacobian of their variation), and ``"curves"`` extra named curves, each with its
own two arrays in the law's entry (``power_law_b``). The full grammar is in ``solver/tools/generate_registries.py``.

.. _superconductor-settings:

Superconductors: the nonlinear settings
---------------------------------------

Measured on the ``formulations`` benchmarks (``powerlaw_hts`` at n = 20 to 100, ``bean_slab``, ``norris_wire``) and
the examples ``wire_2D_tphi`` and ``bulk_2D_tphi_axisym`` (:doc:`../tutorials/tutorial_superconductors`):

* **Use an absolute residual tolerance** scaled to the problem, with the relative one off:
  ``{"absolute_residual_tolerance": tol, "relative_residual_tolerance": 1e-30}``, ``tol`` a small fraction of the
  residual of a T row at the critical state, about :math:`E_c\,|\nabla\times\mathbf{w}|\,V_\text{cell}` (in V m:
  :math:`E_c` times a face area). The relative tolerance is relative to the first residual of the step, which is huge
  when the starting guess has J above Jc anywhere (:math:`E\propto(J/J_c)^n`): at n = 100 a relative 1e-8 stopped on
  states with :math:`E\approx10^{20}E_c`. A loose absolute tolerance stops at once and gives a wrong loss (``bean_slab``:
  1e-6 against a T-row scale of 4e-10 gave 2500 times the loss; 1e-12 to 1e-16 all give the same result).
  Check the result: :math:`\max|E|` against :math:`E_c` (the ``E`` output) tells a wrong state at once.
* **Require the residual to decrease** (``"require_residual_reduction": True``) when the loads change a lot in one
  step at a few cells, e.g. a transport current through a cut: the first iterate has J far above Jc next to the cut, and
  a full Newton step without that check reaches states where the factorization fails (``norris_wire``; the example
  ``wire_2D_tphi`` without it, or with the default settings, stops at its third step on the iteration limit).
  **But not** when a step starts from a current-free conductor in a changing field (the first transient step after a
  field-cooled static step, example ``bulk_2D_tphi_axisym``): the tangent at J = 0 is the floor :math:`\rho_0`, the
  first full step overshoots (the residual grows 1e4-1e5 times) and a line search that must reduce the residual stalls
  on it ("Unrecoverable failure to converge"); the full steps converge (15-18 iterations). Try the default (False)
  first; switch it on when the factorization fails or the iterations run away.
* **End the E-J table with the normal state**, a last segment of the normal resistivity above some :math:`E_\text{max}`
  (current sharing with the stabilizer; ``power_law_table(..., e_max=, rho_n=)`` of the benchmarks): the curve is
  extended along its last segment, so a pure power law extrapolated far above :math:`J_c` gives tangents of
  :math:`10^{30}` and more where an iterate overshoots.
* **Keep the floor negligible**: its own diffusion time :math:`\mu_0 L^2/\rho_0` must be much longer than the
  simulated time, else the field soaks in through it (``bean_slab`` with :math:`\rho_0=10^{-14}` over a 2e5 s cycle:
  a resistor, not a superconductor). In T-phi :math:`\rho_0` can be a numerical zero (1e-18): the mass term keeps the
  system regular.
* Newton with the Brent line search (default) and plain backtracking give the same results; Picard (``"picard":
  True``) does not converge on the power law (100 iterations). Iterations per step grow with n: at most 8 (n = 20)
  to 36 (n = 100) on ``powerlaw_hts``; halving the time step lowers the maximum (36 to 21 at n = 100).
* The power law tends to Bean's critical state as n grows: the AC loss is 1.22 times the critical-state value at
  n = 30 on the slab and on the wire alike, 1.09 (slab) and 1.11 (wire) at n = 100, the same with half the time step:
  compare with the critical state at the n of the material.

Selecting a B-H curve
---------------------

.. code-block:: python

   solver.model.set_bh_curve(26, [        # rows [H (A/m), B (T)], B strictly increasing
       [0.0, 0.0],
       [100.0, 0.37],
       [1113.8, 1.76],
       [20396.8, 2.43],
   ])
   # same as
   solver.model.set_nonlinear_law(26, "bh_curve", B=[0.0, 0.37, 1.76, 2.43], H=[0.0, 100.0, 1113.8, 20396.8])

In ``model_data.json``:

.. code-block:: json

   "nonlinear_laws": {"26": {"type": "bh_curve", "B": [0.0, 0.37, 1.76, 2.43], "H": [0.0, 100.0, 1113.8, 20396.8]}}

An id may have a list of entries for different laws (``[{"type": ...}, {"type": ...}]``); the same law twice on an id
is an error.

The curve
~~~~~~~~~

* ``B`` and ``H`` have the same length (at least 2), finite values, ``B`` strictly increasing; start at (0, 0).
* The curve is **piecewise linear** between the points and is extended beyond the last point along its last segment
  (and below the first along the first). Give the table up to deep saturation, and make the last segment have the
  slope you want there (close to :math:`\mu_0` in B per H, i.e. :math:`dH/dB\approx1/\mu_0`): adding a far point such
  as ``[2.34e6, 5.0]`` does that.
* The initial reluctivity (:math:`B\to0`) is the slope of the first segment.
* Use a monotone table (H increasing with B): a decreasing segment makes the tangent indefinite and Newton may fail.
* The ``mur`` of the material is still required by validation but is not used where the law applies.

The law is isotropic (it depends on :math:`|B|` only) and stateless (no hysteresis). A nonlinear material can be in any
domain: soft iron (``SOFT_MAGNETIC``), a magnet, a conducting steel (``SOLID_CONDUCTOR``, eddy currents and
saturation together), even a coil.

Newton iterations
-----------------

Each step iterates from the solution of the previous step. The settings are in ``solver_settings["nonlinear"]``
(:doc:`solver_settings`): residual and step tolerances, maximum iterations, line search. The log shows, for each step:

* the Newton history: residual norm, its reduction, the step norm, the relative step and the line-search step length;
* the **operating range** of each law: the largest element mean of :math:`|B|` in each nonlinear material, with a
  warning when it is beyond the last point of the curve (the curve is then extrapolated);
* the outcome (converged, or the reason it stopped).

``solver.steps[-1]["report"]`` holds the same data in Python, and the steps CSV the iteration counts.

Convergence tips:

* start with a static step at a small load, or ramp the load over a few steps;
* smaller time steps help in transients that cross the knee of the curve;
* if the line search fails, try ``"brent_line_search": False`` (backtracking), or a smaller ``minimum_step_length``;
* ``"continue_after_max_iterations": True`` (default) accepts the last iterate after ``max_iterations``: check the
  convergence column of the steps table; set it to False to stop the run instead.

Time Steps and Analysis Types
=============================

An analysis is a sequence of steps. Each ``solver.solve()`` solves one step, defined by ``solver.time_step``:

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Key
     - Default
     - Meaning
   * - ``time``
     - 1.0
     - End time of the step [s], > 0. The step length is ``dt = time - previous time`` (the first step starts at 0).
   * - ``transient_effects``
     - True
     - False (or 0): **static** step, no damping (no eddy currents); True: transient step.
   * - ``theta``
     - 1.0
     - Theta method, in (0, 1]: 1 backward Euler, 0.5 Crank-Nicolson (first-order formulations: all the element types
       of the library).
   * - ``beta``, ``gamma``
     - 0.25, 0.5
     - Newmark parameters of second-order formulations, in (0, 0.5] and [0, 1].

Other keys are ignored. ``time_step`` is read at every solve, so it can change from step to step.

Static analysis
---------------

.. code-block:: python

   solver.time_step["transient_effects"] = False
   solver.time_step["time"] = 1.0       # any positive value
   solver.solve()

Magnetostatics with the imposed currents (``Js``, coil currents): no eddy currents. Time-integrated electric unknowns
(``_timint`` types: V of the conductors, VOLT and EMF of the coils) are fixed to zero on a static step.

Transient analysis
------------------

.. code-block:: python

   solver.time_step["theta"] = 1.0
   for n, t in enumerate(times):
       solver.time_step["time"] = t
       solver.time_step["transient_effects"] = n > 0      # a static first step: start from equilibrium
       set_loads(t)
       solver.solve()
       solver.write_output()

* The fields start at zero. A first transient step with a non-zero load is a **step change** of the load at
  :math:`t=0`; to start from the magnetostatic state of the initial loads, make the first step static (at a small
  time, e.g. ``1e-10``).
* Choosing theta: 1 (backward Euler) is robust and damps everything, first-order accurate: small time steps;
  0.5 (Crank-Nicolson) is second-order accurate but does not damp: after a sudden change of a load the eddy currents
  can oscillate from step to step. Values such as 0.6 - 0.75 are a compromise.
  The current densities written in ``Jt`` use the mean rate over the step, which matches the solution only for
  theta = 1 (:doc:`outputs`): prefer theta = 1 when the eddy currents are post-processed after sudden changes.
* The time step should resolve the shortest time constant of interest, :math:`\tau\approx\mu_0\sigma\ell^2` for a
  conductor of size :math:`\ell` (L/R of a circuit).
* Repeating a solve at the same time re-solves the step (e.g. with other loads, as the plasma algorithm does at
  every Picard iteration); theta, beta, gamma and transient_effects must then be unchanged. A time smaller than the
  previous one is an error (the solver stops).

Matrix reuse
------------

The element matrices are computed at the first solve of a session; the system matrix and its factorization are
reused at the following steps as long as **dt, theta (beta, gamma), transient_effects and the set of constrained DOFs**
do not change. The step report says whether the matrix was "assembled and factorized" or "reused". For long
transients keep a constant dt where possible; a change of dt costs one factorization.

Steps and their records
-----------------------

After each solve, ``solver.steps`` holds a record: ``index``, ``time``, ``nonlinear``, ``newton_iterations``,
``linear_iterations``, ``convergence``, ``true_residual``, ``max_law_input``, ``wall_seconds``, ``output`` (the file
written after it) and ``report`` (the full C++ step report, a dict). When the solver is closed the records are written
to ``output/<name>_steps.csv`` and summarized in the log.

Second-order formulations
-------------------------

The solver supports formulations of second order in time (``time_order = 2``, a mass matrix M) with the Newmark
scheme. The element types of the current library are all first order; ``beta`` and ``gamma`` are used only by an
element type that declares ``time_order: 2`` (:ref:`theory-time-integration`).

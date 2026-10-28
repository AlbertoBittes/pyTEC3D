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
* Choosing theta: 1 after sudden changes, a lower value only for smooth excitations (:ref:`choosing-theta`).
* The time step should resolve the shortest time constant of interest, :math:`\tau\approx\mu_0\sigma\ell^2` for a
  conductor of size :math:`\ell` (L/R of a circuit).
* Repeating a solve at the same time re-solves the step (e.g. with other loads, as the plasma algorithm does at
  every Picard iteration); theta, beta, gamma and transient_effects must then be unchanged. A time smaller than the
  previous one is an error (the solver stops).

.. _choosing-theta:

Choosing theta
--------------

A component of the solution that decays with time constant :math:`\tau` (an eddy current mode, the L/R of a circuit)
is multiplied at each step by

.. math::

   R = \frac{1-(1-\theta)\,\Delta t/\tau}{1+\theta\,\Delta t/\tau}
   \;\xrightarrow{\ \Delta t\gg\tau\ }\; -\frac{1-\theta}{\theta}.

For a time step much longer than :math:`\tau` (fast modes: thin skin layers, small conductors, the first steps after a
jump), theta = 1 gives :math:`R\to 0`: the mode is damped at once. Theta = 0.5 gives :math:`R\to -1`: the mode
is not damped and changes sign at every step. The fields, and the current densities computed from them, then oscillate
from step to step (``Jt`` reports the rates of the equations, so it shows the oscillation as it is,
see :doc:`outputs`).

A sudden change also harms theta < 1 when :math:`\Delta t\ll\tau`. The step update
:math:`\mathbf{u}^{n+1}=\mathbf{u}^n+\Delta t\,[(1-\theta)\mathbf{v}^n+\theta\mathbf{v}^{n+1}]` assumes that the rate
varies linearly over the step, so a fraction :math:`1-\theta` of the rate before the change still acts after it. In
``coil_2D_avolt`` (100 A imposed at the first step, then a voltage drive of 1 mV) the driven ring keeps about 100 A with
theta = 1, as it should, but settles near 200 A with theta = 0.5: half of the voltage that built up the 100 A at step 1
is applied again at step 2.

**Use theta = 1 (backward Euler) after any sudden change:**

* a step change of a load (``Js``, ``AMPS``, ``CURR``), including the first transient step from zero fields with a
  non-zero load;
* a switch of excitation (current driven to voltage driven, a Dirichlet value or a coupling added or removed);
* a time step much longer than the time constants of some conductors of the model, which then react as to a jump;
* fast transients such as a plasma current quench.

It is first-order accurate: resolve the transient with time steps short against the time constants of interest.

**Use a lower theta only for smooth excitations, with time steps that resolve the time constants:**

* theta = 0.5 (Crank-Nicolson, second-order accurate) for continuous loads (ramps, sinusoids, slowly varying
  plasma currents) when :math:`\Delta t\lesssim\tau` for every conductor that matters: the same accuracy with
  fewer steps than theta = 1;
* 0.6 - 0.75 as a compromise where some fast modes remain (:math:`R\to -0.67` and :math:`-0.33`): some damping, in
  practice more accurate than theta = 1, but formally first order.

Theta can be changed between steps: a common recipe is one or two steps with theta = 1 right after each sudden change,
then back to 0.5. A change of theta costs one refactorization of the system matrix (see below).

Matrix reuse
------------

The element matrices are computed at the first solve of a session; the system matrix and its factorization are
reused at the following steps as long as **dt, theta (beta, gamma), transient_effects and the constraint rows** (the
constrained DOFs, their masters and coefficients; not the prescribed values) do not change. The step report says
whether the matrix was "assembled and factorized" or "reused". For long transients keep a constant dt where possible;
a change of dt costs one factorization.

A factorization has a symbolic part (ordering and analysis of the sparsity pattern) and a numeric one. The solve
keeps the reduced system of the unconstrained DOFs between steps, so when only the values change (a new dt, theta or
transient_effects, each Newton iteration) only the numeric factorization is redone; a full one only when the
constraint rows change. The "Factorizations" row of the step report (``full_factorizations``,
``numeric_factorizations`` of ``report``) counts them.

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

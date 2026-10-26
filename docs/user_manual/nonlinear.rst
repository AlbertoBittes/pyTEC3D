Nonlinear Materials
===================

Nonlinear laws
--------------

A nonlinear law replaces coefficients of the element types by values computed, at each quadrature point, from a field
of the solution and from the data of the law. The laws are defined in ``element_library/pytec3d_element_library/data/nonlinear_laws.json``
and compiled into the solver. The library has one law:

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

The **source** tells which id selects the law: ``mat`` a MATERIAL id, ``rc`` a CONSTANT id, ``me`` a TYPE id. A law
applies on the cells with that id, in every element type and domain that has a coefficient it replaces (for
``bh_curve``: all of them). A model with at least one such cell is solved with Newton's method at every step
(:ref:`theory-newton`).

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

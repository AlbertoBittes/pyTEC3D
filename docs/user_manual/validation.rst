Validation, Errors and Troubleshooting
======================================

Checks
------

pyTEC3D checks the input in three stages; each collects every problem before raising one
:class:`~pytec3d.errors.ValidationError` that lists them all.

#. **At load** (``load_model_into_solver()``, or ``solver.validate()``):

   * solver settings: names and ranges of every key (:doc:`solver_settings`);
   * mesh arrays: ``TYPE``, ``MATERIAL``, ``CONSTANT``, ``COORDINATE_SYSTEM``, ``node_angles`` (3 components) and the
     field arrays of the element types, with the right number of components and tuples;
   * model data against the registries compiled into the solver: element type and domain names / indices, option
     values, one element type and one axisymmetric flag for the whole model, cell dimension against element type,
     material and real constant names and shapes, coordinate systems, nonlinear laws and their curves, DOF coupling
     fields; ids used by the cells against the rows (missing ``TYPE``, ``MATERIAL`` and ``COORDINATE_SYSTEM`` rows are
     errors, missing ``CONSTANT`` rows warnings; unused rows warnings); every positive id of the ``CP_k`` arrays in
     ``dof_couplings``.

   The C++ solver then checks that every parameter an element needs in its domain is defined (e.g. ``rho`` of a
   conductor, ``turns`` of a coil) and stops with a message naming the TYPE, MATERIAL or CONSTANT id.

#. **Before each solve**: the time step (``time > 0``, ``theta`` in (0, 1], ``beta`` in (0, 0.5], ``gamma`` in
   [0, 1], ``transient_effects`` a bool), the arrays, the arrays that define the model in the session (``TYPE``,
   ``MATERIAL``, ``CONSTANT``, ``COORDINATE_SYSTEM``, ``node_angles``, ``CP_k``: unchanged since the load, otherwise
   ``ModelLockedError``), and two warnings: every load array zero ("No loads"), every Dirichlet array NaN ("No
   boundary conditions").

#. **During the solve**: linear convergence, the true residual (warning above 1e-6), MUMPS status and null pivots,
   Newton convergence, laws used beyond their curve.

Exceptions
----------

All derive from :class:`pytec3d.errors.PyTEC3DError` and from the matching built-in exception:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Exception
     - Raised when
   * - :class:`~pytec3d.errors.ValidationError` (ValueError)
     - The model, the settings or the time step are invalid (the message lists every problem).
   * - :class:`~pytec3d.errors.MissingArrayError` (ValueError)
     - An array is missing at the given location (often a missing ``location="elem"``).
   * - :class:`~pytec3d.errors.ModelLockedError` (RuntimeError)
     - The mesh or the model data are changed while the model is loaded.
   * - :class:`~pytec3d.errors.SolverStateError` (RuntimeError)
     - A call in the wrong state (solve before load, closed solver, no mesh).
   * - :class:`~pytec3d.errors.ConvergenceError` (RuntimeError)
     - The plasma algorithm did not converge.
   * - ``RuntimeError``
     - An error of the C++ solver (in a parallel run the solver aborts instead).

Troubleshooting
---------------

``MissingArrayError: ... 'Js' at 'node'``
   A cell array used with the default ``location="node"``: add ``location="elem"``.

``ValidationError: CP_k arrays: coupling ids not in data['dof_couplings']``
   CP ids painted in the mesh (possibly stored in a saved mesh) without a coupling: add them to ``dof_couplings``, or
   reset the array (``mesh.remove_array`` + ``mesh.add_array``) before building new couplings.

``element_types: all the cells must have the same element_type``
   A model uses one formulation: change the ``element_type`` of the rows, keep different domains.

``'<name>' is not defined for element type ... in domain ...`` (C++)
   A material property or real constant needed by a cell is missing from its row.

``Can only use Cholesky decomposition with symmetric positive definite matrices`` (C++)
   The mass matrix of the projection of the nodal fields is singular: in 2D, the out-of-plane length is zero (a
   ``thickness`` of 0 in a constants row). Leave ``thickness`` out (default 1) or give the real value.

``NEDELEC_ONE on geometric element HEX8 only supports ...`` (libMesh)
   Edge element types need a quadratic mesh.

The results are zero
   Check the warnings "No loads" / "No boundary conditions", the *Constraints* table (were the conditions applied to
   the right field?), and, for time-integrated unknowns, that the step is not static.

The solution is noisy or the true residual is large
   A singular system: a floating potential V, a missing gauge (edge elements without tree gauge), a 2D model without
   any A condition. Enable ``"-mat_mumps_icntl_24": "1"`` and read the number of null pivots in the step report.

Eddy currents oscillate from step to step
   Crank-Nicolson (theta = 0.5) after a sudden load change: use a larger theta or smaller steps.

Newton does not converge
   See :doc:`nonlinear` (load ramps, smaller steps, line search settings, the curve beyond saturation).

Out of memory in MUMPS
   ``"-mat_mumps_icntl_14": "50"`` (more workspace), more MPI ranks, or an iterative solver.

Debug output
   ``pytec3d.set_log_level("DEBUG")`` or ``PYTEC3D_LOG_LEVEL=DEBUG``.

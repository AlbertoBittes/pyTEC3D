Solver Settings
===============

``solver.solver_settings`` has three groups. Keys left out keep their default
(:data:`pytec3d.solver.DEFAULT_SOLVER_SETTINGS`); unknown keys are errors, so a misspelled setting is never ignored
silently. The settings are read at every ``solve()``, so they can change between two steps: a change of the
``linear`` group or of the PETSc options rebuilds and refactorizes the system matrix at the next solve, and a PETSc
option removed (set to ``None``) is cleared.

The linear systems are solved with `PETSc <https://petsc.org/>`_ [PETScManual]_: the ``linear`` group configures a
KSP (Krylov solver) and its PC (preconditioner), described in the `KSP chapter
<https://petsc.org/release/manual/ksp/>`_ of the PETSc manual; the ``nonlinear`` group configures the Newton solver of
libMesh (:libmesh:`NewtonSolver`) [Kirk2006]_.

.. code-block:: python

   solver.solver_settings = {
       "linear": {                       # PETSc KSP, alone or inside each Newton iteration
           "ksp_type": "preonly",        # a supported PETSc KSP name
           "pc_type": "lu",              # a supported PETSc PC name
           "relative_tolerance": 1e-6,   # finite, > 0
           "absolute_tolerance": 1e-10,  # finite, > 0, <= relative_tolerance
           "divergence_tolerance": None, # None (PETSc default) or > 0
           "max_iterations": 100,        # integer >= 1
           "initial_guess_nonzero": False,
           "factor_solver_type": "mumps",  # package of lu / cholesky; None for the PETSc default
       },
       "nonlinear": {                    # libMesh Newton solver (models with nonlinear laws)
           "max_iterations": 50,
           "absolute_residual_tolerance": 1e-35,
           "relative_residual_tolerance": 1e-8,
           "absolute_step_tolerance": 1e-8,
           "relative_step_tolerance": 1e-8,
           "require_residual_reduction": False,
           "require_finite_residual": True,
           "brent_line_search": True,    # False: plain backtracking
           "minimum_step_length": 1e-5,  # in (0, 1]
           "linear_tolerance_multiplier": 1e-3,
           "continue_after_max_iterations": True,
           "continue_after_backtrack_failure": False,
           "verbose": True,
       },
       "petsc_options": {                # set before each linear solve; "" for a flag, None removes a default
           "-ksp_monitor": "", "-ksp_view": "", "-ksp_converged_reason": "", "-log_view": "",
           "-ksp_monitor_true_residual": "", "-ksp_monitor_solution": "", "-mat_mumps_icntl_4": "2",
       },
   }

Changing one setting only:

.. code-block:: python

   solver.solver_settings["nonlinear"] = {"max_iterations": 20}            # the other nonlinear keys keep their default
   solver.solver_settings["petsc_options"] = {"-ksp_view": None, "-log_view": None}   # remove two defaults

``solver.effective_solver_settings()`` returns the merged settings; the log prints the *Solver configuration* table
and lists the settings that differ from the defaults.

Linear solver
-------------

The default is a **direct** solve: ``preonly`` + ``lu`` (:petsc:`KSPPREONLY <KSP/KSPPREONLY>`, :petsc:`PCLU
<PC/PCLU>`) with **MUMPS** (`MUMPS <https://mumps-solver.org/>`_, the parallel multifrontal solver [Amestoy2001]_,
through :petsc:`MATSOLVERMUMPS <Mat/MATSOLVERMUMPS>`), robust for the indefinite and ill-conditioned systems of
eddy-current problems (curl-curl operators, coupled electric unknowns, edge elements), in serial and in parallel. Its
memory grows quickly with the size of 3D models.

Iterative solvers (``ksp_type`` ``gmres`` [SaadSchultz1986]_, ``fgmres``, ``bcgs``... with ``pc_type`` ``ilu``,
``asm``, ``bjacobi``, ``hypre``, ``gamg``...) use less memory but need a preconditioner suited to the problem; generic
preconditioners converge poorly on curl-curl operators, which need auxiliary-space preconditioners such as the AMS of
`hypre <https://hypre.readthedocs.io/>`_ [HiptmairXu2007]_ (:petsc:`PCHYPRE <PC/PCHYPRE>`). Check
``linear_iterations`` and the true residual of each step; the tolerances are those of :petsc:`KSPSetTolerances
<KSP/KSPSetTolerances>`.

Accepted ``ksp_type`` (:petsc:`KSPType <KSP/KSPType>`): ``preonly gmres fgmres bcgs bcgsl cgs cg minres symmlq lgmres
ibcgs gcr tfqmr tcqmr cr lsqr bicg richardson chebyshev``. Accepted ``pc_type`` (:petsc:`PCType <PC/PCType>`): ``lu
cholesky ilu icc asm gasm bjacobi jacobi sor eisenstat hypre ml gamg mg ksp none``.

After every linear solve the solver computes the **true residual** :math:`\|b-Ax\|/\|b\|` of the reduced system and
warns above :math:`10^{-6}`.

Useful PETSc / MUMPS options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every ``-mat_mumps_icntl_<n>`` option of PETSc is listed in :petsc:`MATSOLVERMUMPS <Mat/MATSOLVERMUMPS>`, the meaning
of each ICNTL parameter in the `MUMPS user guide <https://mumps-solver.org/index.php?page=doc>`_. The monitors are
described in :petsc:`KSPMonitor <KSP/KSPMonitor>`, :petsc:`KSPView <KSP/KSPView>` and the `profiling chapter
<https://petsc.org/release/manual/profiling/>`_ (``-log_view``).

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Option
     - Effect
   * - ``"-mat_mumps_icntl_4": "2"``
     - MUMPS verbosity (default 2; ``"0"`` for silence).
   * - ``"-mat_mumps_icntl_24": "1"``
     - Null pivot detection: the step report then gives the number of null pivots; non-zero means a singular matrix
       (a missing boundary condition, a floating potential).
   * - ``"-mat_mumps_icntl_14": "50"``
     - Extra workspace (percent) when MUMPS stops with a memory error (INFOG(1) = -9).
   * - ``"-mat_mumps_icntl_7": "5"``
     - Ordering (5 = METIS [KarypisKumar1998]_).
   * - ``"-ksp_monitor": None``, ``"-ksp_view": None``, ``"-log_view": None``
     - Remove the default monitors (they print to the standard output only, not to the log).

Nonlinear solver
----------------

Used when at least one cell has a nonlinear law (:doc:`nonlinear`, :ref:`theory-newton`).

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Setting
     - Meaning
   * - ``max_iterations``
     - Maximum Newton iterations per step.
   * - ``absolute_residual_tolerance`` / ``relative_residual_tolerance``
     - Converged when the residual norm, or its ratio to the initial one, is below.
   * - ``absolute_step_tolerance`` / ``relative_step_tolerance``
     - Converged when the update norm, or its ratio to the solution norm, is below.
   * - ``require_residual_reduction``
     - Reject iterates that do not reduce the residual (with the line search).
   * - ``require_finite_residual``
     - Stop on NaN / Inf residuals.
   * - ``brent_line_search``
     - Line search by Brent's method [Brent1973]_ (True) or plain backtracking (False).
   * - ``minimum_step_length``
     - Smallest step of the line search, in (0, 1].
   * - ``linear_tolerance_multiplier``
     - Tolerance of each linear solve, relative to the current residual (iterative solvers).
   * - ``continue_after_max_iterations``
     - Accept the last iterate when ``max_iterations`` is reached (True) or fail.
   * - ``continue_after_backtrack_failure``
     - Accept the iterate when the line search fails.
   * - ``verbose``
     - Print the Newton iterations (standard output).

Threads and MPI
---------------

``Solver(num_threads=n)`` sets the threads per MPI process of the C++ solver (the first session of a process sets it
for the process). Parallel runs use ``mpirun -np <n> python script.py``; the memory used by the direct solver is reported
in the step report of every factorization (total and on the largest rank).

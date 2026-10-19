Results and Outputs
===================

The solution in the mesh arrays
-------------------------------

After each ``solve()`` the solution is written back into the arrays of ``solver.model.mesh``, so it can be read in
Python right away:

.. list-table::
   :header-rows: 1
   :widths: 14 12 74

   * - Array
     - Where
     - Content
   * - ``A``
     - nodes
     - Vector potential [Wb/m], global components (2D: component 2). With edge elements: the edge field evaluated at
       the nodes (3 components at corner and midside nodes).
   * - ``V``
     - nodes
     - Electric potential of the 3D conductors [V] (time integral [V s] where time-integrated).
   * - ``VOLT``, ``EMF``
     - nodes
     - Terminal and induced voltage of the conductors and coils (time integrals in the ``_timint`` types).
   * - ``B``, ``H``
     - nodes
     - Flux density [T] and field [A/m], :math:`L^2`-projected per cell and averaged at the nodes.
   * - ``B``, ``H``
     - cells
     - The same, mean over the quadrature points of the cell.
   * - ``Jt``
     - cells
     - Total current density [A/m^2]: imposed (``Js``), eddy / conduction (:math:`-\sigma(\partial_t\mathbf{A}+\nabla
       V)`) and coil current density, depending on the domain (:doc:`element_library`).

How these fields are computed (projection, quadrature mean, time derivative as the mean rate over the step) is
explained in :ref:`theory-postprocessing`.

.. warning::

   Two limits of the derived fields, to keep in mind when post-processing:

   * **Element means.** The element values are the plain mean of the quadrature-point values. In a conductor whose
     current flows in a skin layer thinner than the cells (short time steps, coarse mesh) the sum of ``Jt`` times
     the cell area over a cross-section can be far from the total current (e.g. 155 A for 100 A imposed in the first
     tutorial). Refine the mesh, or integrate the nodal solution, for integral quantities.
   * **theta < 1.** The time derivatives in the derived fields are the mean rates over the step,
     :math:`(\mathbf{u}^{n+1}-\mathbf{u}^n)/\Delta t`, while the equations of the theta method use the theta
     velocity. They coincide for theta = 1; for theta < 1, after a sudden change of a load or of a constraint, the eddy
     and conduction current densities in ``Jt`` can differ widely from the solution (with Crank-Nicolson they may
     alternate in sign from step to step). Use theta = 1, or theta close to 1, when currents are post-processed.

   The part of a field that does not depend on the solution (``Jt = Js`` in a ``STRANDED_CONDUCTOR``) is evaluated
   at every output from the current loads: ``Jt`` follows the changes of ``Js`` between steps.

Reading results in Python:

.. code-block:: python

   from pytec3d import mesh as meshes, selection

   grid = solver.model.mesh
   b = meshes.array_values(grid, "B", "node")                      # (points, 3)
   jt = meshes.array_values(grid, "Jt", "elem")                    # (cells, 3)
   area = meshes.array_values(grid, "cell_volume", "elem")         # cell areas (2D) / volumes (3D)

   # Net current through the cross-section of a 2D conductor (approximate: mean J times the cell area)
   selection.select(grid, "MATERIAL", 24, location="elem")
   ids = selection.selected_ids(grid, "elem")
   current = float((jt[ids, 2] * area[ids]).sum())

Files
-----

``write_output()``
   Saves ``output/<name>_out_<n>.vtu`` (rank 0): the whole model (mesh, every array, the model data). ``n`` is the
   number of steps solved by this solver so far (writing again after the same step overwrites the file). On the first
   call of a run, the outputs of a previous run of the same model (``<name>_out_*.vtu``, ``<name>_steps.csv``) in
   ``output`` are deleted, so the series only holds the current run. The plasma algorithm solves many times per step
   (one solve per Picard iteration): ``n`` then counts those solves.

``<name>.log``
   The job log (``Solver(log_file=...)``), started when the model is read and overwritten by a new run. It holds every
   Python and C++ message and the report tables:

   * **Run**: date, command line, working directory, host, platform, Python / numpy / pytec3d versions, and the
     SHA-256 of the input files (mesh, model data, script);
   * **Solver banner**: version, build mode and stamp, git commit (``-dirty`` with uncommitted changes), compiler,
     libMesh and PETSc versions, MPI ranks, threads;
   * **Mesh** and **Mesh quality** (per cell shape: largest edge ratio, smallest and largest angle, smallest scaled
     Jacobian with the worst cell, distorted and degenerate cells, cells reoriented on import);
   * **Model**: the analysis (2D planar, 2D axisymmetric, 3D), units, the element types of each TYPE id (domain, order,
     out-of-plane length, cells, laws), materials, real constants, nonlinear laws, DOF couplings;
   * **Solver configuration** and **Systems** (variables, FE families and orders, DOFs per rank);
   * for each step: **Time Step** (scheme and coefficients), **Constraints** (static-step, Dirichlet and coupling
     constraints by variable) and the **step report**: kind, DOFs, whether the element matrices and the system matrix
     were computed or reused, linear iterations and residuals, the true residual, the MUMPS statistics (status,
     memory, entries, fill-in, operations, delayed and null pivots), the Newton iterations and the operating range of
     the nonlinear laws, the wall time of each phase;
   * **Run summary** at ``close()``: one row per step and the totals (wall and CPU time, warnings, errors).

   The raw output of PETSc and libMesh (``-ksp_monitor``, ``-log_view``, MUMPS, the Newton lines) goes to the standard
   output only.

``output/<name>_steps.csv``
   One row per step: index, time, nonlinear, Newton and linear iterations, convergence, true residual, largest law
   input, wall time, output file.

``emat/``
   The element matrices of the session (binary, one file per rank), reused between steps. They can be deleted after
   the run.

Post-processing in ParaView
---------------------------

Open the ``output/<name>_out_*.vtu`` series in ParaView (it is recognised as a time series). Useful filters: *Cell
Data to Point Data*, *Glyph* (B, Jt as arrows), *Stream Tracer* (field lines), *Slice*, *Threshold* on ``TYPE`` or
``MATERIAL`` (show a region), *Plot Over Line*, *Integrate Variables* (integrals of ``Jt`` over a slice give
currents). In 2D axisymmetric models remember that the integration measure is :math:`2\pi r\,dr\,dz`. The pyTEC3D
plugin (:doc:`paraview`) adds the selection filters.

Derived quantities in Python
----------------------------

:mod:`pytec3d.model_edit` provides a few helpers on the mesh: ``element_solution`` (point array averaged on the cells,
as ``<name>e``), ``poloidal_flux`` (the flux function :math:`\psi = rA_\theta` of an axisymmetric field),
``add_component`` (one component of an array as a new array), ``to_local_coordinates`` and
``map_vector_field_to_local_components`` (e.g. B in cylindrical components).

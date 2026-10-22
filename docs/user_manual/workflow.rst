The Analysis Workflow
=====================

The objects
-----------

:class:`pytec3d.solver.Solver`
   The analysis. It owns a :class:`~pytec3d.model.Model`, the time step settings (``solver.time_step``), the solver settings
   (``solver.solver_settings``), the C++ solver session and the run log. Create one per analysis, and close it when
   done (``solver.close()``, or use it as a context manager: ``with Solver() as solver: ...``). Arguments:
   ``num_threads`` (threads per MPI process), ``root_only_output`` (silence the standard output of the MPI ranks
   other than 0, default True), ``log_file`` (True: ``<mesh_name>.log`` in the working directory; a path; or False).

:class:`pytec3d.model.Model`
   The mesh (``model.mesh``: a ``vtkUnstructuredGrid``) and the model data (``model.data``: a dict with the sections
   ``element_types``, ``materials``, ``constants``, ``dof_couplings``, ``nonlinear_laws``, ``coordinate_systems``).
   The model is named after the mesh file (``model.mesh_name``); output files use that name.

:mod:`pytec3d.selection`
   Selection of nodes and cells and assignment of values to arrays: the way loads, boundary conditions, couplings,
   element types and materials are applied to regions of the mesh (:doc:`selection`).

:class:`pytec3d.algorithm.Algorithm`
   The plasma equilibrium and evolution algorithm, which drives a ``Solver`` (:doc:`plasma_algorithm`).

The life cycle of an analysis
-----------------------------

.. code-block:: python

   from pytec3d import Solver, selection

   solver = Solver()

   # 1. Model: mesh + model data
   solver.read_model("model.vtu")            # adds the arrays of the element types, resets them, clears the selection
   solver.model.read_json("model_data.json")

   # 2. Preparation: everything that changes the model (TYPE / MATERIAL of regions, couplings, coordinate systems,
   #    nonlinear laws) must be done now ...
   selection.select(solver.model.mesh, "MATERIAL", 2, location="elem")
   ...

   # 3. Load the model into the solver: validation, mesh distribution, DOF numbering
   solver.load_model_into_solver()

   # 4. Steps: loads and Dirichlet values may change between steps
   for t in times:
       solver.time_step["time"] = t
       ...                                   # set loads, Dirichlet values
       solver.solve()
       solver.write_output()

   solver.close()                            # run summary, steps CSV, closes the log

Reading the model
   ``Solver.read_model(path)`` reads a ``.vtu`` (or ``.pvtu``) file. With the default ``reset_arrays=True`` the
   arrays of the element types (solution, Dirichlet values, loads) are set back to their defaults (0, or NaN for the
   Dirichlet arrays), the selection is cleared and the midside-node mask ``MidsideNodes`` is computed: a mesh saved
   by a previous run starts clean. Pass ``reset_arrays=False`` to keep the values stored in the file (for instance to
   continue from a saved, prepared model). If the file holds model data (a model saved by pyTEC3D), it is read too;
   ``Model.read_json()`` replaces it.

Preparing the model
   Before the model is loaded, everything can change: the cell arrays ``TYPE``, ``MATERIAL``, ``CONSTANT``,
   ``COORDINATE_SYSTEM``, the model data (``model.data``, ``set_bh_curve()``, ``set_coordinate_system()``,
   ``periodic_boundaries()``...), and the arrays. ``model.update_grid()`` (called by ``read_json()`` and by the load)
   adds the arrays needed by the element types named in the model data and removes those of the other element types
   of the library.

Loading the model
   ``Solver.load_model_into_solver()`` checks the solver settings and the model (:doc:`validation`), prints the
   provenance table (input files and their SHA-256), and creates the C++ session: the mesh is read from the VTK grid,
   distributed over the MPI ranks, the element types, materials and real constants are resolved per cell, the DOFs
   are numbered. The session **copies the model data and the solver settings**: later changes of ``model.data`` or
   ``solver_settings`` only take effect at the next load. While loaded, the model is **locked**: the mesh and the
   model data cannot change (``model.data`` becomes a read-only copy: any edit raises
   :class:`~pytec3d.errors.ModelLockedError`); the values of the mesh arrays can, except those that define the model in
   the session (table below).

Solving
   ``Solver.solve()`` solves the step that ends at ``time_step["time"]`` (:doc:`analysis`). Before each solve the time
   step and the arrays are checked, and warnings are logged when every load array is zero or every Dirichlet array
   is NaN. The C++ solver reads the current values of the loads and of the Dirichlet arrays, solves, and writes the
   solution back into the mesh arrays (``A``, ``V``, ``VOLT``, ``EMF``, ``B``, ``H``, ``Jt``...).

Writing the output
   ``Solver.write_output()`` saves the whole model (mesh, arrays, model data) as ``output/<name>_out_<n>.vtu``
   (:doc:`outputs`). ``Solver.save_model(path=None)`` saves it as ``<name>.vtu`` in the working directory, or at a
   given path: note that, run from the folder of the input mesh, it overwrites the input file.

Changing the model after a load
   Close the session (``solver.close_session()``), change the model, and load it again. A new session starts from
   zero fields.

What can change between steps
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Change
     - Allowed while loaded
     - Effect
   * - Values of loads (``Js``, ``CURR``, ``AMPS``...)
     - yes
     - Read at the next solve. The system matrix is reused.
   * - Values of Dirichlet arrays (``_A``, ``_V``, ``_VOLT``...)
     - yes
     - If the set of constrained DOFs changes (a NaN becomes a number or the opposite), the system matrix is rebuilt
       and refactorized; if only the values change, it is reused.
   * - ``time_step`` (time, theta, transient_effects)
     - yes
     - A change of dt, theta, beta, gamma or of transient_effects rebuilds the system matrix.
   * - ``TYPE``, ``MATERIAL``, ``CONSTANT``, ``COORDINATE_SYSTEM``, ``node_angles``, ``CP_k`` arrays
     - no
     - The element matrices, the element coordinate systems, the nodal frames and the couplings are built when the
       model is loaded: every ``solve()`` compares these arrays with their values at the load and stops with an error
       if one changed. Close the session,
       change, load again.
   * - ``model.data``
     - no
     - Read-only while loaded: an edit raises at once. Close the session, change, load again.
   * - solver settings
     - no
     - Read at the load: close the session, change, load again.

Files of a run
--------------

Run in a working directory, an analysis reads and writes:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File
     - Content
   * - ``<name>.vtu``, ``model_data.json``
     - The inputs: mesh and model data.
   * - ``<name>.log``
     - The job log: every Python and C++ message and all the report tables (:doc:`outputs`). Overwritten by a new
       run.
   * - ``output/<name>_out_<n>.vtu``
     - The model after step ``n`` (one per ``write_output()`` call).
   * - ``output/<name>_steps.csv``
     - One row per solved step, written when the solver is closed.
   * - ``output/<name>_iteration_diagnostics.csv``
     - Plasma algorithm only: one row per Picard iteration.

MPI parallel runs
-----------------

The same script runs in parallel with ``mpirun -np <n> python script.py``. Every rank executes the whole script and
holds the whole VTK grid (the Python preparation is replicated); the C++ solver distributes the finite-element mesh
and the system. Files are written by rank 0 only. Use ``solver.rank`` to restrict Python work (prints, plots) to one
rank.
